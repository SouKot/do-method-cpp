/*********************************************************************
 * Copyright by Jonas Thies, Univ. of Groningen 2006/7/8.             *
 * Permission to use, copy, modify, redistribute is granted           *
 * as long as this header remains intact.                             *
 * contact: jonas@math.rug.nl                                         *
 **********************************************************************/
#include "NOX_Config.h"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "Teuchos_RCP.hpp"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_LinearSystem_Hymls.hpp"
#include "NOX_Epetra_LinearSystem_Amesos.hpp"
#include "NOX.H"
#include "BelosTypes.hpp"
#include "Ifpack_Preconditioner.h"
#include "NOX_Epetra_LinearSystem_Belos.H"
//#include "TRIOS_SolverFactory.H"
#include "FVM_LocaInterface.H"
#include "FVM_Domain.H"
#include "FVM_model_interface.h"
#include "ImplicitTimeStepper.hpp"
#include "ThetaStepperEvaluator.H"
#include "globdefs.H"
#include "HYMLS_Tools.hpp"
#include "HYMLS_MatrixUtils.hpp"
#include "HYMLS_Preconditioner.hpp"
#include "HYMLS_Solver.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
using namespace std;
// this class needs a serious makeover or possibly should be replaced by Rythmos

///////////////////////////////////////////////
// Construct from OceanModel and NOX solver //
/////////////////////////////////////////////
ImplicitTimeStepper::ImplicitTimeStepper(
    Teuchos::RCP<FVM::LocaInterface> model_,
    Teuchos::ParameterList& stepperParams_, Teuchos::ParameterList& nlParams_)
: model(model_), paramList(stepperParams_), nlParams(nlParams_)
{
  DEBUG("ImplicitTimeStepper constructor done");
  // cout<<"\n"<<__FILE__<<__LINE__;
  out = paramList.get("Output Stream",info);
  predictor = paramList.get("Predictor","Constant");
  scheme = paramList.get("Scheme","Theta");
  if (scheme == "Theta")
  {
    theta = paramList.get("Theta",1.0); // backward Euler by default
  }
  else if (scheme!="Adaptive Theta")
  {
    Error("only 'Theta' and 'Adaptive Theta' Schemes are implemented!",__FILE__,__LINE__);
  }
  // initial step size
  double dt = paramList.get("Step Size",1.0);
  double t  = paramList.get("Start Time",0.0);
  haveNormEst = true; // currently we have only theta-schemes that support LTE estimates 

  //Create the Epetra_RowMatrix for the Jacobian/Preconditioner
  Teuchos::RCP<Epetra_CrsMatrix> A = model->getJacobian();

  thetaModel = Teuchos::rcp(new ThetaStepperEvaluator(model,model->getJacobian(),model->getSolution(),t,theta,dt));
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = thetaModel;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = thetaModel;

  if (iReq==Teuchos::null ||iJac==Teuchos::null)
  {
    Error("null pointer detected!",__FILE__,__LINE__);
  }

  // register pre/post operator (TODO: do we need it here?)
  Teuchos::RCP<NOX::Abstract::PrePostOperator> iPrePost = model;

  //Create the "Solver Options" sublist
  Teuchos::ParameterList& optParams = nlParams.sublist("Solver Options");

  // we need this so that the vmix_fix flag is handled correctly in THCM
  optParams.set("User Defined Pre/Post Operator",iPrePost);

  double TolNewton = nlParams.get("Convergence Tolerance",1.0e-6);
  Teuchos::ParameterList& searchParams=nlParams.sublist("Line Search");
  int    MaxIt  = searchParams.get("Max Iters",10);

  // Set up the Solver Convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormF(TolNewton));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters
    = Teuchos::rcp(new NOX::StatusTest::MaxIters(MaxIt));
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

  convTests = combo;
#ifdef DEBUGGING    
  convTests->print(*debug);
#endif

  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.set("MyPID", model->get_x_map()->Comm().MyPID());
  nlPrintParams.set("Output Stream",out);
  nlPrintParams.set("Error Stream",out);
  nlPrintParams.set("Output Process",0);
    //nlPrintParams.set("Output Information",
      //NOX::Utils::Details +
      //NOX::Utils::OuterIteration +
      //NOX::Utils::InnerIteration +
      //NOX::Utils::OuterIterationStatusTest +
      //NOX::Utils::LinearSolverDetails +
      //NOX::Utils::Debug +
      //NOX::Utils::Warning +
      //NOX::Utils::StepperDetails +
      //NOX::Utils::StepperIteration +
      //NOX::Utils::StepperParameters);


  uprime = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));
  ustart = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));
  //cout<<"\n length of ustart is: "<<ustart->MyLength()<<"\n";
  f_nm1 = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));
  f_n = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));
  f_np1 = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));

  DEBUG("Construct Linear System...");
  Teuchos::ParameterList& lsParams = nlParams.
    sublist("Direction").
    sublist("Newton").
    sublist("Linear Solver");

  // note that we want the jacobian/rhs to be evaluated by the ThetaStepperEvaluator,
  // not the original OceanModel
  // linSys  = THCM::Instance().createLinearSystem(model, lsParams, nlPrintParams, out,
  //               thetaModel, thetaModel);

  // get a pointer to the scaling arrays for the linear system
  Teuchos::RCP<NOX::Epetra::Scaling> scaling = model->getScaling();

  // the model serves as Preconditioner factory if "User Defined" is selected
  std::string PrecType = lsParams.get("Preconditioner","Ifpack");

  Teuchos::RCP<Epetra_Vector>soln= model->getSolution();
  if (PrecType == "User Defined")
  {
    DEBUG("user defined preconditioning");
    Teuchos::RCP<Epetra_Operator> myPrecOperator = model->getPreconditioner();
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = model;

    if (Teuchos::is_null(myPrecOperator))
    {
      HYMLS::Tools::Error("preconditioner is null!",__FILE__,__LINE__);
    }

    //      linsys = Teuchos::rcp(new NOX::Epetra::LinearSystemBelos(printParams,
    //                        lsParams, iJac, A, iPrec, myPrecOperator, soln, scaling));
/*     linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemHymls(nlPrintParams,
 * 	  lsParams, iJac, model->getJacobian(), iPrec, myPrecOperator,soln , scaling,
 * 	  model->getMassMatrix()));
 */
    /*LinearSystemAztecOO(nlPrintParams, lsParams, 
      iJac, jacobian, iPrec, preconditioner, cloneVector, s),
      massMatrix_(massMatrix)*/
  }
  else
  {
    DEBUG("Trilinos preconditioning");
    std::string solver_type = lsParams.get("Solver Type","Direct");
    if (solver_type == "Direct")
    {
      std::cout<<"\n===================================================\n";
      std::cout<<"USING AMESOS SOLVER FOR MEAN";
      std::cout<<"\n===================================================\n";
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAmesos(nlPrintParams,
	    lsParams, iReq, iJac, model->getJacobian(), soln, scaling));
    }
    else
    {
      std::cout<<"\n===================================================\n";
      std::cout<<"USING BELOS SOLVER FOR MEAN";
      std::cout<<"\n===================================================\n";
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemBelos(nlPrintParams, 
	    lsParams, iReq, iJac, model->getJacobian(), soln,scaling));
    }
  }


  NOX::Epetra::Vector noxSoln(*ustart);

  DEBUG("ImplicitTimeStepper: create NOX Group");
  // Create the Group
  curGroup = Teuchos::rcp(new NOX::Epetra::Group(nlPrintParams, 
	iReq, noxSoln, linSys));


  NOX::Solver::Factory nox_factory;
  DEBUG("ImplicitTimeStepper: create NOX Solver");
  nlSolver = nox_factory.buildSolver(curGroup, convTests, Teuchos::rcp(&nlParams,false));
  DEBUG("solver constructed successfully");

  DEBUG(" call reset()...");
  this->reset();

  (*out) << "Stepper Parameters: \n";
  (*out) << paramList;
  restarted = true;
  DEBUG("ImplicitTimeStepper constructor done");
}

////////////////////////////////////////////
// Destructor //////////////////////////////
////////////////////////////////////////////

ImplicitTimeStepper::~ImplicitTimeStepper()
{
}

//////////////////////////////////////////////////////////////
// reset stepper after failed step or at start of new run   //
//////////////////////////////////////////////////////////////
void ImplicitTimeStepper::reset()
{
  DEBUG("ImplicitTimeStepper: reset");
  restarted = true;
  lteEst = -1.0;
  effortEst = 1.0;
  dt_nm1 =  0.0;
  CHECK_ZERO(uprime->PutScalar(0.0));
  CHECK_ZERO(f_nm1->PutScalar(0.0));
  CHECK_ZERO(f_n->PutScalar(0.0));
  DEBUG("ImplicitTimeStepper: reset done");
}

//////////////////////////////////////////////////////////////

bool ImplicitTimeStepper::hasErrorEstimate(void) const
{
  return (haveNormEst && (lteEst!=-1.0));
}


///////////////////////////////////////////////////////////////////////////
// advance one step, Bu_np1 = Bu_n + dt*(theta*F(u_n+1)+(1-theta)F(u_n)) //
///////////////////////////////////////////////////////////////////////////

// after a step fails, the user has to call reset(). In that case we need to recompute
// f(u_n) before solving for u_{n+1}, and we are not able to provide an error stimate 
// for u_{n+1}.
bool ImplicitTimeStepper::Step(const Epetra_Vector& u_n, double t_n, 
    Epetra_Vector& u_np1, double dt)
{
  if (restarted)    // first step or restart after failed step:
  {               // need to evaluate f_n=F(u_n)
    thetaModel->reset(t_n,u_n);
    *f_n = thetaModel->get_f();
    // destroy preconditioner, regardless of what our strategy may be. This may mean
    // that the new precond is reused for a shorter time.
    linSys->destroyPreconditioner();
  }

  bool success = true; //indicates wether Newton converged or not

  double theta_old = theta;

  if (scheme=="Adaptive Theta")
  {
    double kappa; // kappa is used to select theta somewhere between 0.5 and 1.0, 
    // depending on the time step-size dt
    // TODO: check the correct scaling of dt, it should be 0<ds<=1 here!
    kappa = min(1.0,0.5/dt);
    theta = 0.5 + kappa*dt;
    (*out) << "Stepper: using theta = "<<theta<<std::endl;
  }

  if (theta!=theta_old)
  {
    thetaModel->set_theta(theta);
  }
  thetaModel->set_dt(dt);

  // predictor: we use the secant method or no predictor at all (constant)
  if (predictor=="Constant")
  {
    *out << "Constant Predictor..."<<std::endl;
    CHECK_ZERO(ustart->Update(1.0,u_n,0.0));
  }
  else if (predictor=="Secant")
  {
    *out << "Secant Predictor..."<<std::endl;
    CHECK_ZERO(ustart->Update(1.0,u_n,dt,*uprime,0.0));
  }
  else
  {
    Error("Invalid Predictor for implicit time-integration ",__FILE__,__LINE__);
  }

  //  nlSolver->reset(curGroup,convTests,nlParams);
  DEBUG("reset the solver...");
  NOX::Epetra::Vector noxvec(*ustart);
  nlSolver->reset(noxvec);

  // Compute next point on continuation curve
  DEBUG("solve using nox...");
  *out << "Newton Corrector..."<<std::endl;

  meanSolverTime->start();
  
  NOX::StatusTest::StatusType solverStatus = nlSolver->solve();  
  
  meanSolverTime->stop();
  meanSolverTime->incrementNumCalls();

  // Check solver status
  if (solverStatus == NOX::StatusTest::Failed) 
  {
    DEBUG("WARNING: Nonlinear Solve failed!");
    success=false;
    lteEst=-1;
    // skip the rest of the step. The caller has to call reset()
    // and try with a different step-size
    return success;
  }

  // Copy solution out of solver
  *curGroup = dynamic_cast<const NOX::Epetra::Group&>(nlSolver->getSolutionGroup());
  const NOX::Epetra::Vector& noxSoln = dynamic_cast<const NOX::Epetra::Vector&>(curGroup->getX());
  u_np1 = noxSoln.getEpetraVector();

  int numIters = nlSolver->getNumIterations();
  Teuchos::ParameterList& searchParams=nlParams.sublist("Line Search");
  int maxIters = searchParams.get("Max Iters",10);

  // the user can then ask how 'difficult' this step was:
  effortEst = (double)numIters/(double)maxIters;

  // we can normalize the pressure here, but this
  // is not the same as in the THCM implementation
  // TODO: do we need this???
  //  THCM::Instance().normalizePressure(u_np1);

  // compute secant to approximate u' for the next predictor
  if (predictor=="Secant")
  {
    uprime->Update(1.0/dt,u_np1,-1.0/dt,u_n,0.0);
  }

  // prepare evaluator for next step
  thetaModel->reset(t_n+dt,u_np1);
  // get f(u_new) from the theta evaluator
  *f_np1 = thetaModel->get_f();

  double nrm;
  f_np1->Norm2(&nrm);
  *out << "t = "<<t_n+dt<<": NORM OF F = "<<nrm<<std::endl;

  if (!restarted) // after a reset() we can't do this estimate
  {
    lteEst = estimateError(theta,
	*f_nm1,*f_n,*f_np1,
	dt_nm1,dt);
    //(*out) << "estimated truncation error: "<<lteEst<<std::endl;
  }

  dt_nm1 = dt;
  restarted = false;

  // avoid deleting the Teuchos::rcp f_nm1:
  Teuchos::RCP<Epetra_Vector> tmp = f_nm1;

  f_nm1 = f_n;
  f_n = f_np1;
  f_np1 = tmp;
  f_np1->PutScalar(0.0);

  curPrecAge++;

  return success;
}//Step

// this error stimate was taken from the THCM code (file time.f version 6.0)
// rh0/1/2 are the value of F at u_{n-1}, u_n and u_{n+1} respectively.
double ImplicitTimeStepper::estimateError(double theta,
    const Epetra_Vector& rh0,
    const Epetra_Vector& rh1, 
    const Epetra_Vector& rh2,
    double dt_old, double dt) const
{
  Epetra_Vector D(rh0.Map());

  DEBUG("estimating local truncation error...");
  DEBVAR(theta);
  DEBVAR(dt_old);
  DEBVAR(dt);

#if 0 //{
  // this version doesn't treat adaptive stepsizes, it is only
  //   retained for documentation purposes:

  // THCM version:
  //  D = kapp*(rh2-rh0)/(2*ds)-(rh2-2*rh1+rh0)/(12*ds**2)
  //         error = linrm(D*ds**3,ndim);

  Epetra_Vector D1(rh0.Map(),Copy);

  // approximate t-derivative of u at time level n (2nd order accurate)
  CHECK_ZERO(D1.Update(1.0/(2*ds),rh2,-1.0/(2*ds),rh0,0.0));

  // approximate tt-derivative of u at time level n (2nd order accurate)
  Epetra_Vector D2 = rh0;
  CHECK_ZERO(D2.Update(1/(ds*ds), rh2, -2/(ds*ds),rh1, 1.0));

#else // in the version below I have added proper finite differences for variable dt
  //}{
  Epetra_Vector D1 = rh1;

  double alpha = dt_old/dt;
  double denom = dt_old + alpha*alpha*dt;

  // approximate t-derivative of u at time level n (2nd order accurate)
  //      f_n-1 - (1-a^2) f_n - a^2 f_n+1
  // f' = --------------------------- + O(dt^2), a = (dt_n-1) / dt_n
  //             dt_n-1 + a^2 dt_n
  CHECK_ZERO(D1.Update(1.0/denom,rh0,-(alpha*alpha)/denom,rh2,(alpha*alpha-1.0)/denom));

  // approximate tt-derivative of u at time level n (2nd order accurate)
  //
  //         f_n-1 - (1+a) f_n + a f_n+1   
  // f'' = ------------------------------- 
  //       (1/2) ((dt_n-1)^2 + a (dt_n)^2) 
  //
  Epetra_Vector D2 = rh1;
  denom = dt_old*dt_old/2 + alpha*dt*dt/2;
  CHECK_ZERO(D2.Update(1.0/denom,rh0,alpha/denom,rh2,-(1.0+alpha)/denom));

#endif  //}
  // estimate LTE
  double c1 = dt*dt*(0.5-theta);
  double c2 = dt*dt*dt*(1.0/6.0 - theta/2.0);
  CHECK_ZERO(D.Update(c1,D1,c2,D2,0.0));

  //DEBVAR(D);

  double error;
  CHECK_ZERO(D.Norm1(&error));
  DEBVAR(error);
  return error;
}

