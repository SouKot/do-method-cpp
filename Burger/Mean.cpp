/*
 * =====================================================================================
 *
 *       Filename:  Mean.cpp
 *
 *    Description:
 *
 *    Calculates mean of Stochastic Burgers Equation
 *
 *        Version:  1.0
 *        Created:  04/14/2018 07:55:56 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#include "Mean.hpp"
#include "Amesos_Klu.h"
#include "Amesos_Scalapack.h"
#include "EpetraExt_RowMatrixOut.h"
#include "HYMLS_MatrixUtils.hpp"
#include <AztecOO.h>
#include <EpetraExt_MatrixMatrix.h>
//#include "Galeri_Maps.h"

// using namespace Galeri;
using namespace Teuchos;
const double PI = 3.14159265358979323846264338328;
/*
 * ===  FUNCTION
 * ======================================================================
 *         Name:  Mean(constructor)
 *  Description:
 * =====================================================================================
 */
Mean::Mean(RCP<Teuchos::ParameterList> PrmLst, double* t, double* dt,
           RCP<Epetra_Comm> Comm)
  : t_(t)
  , dt_(dt)
  , toleranceRHS_(10e-12)
  , maxNumIterations_(30)
  , numBackTrackingSteps_(5)
  , backTracking_(true)
{
  Teuchos::RCP<Epetra_Map> Map;

  m = PrmLst->get("nx", 200);  // global number of elements in the map
  mu = PrmLst->get("mu", 0.1); // global number of elements in the map
  theta = PrmLst->get("theta", 0.5);
  NumStchFrcVec_ = PrmLst->get("No. of vectors in stoch. forcing",
                               2); // global number of elements in the map
  solver_type = PrmLst->get("Solver Package", "Amesos2");
  test_ = PrmLst->get("Testing", true);
  debug_ = PrmLst->get("Debugging", true);
  if (debug_) {
    std::cout << "\n nx = " << m << std::endl;
    std::cout << "\n mu = " << mu << std::endl;
    std::cout << "\n numvecV = " << NumStchFrcVec_ << std::endl;
  }
  Map = rcp(new Epetra_Map(m, 0, *Comm)); // =CreateMap(MapType, *Comm, List);
  LinOp_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  Op2x_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 2));
  Op1x_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
  eye_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
  u_ = rcp(new Epetra_Vector(*Map));
  u0 = rcp(new Epetra_Vector(*u_));
  x_ = rcp(new Epetra_Vector(*u_));
  createLinOp(Comm);
  createW();
  LinOp_->FillComplete();
  Op2x_->FillComplete();
  Op1x_->FillComplete();
  eye_->FillComplete();
  MyPID = LinOp_->Comm().MyPID();
  Teuchos::ParameterList& SolverSublist = PrmLst->sublist(solver_type);
  std::string solver_name =
    SolverSublist.get("Solver name", "Amesos_Scalapack");
  const Teuchos::RCP<Teuchos::ParameterList> SolParams =
    Teuchos::getParametersFromXmlFile(solver_name + ".xml");
  if (MyPID == 0)
    std::cout << solver_type << " solver(" << solver_name
              << ") has been chosen as stoch basis solver";

  LinJac = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  NlinJac = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  ThetaJac = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  jac = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  SpDiag = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
  *SpDiag = *eye_;
  SpDiag->FillComplete();
  RHS = rcp(new Epetra_Vector(*u_));
  dx_ = rcp(new Epetra_Vector(*u_));
  Op1x_u = rcp(new Epetra_Vector(*RHS));
  Op2x_u = rcp(new Epetra_Vector(*RHS));
  rhs = rcp(new Epetra_Vector(*RHS));
  ExpVyVy_ = rcp(new Epetra_Vector(*RHS));
  tmp2 = rcp(new Epetra_Vector(*RHS));
  uu = rcp(new Epetra_Vector(*RHS));
  ThetaRHS = rcp(new Epetra_Vector(*RHS));

  TmpMat1 = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
  TmpMat2 = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
  // theta=0.5;

  Epetra_Vector rndmvct(*u_);

  EpetraExt::MatrixMatrix::Add(*LinOp_, false, 1.0, *LinJac, 0.0);
  LinJac->FillComplete();
  rndmvct = *u_;
  Op1x_->Multiply(false, rndmvct, *Op1x_u);
  Op2x_->Multiply(false, rndmvct, *Op2x_u);
  SpDiag->ReplaceDiagonalValues(*Op1x_u);
  EpetraExt::MatrixMatrix::Multiply(*SpDiag, false, *Op2x_, false, *TmpMat1,
                                    false);
  if (debug_) {
    std::cout << "\n frob. norm of spdiag(Op1 * u) = "
              << SpDiag->NormFrobenius() << std::endl;
  }
  SpDiag->ReplaceDiagonalValues(*u_);
  EpetraExt::MatrixMatrix::Multiply(*Op2x_, false, *SpDiag, false, *TmpMat2);
  if (debug_) {
    std::cout << "\n frob. norm of spdiag(Op2 * u) = "
              << SpDiag->NormFrobenius() << std::endl;
  }
  int err = EpetraExt::MatrixMatrix::Add(*TmpMat2, false, 1.0, *NlinJac, 0.0);
  if (debug_) {
    std::cout << "\n error in Tmpmat2+nlinjac: " << err << std::endl;
    std::cout << "\n if fillcommplete has been called for TmpMat1: "
              << TmpMat1->Filled() << std::endl;
  }
  if (TmpMat1->Filled() == false)
    TmpMat1->FillComplete();
  // err = EpetraExt::MatrixMatrix::Add(*TmpMat1, false, 1.0, *NlinJac, 1.0);
  if (debug_) {
    std::cout << "\n if fillcommplete has been called for NlinJac: "
              << NlinJac->Filled() << std::endl;
    std::cout << "\n error in Tmpmat1+nlinjac: " << err << std::endl;
  }
  NlinJac->FillComplete();
  //  Epetra_CrsMatrix *tmpjac, *tmpthetajac;
  // tmpjac=new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map,  1);
  //  tmpjac=jac.get();
  EpetraExt::MatrixMatrix::Add(*LinJac, false, 1.0, *jac, 0.0);
  EpetraExt::MatrixMatrix::Add(*NlinJac, false, -1.0, *jac, 1.0);
  // EpetraExt::MatrixMatrix::Add(*LinJac,false,1.0,*NlinJac,false,-1.0,tmpjac);
  jac->FillComplete();
  //  delete tmpjac;

  //  tmpthetajac=new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map,  1);
  //  tmpthetajac=ThetaJac.get();
  EpetraExt::MatrixMatrix::Add(*jac, false, *dt_, *ThetaJac, 0.0);
  EpetraExt::MatrixMatrix::Add(*eye_, false, 1.0, *ThetaJac, -1.0);
  //  EpetraExt::MatrixMatrix::Add(*jac,false,-*dt_,*eye_,false,1.0,tmpthetajac);
  ThetaJac->FillComplete();
  //  delete tmpthetajac;
  if (debug_) {
    std::string flnm = "Op2x.dat";
    EpetraExt::RowMatrixToMatlabFile(flnm.c_str(), *Op2x_);
    flnm = "Lplc.dat";
    EpetraExt::RowMatrixToMatlabFile(flnm.c_str(), *LinOp_);
    double nrm;
    rndmvct.Norm2(&nrm);
    std::cout << "\n norm of u = " << nrm << std::endl;
    Op2x_u->Norm2(&nrm);
    std::cout << "\n norm of Op2x_u = " << nrm << std::endl;
    std::cout << "\n frob. norm of Op2x = " << Op2x_->NormFrobenius()
              << std::endl;
    std::cout << "\n frob. norm of NlinJac = " << NlinJac->NormFrobenius()
              << std::endl;
    std::cout << "\n frob. norm of Jac = " << jac->NormFrobenius() << std::endl;
    std::cout << " frob. norm of ThetaJac = " << ThetaJac->NormFrobenius()
              << std::endl;
    getchar();
  }
  Prblm = Teuchos::rcp(new Epetra_LinearProblem);
  if (solver_type == "Amesos") {
    Prblm->SetOperator(ThetaJac.get());
    Amesos Factory;
    v_solve = Teuchos::rcp(Factory.Create(solver_name.c_str(), *Prblm), false);
    // v_solve = Teuchos::rcp(new Amesos_Klu(*v_prob));
    v_solve->SetUseTranspose(false);
    // Teuchos::ParameterList List;
    v_solve->SetParameters(*SolParams);
    Epetra_Time timer(ThetaJac->Comm());
    v_solve->SymbolicFactorization();
    v_solve->NumericFactorization();
    double time;
    if (test_) {
      double LocTime = timer.ElapsedTime();
      LinOp_->Comm().SumAll(&LocTime, &time, 1);
      if (ThetaJac->Comm().MyPID() == 0) {
        std::cout << "\n*******************************************\n";
        std::cout << "Average time per proc. for factorization"
                  << " = " << (time) / (LinOp_->Comm().NumProc()) << " sec\n";
        std::cout << "*******************************************\n"
                  << std::endl;
      }
    }
  }
  else if (solver_type == "Amesos2") {
    amesos2_solve = Amesos2::create<MAT,MV>(solver_name,ThetaJac);
    amesos2_solve->setParameters(SolParams); 
    Epetra_Time timer(ThetaJac->Comm());
    amesos2_solve->symbolicFactorization();
    amesos2_solve->numericFactorization();
    double time;
    if (test_) {
      double LocTime = timer.ElapsedTime();
      if (ThetaJac->Comm().MyPID() == 0) {
        std::cout << "\n*******************************************\n";
        std::cout << "*******************************************\n"
                  << std::endl;
      }
    }
  }
  /*   else
   *   {
   *     v_solve_iter = Teuchos::rcp(new AztecOO(*Prblm));
   *     v_solve_iter->SetUserMatrix(ThetaJac.get());
   *     v_solve_iter->SetParameters(*SolParams, true);
   *     MaxIter_ = SolParams->get("Max_Iter", 500);
   *     Tol_ = SolParams->get("Tol", .00000001);
   *     //v_solve_iter->SetAztecOption(AZ_precond, AZ_Jacobi);
   *     // v_solve_iter->CheckInput();getchar();
   *   }
   */

  //    eye_->Print(std::cout);
  //    std::cout<<"\n mass matrix above"<<std::endl;
  //    getchar();

} /* end of constructor */

void
Mean::createLinOp(Teuchos::RCP<Epetra_Comm> Comm)
{
  Teuchos::RCP<Epetra_Map> Map; // the map to be created
  double x_end = 1, x_in = 0;

  double dx = (x_end - x_in) / m;
  // std::string MapType = "Linear";
  Map = rcp(new Epetra_Map(m, 0, *Comm)); // =CreateMap(MapType, *Comm, List);
  Epetra_Vector x(*Map), ICom(*Map);
  int GlblRw;
  for (int i = 0; i < x.MyLength(); ++i) {
    GlblRw = x.Map().GID(i);
    x[i] = GlblRw * dx;
  }
  for (int i = 0; i < ICom.MyLength(); ++i) {
    ICom[i] =
      0.5 * (exp(cos(2 * PI * x[i])) - 1.5) * sin(2 * PI * (x[i] + 0.37));
  }
  Epetra_CrsMatrix* LplcOp =
    new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3);
  Epetra_CrsMatrix* GradOp =
    new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 2);
  Epetra_CrsMatrix eye(Epetra_DataAccess::Copy, *Map, 1);
  int NumGlobalElements = Map->NumGlobalElements();
  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = nullptr;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  std::vector<double> Values(2), GradVal(2);
  std::vector<int> Indices(2);
  int NumEntries;
  double val = 1.0;

  for (int i = 0; i < NumMyElements; ++i) {
    if (MyGlobalElements[i] == 0) {
      // off-diagonal for first row
      Indices[0] = 1;
      Values[0] = 1.0;
      Indices[1] = NumGlobalElements - 1;
      Values[1] = 1.0;
      NumEntries = 2;
      GradVal[0] = 1.0;
      GradVal[1] = -1.0;
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements - 2;
      Values[0] = 1.0;
      Indices[1] = 0;
      Values[1] = 1.0;
      NumEntries = 2;
      GradVal[0] = -1.0;
      GradVal[1] = 1.0;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1] = 1.0;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0] = 1.0;
      NumEntries = 2;
      GradVal[0] = -1.0;
      GradVal[1] = 1.0;
    }
    LplcOp->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0],
                               &Indices[0]);
    GradOp->InsertGlobalValues(MyGlobalElements[i], NumEntries, &GradVal[0],
                               &Indices[0]);
    // Put in the diagonal entry
    Values[0] = -2.0;
    LplcOp->InsertGlobalValues(MyGlobalElements[i], 1, &Values[0],
                               MyGlobalElements + i);

    eye.InsertGlobalValues(MyGlobalElements[i], 1, &val, MyGlobalElements + i);
  }

  LplcOp->FillComplete();
  LplcOp->OptimizeStorage();
  GradOp->FillComplete();
  GradOp->OptimizeStorage();
  eye.FillComplete();
  eye.OptimizeStorage();

  LplcOp->Scale(mu / pow(dx, 2));
  GradOp->Scale(1.0 / (2 * dx));

  *LinOp_ = *LplcOp;
  *Op1x_ = eye;
  *Op2x_ = *GradOp;
  *eye_ = eye;
  *u_ = ICom;
  *u0 = ICom;
  *x_ = x;
}

void
Mean::createW()
{
  // NumStchFrcVec_=1;
  W_ = rcp(new Epetra_MultiVector(x_->Map(), NumStchFrcVec_));
  for (int i = 0; i < W_->NumVectors(); ++i) {
    for (int j = 0; j < W_->MyLength(); ++j) {
      (*(*W_)(i))[j] = 0.5 * cos(4 * PI * (*x_)[j]);
    }
  }
}

void
Mean::BilinearTerm(RCP<Epetra_Vector> u1, RCP<Epetra_Vector> u2,
                   RCP<Epetra_Vector> u3)
{
  Op1x_->Multiply(false, *u1, *Op1x_u);
  Op2x_->Multiply(false, *u2, *Op2x_u);
  u3->Multiply(1.0, *Op1x_u, *Op2x_u, 0.0);
  // u3->Scale(-1.0);
}

/*
 * ===  FUNCTION
 * ======================================================================
 *         Name:  createRHS
 *  Description: create RHS such that du = dt * rhs  and rhs is calculated at u
 * =====================================================================================
 */
void
Mean::createRHS(RCP<Epetra_Vector> u)
{
  // EpetraExt::MatrixMatrix::Add(*LinOp_,false,1.0,*LinJac,0.0);
  LinOp_->Multiply(false, *u, *rhs); /* rhs = LinJac * u */
  uu->Multiply(1.0, *u, *u, 0.0);
  Op2x_->Multiply(false, *uu, *tmp2);
  tmp2->Scale(0.5);              /* tmp2 = 0.5 * Op2x * (u .* u) */
  rhs->Update(-1.0, *tmp2, 1.0); /* rhs = rhs - tmp2 - ExpVyVy */
  if (debug_) {
    double nrm;
    rhs->Norm2(&nrm);
    std::cout << "\n Mean.cpp:: norm of linJac*X + 0.5 * Op2x * (u .* u)  = "
              << nrm;

    //    getchar();
  }
  rhs->Update(-1.0, *ExpVyVy_, 1.0); /* rhs = rhs - tmp2 - ExpVyVy */
  // ExpVyVy_->Scale(0.0); /* used for testing and to check one way coupling */
  //  rhs->Update(-1.0,*tmp2,-1.0,*ExpVyVy_,1.0);/* rhs = rhs - tmp2 - ExpVyVy
  //  */
  rhs->Scale(*dt_); /* rhs = dt * rhs */

} /* end of createRHS() */
/*************************************************/
void
Mean::computeJac(RCP<Epetra_Vector> u)
{
  /* Calculate Op1x_u = Op1x * u ans Op2x_u = Op2x*u */
  Op1x_->Multiply(false, *u, *Op1x_u);
  Op2x_->Multiply(false, *u, *Op2x_u);
  if (debug_) {
    // Op2x_->Print(std::cout);
    // std::cout<<"\nOperator2 above!!!"<<std::endl; getchar();
  }
  // jac->Scale(0.0);
  /* calculate NlinJac= D( Op1x_u , Op2x_u )/Du */
  SpDiag->ReplaceDiagonalValues(*Op1x_u);
  EpetraExt::MatrixMatrix::Multiply(*SpDiag, false, *Op2x_, false, *TmpMat1);
  SpDiag->ReplaceDiagonalValues(*u);
  EpetraExt::MatrixMatrix::Multiply(*Op2x_, false, *SpDiag, false, *TmpMat2);
  EpetraExt::MatrixMatrix::Add(*TmpMat2, false, 1.0, *NlinJac, 0.0);
  //  EpetraExt::MatrixMatrix::Add(*TmpMat1, false, 1.0, *NlinJac, 1.0);
  /* jac = D( Op1x_u , Op2x_u )/Du + LinOp = LinOp - NlinJac */
  EpetraExt::MatrixMatrix::Add(*LinOp_, false, 1.0, *jac, 0.0);
  EpetraExt::MatrixMatrix::Add(*NlinJac, false, -1.0, *jac, 1.0);
  if (debug_) {
    std::cout.precision(18);
    double nrm;
    u->Norm2(&nrm);
    std::cout << "\nnorm of X : " << nrm;
    std::cout << "\nfrob norm of det-jac linear part : "
              << LinOp_->NormFrobenius();
    std::cout << "\nfrob norm of det-jac non-linear part : "
              << NlinJac->NormFrobenius();
  }
  /* ThetaJac = I - dt*jac  */
  EpetraExt::MatrixMatrix::Add(*jac, false, (*dt_) * theta, *ThetaJac, 0.0);
  EpetraExt::MatrixMatrix::Add(*eye_, false, 1.0, *ThetaJac, -1.0);
}
void
Mean::ThetaStepper()
{
  /* calculate rhs(x_t) */
  createRHS(u_);
  ThetaRHS->Update(theta, *rhs, 0.0);
  if (debug_) {
    std::cout.precision(18);
    double nrm;
    u_->Norm2(&nrm);
    std::cout << "\nnorm of X : " << nrm;
    u0->Norm2(&nrm);
    std::cout << "\nnorm of X0 : " << nrm;
    rhs->Norm2(&nrm);
    std::cout << "\nnorm of eq-rhs : " << nrm;
  }

  createRHS(u0);
  ThetaRHS->Update(1 - theta, *rhs, 1.0);
  ThetaRHS->Update(1.0, *u_, -1.0, *u0, -1.0);

  if (debug_) {
    std::cout.precision(18);
    double nrm;
    ThetaRHS->Norm2(&nrm);
    std::cout << "\nnorm of rhs : " << nrm;
    //    ThetaRHS->Print(std::cout<<"\n printing RHS\n");
    getchar();
  }
}
int
Mean::LinSolve(Epetra_Vector* LHS, Epetra_Vector* RHS)
{
  if (solver_type == "Amesos")
  {
  Prblm->SetLHS(LHS);
  Prblm->SetRHS(RHS);
  AMESOS_CHK_ERR(v_solve->NumericFactorization());
  AMESOS_CHK_ERR(v_solve->Solve());
  }
  else if (solver_type == "Amesos2")
  {
    amesos2_solve->numericFactorization();
    amesos2_solve->solve(LHS,RHS);
  }
  return 0;
}
/* data */
bool
Mean::NewtonSolver()
{
  isConverged_ = false;
  dx_->PutScalar(0.0);
  // y_jac();
  ThetaStepper();
  ThetaRHS->Scale(-1.0);
  ThetaRHS->Norm2(&NormRHS_);
  if (debug_)
    std::cout << "\nMean:: Newton:      initial norm: " << NormRHS_;
  for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
    computeJac(u_);
    LinSolve(dx_.get(), ThetaRHS.get());
    u_->Update(1.0, *dx_, 1.0);
    ThetaStepper();
    ThetaRHS->Scale(-1.0);
    ThetaRHS->Norm2(&NormRHStest_);
    if (debug_) {
      std::cout << "\nNewton:      iter: " << iter_;
      std::cout << "\nNewton:      norm: " << NormRHStest_;
    }
    if (NormRHStest_ < toleranceRHS_) {
      if (debug_) {
        std::cout << "\nSuccess...";
      }
      break;
    }
    if (backTracking_ and (NormRHS_ < NormRHStest_))
      RunBackTracking();
    NormRHS_ = NormRHStest_;
  }
  if (iter_ == maxNumIterations_) {
    std::cout << "\nNewton: ---> TROUBLE" << __FILE__ << __LINE__;
    std::cout << "\nNewton not converged after Max number of iter!!!!!!\n";
    std::cout << "\nrhs norm = " << NormRHStest_;
    getchar();
  } else {
    *u0 = *u_;
    isConverged_ = true;
  }
  return isConverged_;
}
//======================================================================
void
Mean::RunBackTracking()
{
  // Initialize reduction with -1/2
  double reduction = -1.0 / 2;
  for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_) {
    if (NormRHStest_ < NormRHS_) {
      if (debug_) {
        std::cout << "\nSuccess...";
      }
      break;
    }
    // Apply reduction to the state vector
    u_->Update(reduction, *dx_, 1.0);
    ThetaStepper();
    rhs->Norm2(&NormRHStest_);
    //		std::cout<<"Newton: --> backtracking:\n "
    //			 <<	" step: "      << backTrack_
    //			 << "\n reduction: " << reduction
    //			 << "\n norm: "      << NormRHStest_;
    // Update reduction
    reduction /= 2.0;
  }
  if (backTrack_ == numBackTrackingSteps_)
    std::cout << "\nNewton: --> BACKTRACKING FAILED" << __FILE__ << __LINE__;
}
void
Mean::WriteSolution(std::string filename, double param,
                    const Epetra_Vector& soln)
{
  Teuchos::RCP<std::ostream> out;
  if (soln.Comm().MyPID() == 0) {
    out = Teuchos::rcp(new std::ofstream(filename.c_str()));
  } else { // dummy stream
    out = Teuchos::rcp(new Teuchos::oblackholestream());
  }
  (*out) << std::setw(15) << std::setprecision(15);
  out->setf(std::ios::scientific);
  (*out) << param;
  (*out) << *(HYMLS::MatrixUtils::Gather(soln, 0));
}
