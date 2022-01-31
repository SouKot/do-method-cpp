/*
 * =====================================================================================
 *
 *       Filename:  Mean.cpp
 *
 *    Description:
 *
 *    Calculates mean of QG model
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
#include "EpetraExt_MultiVectorOut.h"
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
  , toleranceRHS_(1.0e-10)
  , maxNumIterations_(50)
  , numBackTrackingSteps_(20)
  , backTracking_(true)
{
  Teuchos::RCP<Epetra_Map> Map;

  nx = PrmLst->get("nx", 200); // global number of elements in the map
  ny = PrmLst->get("ny", 200); // global number of elements in the map
  rynldsNum = PrmLst->get("Reynolds Number",0.0); // Reynolds Numbers
  double topography = PrmLst->get("Topography",0.0); // Reynolds Numbers
  n = nx * ny * 2;
  qg = Teuchos::rcp(new QG::QG(nx, ny));
  double zeta = 1 / n;
  qg->set_par(20, zeta);
  qg->set_par(21, topography);
  qg->set_par(5, rynldsNum);
  qg->set_par(11, 1.0);
  theta = PrmLst->get("theta", 0.5);
  NumStchFrcVec_ = PrmLst->get("No. of vectors in stoch. forcing",
                               2); // global number of elements in the map
  solver_type = PrmLst->get("Solver Package", "Amesos2");
  test_ = PrmLst->get("Testing", true);
  debug_ = PrmLst->get("Debugging", true);
  if (debug_) {
    std::cout << "\n nx = " << n << std::endl;
    std::cout << "\n mu = " << mu << std::endl;
    std::cout << "\n numvecV = " << NumStchFrcVec_ << std::endl;
  }
  Map = rcp(new Epetra_Map(n, 0, *Comm)); // =CreateMap(MapType, *Comm, List);
  mass_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
  u_ = rcp(new Epetra_Vector(*Map));
  double diagmassval[n];
  qg->mass(diagmassval);
  Epetra_Vector diagMass(Epetra_DataAccess::Copy, *Map, diagmassval);
  for (int i = 0; i < n; i++) {
    mass_->InsertGlobalValues(i, 1, &diagMass[i], &i);
  }
  mass_->FillComplete();
  // mass_->ReplaceDiagonalValues(diagMass);
  mass_->Scale(-1.0); // NOTE: mass matrix multiplied by -1.0 !! see QG.m

  u0 = rcp(new Epetra_Vector(*u_));
  x_ = rcp(new Epetra_Vector(*u_));

  createW(diagMass);
  MyPID = mass_->Comm().MyPID();
  Teuchos::ParameterList& SolverSublist = PrmLst->sublist(solver_type);
  std::string solver_name =
    SolverSublist.get("Solver name", "Mumps");
  const Teuchos::RCP<Teuchos::ParameterList> SolParams =
    Teuchos::getParametersFromXmlFile(solver_name + ".xml");
  if (MyPID == 0)
    std::cout << solver_type << " solver(" << solver_name
              << ") has been chosen as mean solver";

  ThetaJac = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  jac = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
  RHS = rcp(new Epetra_Vector(*u_));
  dx_ = rcp(new Epetra_Vector(*u_));
  rhs = rcp(new Epetra_Vector(*RHS));
  ExpVyVy_ = rcp(new Epetra_Vector(*RHS));
  uu = rcp(new Epetra_Vector(*RHS));
  ThetaRHS = rcp(new Epetra_Vector(*RHS));
  Epetra_Vector rndmvct(*u_);
  u_->Random();
  computeJac(u_); // create jac and Thetajac first time.
  ThetaJac->FillComplete();
  u_->PutScalar(0.0);
  *u0 = *u_;
  if (debug_) {
    double nrm;
    u_->Norm2(&nrm);
    std::cout << "\n norm of u = " << nrm << std::endl;
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
      if (ThetaJac->Comm().MyPID() == 0) {
        std::cout << "\n*******************************************\n";
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
//===============================================================
void
Mean::createW(Epetra_Vector& diagmass)
{
  double l = 0.125;
  double C = 1.0;
  double NX = sqrt(n / 2);
  double etax, etay, j, rem, q;
  //NumStchFrcVec_ = 2; //NOTE: Hard coded right now!!
  W_ = rcp(new Epetra_MultiVector(x_->Map(), NumStchFrcVec_));
  for (int col = 0; col < W_->NumVectors(); ++col) {
    etax = C * col;
    etay = C * (1 - col);
    for (int i = 1; i <= W_->MyLength(); i = i + 2) {
      j = i / 2.0;
      rem = fmod(j, NX) / (NX - 1);
      q = floor(j / NX) / (NX - 1);
      (*(*W_)(col))[i - 1] =
        exp(-2. * ((rem - 0.5) * (rem - 0.5) + (q - 0.5) * (q - 0.5)) /
            (4. * l * l)) *
        (etay - 2. * etay * rem + 2. * etax * q - etax) / (2. * l * l);
      if (abs(diagmass[i - 1]) < 1.0e-14)
        (*(*W_)(col))[i - 1] = 0.0;
      (*(*W_)(col))[i] = 0.0;
    }
  }
  const char* flnm="forcing.mm";
  EpetraExt::MultiVectorToMatrixMarketFile(flnm,*W_);
  if (debug_)
    printnormMV(*W_, 2, "norm of Original W :");
}
//========================================================================
void
Mean::BilinearTerm(RCP<Epetra_Vector> u1, RCP<Epetra_Vector> u2,
                   RCP<Epetra_Vector> u3)
{
  double *viewU1, *viewU2, *viewU3;
  u1->ExtractView(&viewU1);
  u2->ExtractView(&viewU2);
  u3->ExtractView(&viewU3);
  qg->bilin(viewU1, viewU2, viewU3);
  // u3->Scale(1.0);
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
  double* viewRHS;
  double* viewX;
  u->ExtractView(&viewX);
  rhs->ExtractView(&viewRHS);
  qg->rhs(viewX, viewRHS);
  rhs->Update(1.0, *ExpVyVy_, -1.0); /* f1=ExpVyVy-rhs!! check stochastic.m */
  // EpetraExt::MatrixMatrix::Add(*LinOp_,false,1.0,*LinJac,0.0);
  if (debug_) {
    double nrm;
    printnormMV(*rhs, 2, "norm of qg.rhs");
    printnormMV(*ExpVyVy_, 2, "Mean.cpp:: norm of ExpVyVy");
  }
  // ExpVyVy_->Scale(0.0); /* used for testing and to check one way coupling */
  //  rhs->Update(-1.0,*tmp2,-1.0,*ExpVyVy_,1.0);/* rhs = rhs - tmp2 - ExpVyVy
  //  */

  rhs->Scale(*dt_); /* rhs = dt * rhs */
}
//======================================================================
void
Mean::setSolution(Epetra_Vector& x)
{
  *u0 = x;
  *u_ = x;
}
//======================================================================
void
Mean::getProblemRHS(Epetra_Vector& u, Epetra_Vector& F)
{
  double* viewRHS;
  double* viewX;
  u.ExtractView(&viewX);
  F.ExtractView(&viewRHS);
  qg->rhs(viewX, viewRHS);

} /* end of createRHS() */
//======================================================================
void
Mean::computeJac(RCP<Epetra_Vector> u)
{
  int mxcolsz = n * 20; // hoping this would be enough
  int numentries;
  double* viewU;
  int* rowmrkr = new int[n + 1];
  int* indexOffSet;
  int* indices;
  double* values;
  int* col = new int[mxcolsz];
  double* val = new double[mxcolsz];
  u->ExtractView(&viewU);
  qg->jacobian(viewU, 0.0, rowmrkr, col, val);
  if (jac->Filled()) {
    jac->PutScalar(0.0);
    jac->ExtractCrsDataPointers(indexOffSet, indices, values);
    for (int i = 0; i < n + 1; ++i) {
      indexOffSet[i] = rowmrkr[i];
    }
    for (int i = 0; i < indexOffSet[n]; ++i) {
      indices[i] = col[i];
      values[i] = val[i];
    }
  } else {
    for (int i = 0; i < n + 1; ++i) {
      numentries = rowmrkr[i + 1] - rowmrkr[i];
      if (numentries == 0)
        continue;
      for (int j = rowmrkr[i]; j < rowmrkr[i + 1]; ++j)
        jac->InsertGlobalValues(i, 1, &val[j], &col[j]);
    }
  }
  if (!jac->Filled())
    jac->FillComplete();

  jac->Scale(-1.0); // NOTE: jac is multiplied by -1.0 !!! see QG.m
  if ((0)) {        // NOTE : switch it off!!
    double nrm;
    u->Norm2(&nrm);
    std::cout << " \n norm of u : " << nrm;
    std::cout << "\n rowmrkr[n] :" << rowmrkr[n];
    int nn = 20;
    std::cout << "\n In" << __PRETTY_FUNCTION__ << " : " << __LINE__
              << "\n values in rowmrkr:  ";
    for (int i = 0; i < nn; ++i) {
      std::cout << rowmrkr[i] << "  ";
    }
    std::cout << "\n values in col:  ";
    for (int i = 0; i < nn; ++i) {
      std::cout << col[i] << "  ";
    }
    std::cout << "\n values in val:  ";
    for (int i = 0; i < nn; ++i) {
      std::cout << val[i] << "  ";
    }
    getchar();
  }
  delete[] rowmrkr;
  delete[] col;
  delete[] val;
  /* Calculate Op1x_u = Op1x * u ans Op2x_u = Op2x*u */
  if (debug_) {
    std::cout.precision(18);
    double nrm;
    u->Norm2(&nrm);
    std::cout << "\nnorm of X : " << nrm;
  }
  /* ThetaJac = I - dt*jac  */
  EpetraExt::MatrixMatrix::Add(*jac, false, (*dt_) * theta, *ThetaJac, 0.0);
  EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *ThetaJac, -1.0);
  if (debug_) {
    std::cout << "\n Frob. norm of det-jac : " << ThetaJac->NormFrobenius()
              << std::endl;
  }
}
//=============================================================================
void
Mean::ThetaStepper(Teuchos::RCP<Epetra_Vector> rhs_u0)
{
  /* FIXME : even when thea=[1 or 0], it calls createRHS twice, which just waste
   * time. */
  createRHS(u_); /* calculate rhs(x_t) */
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

  // createRHS(u0);
  ThetaRHS->Update(1 - theta, *rhs_u0, 1.0);
  Epetra_Vector U(*u_), MU(*u_);
  U.Update(1.0, *u_, -1.0, *u0, 0.0); // U = u - u0
  mass_->Multiply(false, U, MU);      // MU = M * U = M * (u - u0)
  ThetaRHS->Update(1.0, MU, -1.0);

  if (debug_) {
    std::cout.precision(18);
    double nrm;
    ThetaRHS->Norm2(&nrm);
    std::cout << "\nnorm of rhs : " << nrm;
    //    ThetaRHS->Print(std::cout<<"\n printing RHS\n");
    getchar();
  }
}
//========================================================================
int
Mean::LinSolve(Epetra_Vector* LHS, Epetra_Vector* RHS)
{
  if (solver_type == "Amesos")
  {
  Prblm->SetLHS(LHS);
  Prblm->SetRHS(RHS);
  // AMESOS_CHK_ERR(v_solve->SymbolicFactorization());
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
//========================================================================
bool
Mean::NewtonSolver()
{
  isConverged_ = false;
  dx_->PutScalar(0.0);
  createRHS(u0);              // creates rhs for u0
  Epetra_Vector rhs_u0(*rhs); // rhs_u0 = rhs(u0)
  ThetaStepper(Teuchos::rcpFromRef(rhs_u0));
  ThetaRHS->Scale(-1.0);
  ThetaRHS->Norm2(&NormRHS_);
  if (debug_)
    std::cout << "\nMean:: Newton:      initial norm: " << NormRHS_;
  for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
    computeJac(u_);
    LinSolve(dx_.get(), ThetaRHS.get());
    //////////////////////////////////////////////
    double s, s0, etasq;
    Epetra_Vector JU(*dx_);
    ThetaJac->Multiply(false, *dx_, JU);
    dx_->Dot(JU, &s);
    ThetaRHS->Dot(*dx_, &s0);
    etasq = (s0 / s);
    if (debug_) {
      std::cout << "\n value of s0: " << s0;
      std::cout << "\n value of s: " << s;
      std::cout << "\n value of eta: " << etasq;
    }
    //////////////////////////////////////////////
    if (etasq > 0.0)
      u_->Update(sqrt(etasq), *dx_, 1.0);
    else
      u_->Update(sqrt(-etasq), *dx_, 1.0);

    ThetaStepper(Teuchos::rcpFromRef(rhs_u0));
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
    // if (backTracking_ and (NormRHS_ < NormRHStest_))
    // RunBackTracking(rhs_u0);
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
Mean::RunBackTracking(Epetra_Vector& rhs_u0)
{
  // Initialize reduction with -1/2
  std::cout << "\n Running backtracking for solving Mean!!";
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
    ThetaStepper(Teuchos::rcpFromRef(rhs_u0));
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
//===========================================================================
bool
Mean::newtonLineSearchSolve(Epetra_Vector& x0)
{
  int nIter = 0;
  double alpha = 1.0e-4, maxiter = 30, lambdaMax = 0.5, lambdaMin = 0.1;
  double tolX = 1.0e-12, tolFun = 1.0e-8, resnrm, lambda = 1.0;
  double fold = 0.0, dscrmnnt, resnrm0, convergence, stepnorm;
  double slope = 0.0, mxval, f, f2 = 0.0, lambdaPre, lambda2 = 0.0, A, c1, c2,
         a, b;
  isConverged_ = false;
  Epetra_Vector x(x0), gT(x0), xold(x0), xoldMod(x0);
  createRHS(Teuchos::rcpFromRef(x0)); // creates rhs for u0
  Epetra_Vector rhs_u0(*rhs);         // rhs_u0 = rhs(u0)
  ThetaStepper(Teuchos::rcpFromRef(rhs_u0));
  ThetaRHS->Scale(-1.0);
  computeJac(Teuchos::rcpFromRef(x));
  ThetaRHS->Norm2(&resnrm);
  dx_->PutScalar(0.0);
  if (debug_)
    std::cout << "\n ITER. NO.    RES. NORM    STEP NORM     LAMBDA ";
  while ((resnrm > tolFun or lambda < 1) && nIter <= maxiter) {
    if (abs(lambda - 1.0) < 1.0e-10) {
      nIter += 1;
      LinSolve(dx_.get(), ThetaRHS.get());
      ThetaJac->Multiply(true, *ThetaRHS, gT);
      gT.Dot(*dx_, &slope);
      ThetaRHS->Dot(*ThetaRHS, &fold);
      xold = x;
      for (int i = 0; i < xold.MyLength(); ++i) {
        xoldMod[i] = abs((*dx_)[i]) / std::max(abs(xold[i]), 1.0);
      }
      xoldMod.MaxValue(&mxval);
      lambdaMin = tolX / mxval;
      // std::cout << "\n frob. norm of jacobian: " <<
      // ThetaJac->NormFrobenius();
      // printnormMV(*dx_, 2, "norm of dx: ");
      // std::cout << "\n value of fold: " << fold;
    }
    x.Update(1.0, xold, lambda, *dx_, 0.0);
    *u_ = x;
    ThetaStepper(Teuchos::rcpFromRef(rhs_u0));
    ThetaRHS->Scale(-1.0);
    computeJac(Teuchos::rcpFromRef(x));
    ThetaRHS->Dot(*ThetaRHS, &f);
    lambdaPre = lambda;
    // std::cout << "\n value of f: " << f;
    // std::cout << "\n value of slope: " << slope;
    // getchar();
    if (f > fold + alpha * lambda * slope) {
      if (abs(lambda - 1.0) < 1.0e-12) {
        // std::cout << "\n f-fold-slope = " << f - fold - slope;
        // std::cout << "\n f-fold-alpha*lambda*slope = "
        // << f - fold - alpha * lambda * slope;
        lambda = -slope / (2 * (f - fold - slope));
        // std::cout << "\n new lambda = " << lambda;
      } else {
        A = 1.0 / (lambdaPre - lambda2);
        c1 = f - fold - lambdaPre * slope;
        c2 = f2 - fold - lambda2 * slope;
        a =
          A * (1 / pow(lambdaPre, 2.0) * c1 + (-1.0) / pow(lambda2, 2.0) * c2);
        b = A * (-lambda2 / pow(lambdaPre, 2.0) * c1 +
                 lambdaPre / pow(lambda2, 2.0) * c2);
        if (abs(a) < 1.0e-12) {
          lambda = -slope / (2 * b);
        } else {
          dscrmnnt = pow(b, 2.0) - 3 * a * slope;
          if (dscrmnnt < 0)
            lambda = lambdaMax * lambdaPre;
          else if (b <= 0.0)
            lambda = (-b + sqrt(dscrmnnt)) / (3 * a);
          else
            lambda = -slope / (b + sqrt(dscrmnnt));
        }
        lambda = std::min(lambda, lambdaMax * lambdaPre);
      }
    } else if (std::isnan(f) or std::isinf(f)) {
      lambda = lambdaMax * lambdaPre;
    } else {
      lambda = 1.0;
    }
    if (lambda < 1.0) {
      lambda2 = lambdaPre;
      f2 = f;
      lambda = std::max(lambda, lambdaMin * lambdaPre);
      continue;
    }
    resnrm0 = resnrm;
    ThetaRHS->Norm2(&resnrm);
    convergence = log(resnrm0 / resnrm);
    dx_->Norm2(&stepnorm);
    if (debug_) {
      std::cout << "\n"
                << "   " << nIter << "    " << resnrm << "    " << stepnorm
                << "    " << lambda;
      printnormMV(*dx_, 2, "dx = ");
      printnormMV(xold, 2, "xold = ");
      printnormMV(x, 2, "x = ");
    }
  }
  if (resnrm <= tolFun)
    isConverged_ = true;
  else {
    std::cout << "\n WARNING: NEWTON SOLVER FOR MEAN NOT CONVERGED!!";
    getchar();
  }

  setSolution(x);
  return isConverged_;
}
//===========================================================================
RCP<Epetra_CrsMatrix>
Mean::getMassMatrix()
{
  return mass_;
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
void
Mean::printnormMV(Epetra_MultiVector& mv, int normType, string str)
{
  double nrm1[mv.NumVectors()];
  if (normType == 1)
    mv.Norm1(&nrm1[0]);
  else if (normType == 2)
    mv.Norm2(&nrm1[0]);
  std::cout << "\n" << str << std::endl;
  for (int i = 0; i < mv.NumVectors(); i++)
    std::cout << "  " << nrm1[i];
}
