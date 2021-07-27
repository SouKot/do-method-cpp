/*
 * =====================================================================================
 *
 *       Filename:  Mean.hpp
 *
 *    Description:  
 *
 *
 *
 *    class for solving mean of SDE
 *
 *        Version:  1.0
 *        Created:  04/15/2018 03:38:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef mean_solve_h
#define mean_solve_h
//#include <BelosEpetraAdapter.hpp>
//#include <BelosIMGSOrthoManager.hpp>
//#include <BelosMultiVec.hpp>
//#include <Teuchos_SerialDenseMatrix.hpp>
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
//#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
//#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Ifpack_Preconditioner.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_LinearProblem.h"
#include "Amesos.h"
#include "Amesos_ConfigDefs.h"  //#include <random>
#include "Amesos_BaseSolver.h"
// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.
#include "AnasaziEpetraAdapter.hpp"
#include "QG.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
class AztecOO;
using namespace Teuchos;
class Mean
{
  public:
   Mean(
        RCP<Teuchos::ParameterList> PrmLst,
        double* t, double* dt,
	RCP<Epetra_Comm> Comm);
  void createLinOp(Teuchos::RCP<Epetra_Comm> Comm);
void BilinearTerm(RCP<Epetra_Vector> u1,
			RCP<Epetra_Vector> u2,
			RCP<Epetra_Vector> u3);
  void createRHS(RCP<Epetra_Vector> x);
  void computeF(Epetra_Vector& u,Epetra_Vector& F)
  { 
    createRHS(Teuchos::rcpFromRef(u));
    F=*rhs;
  }
  void getProblemRHS(Epetra_Vector& u,Epetra_Vector& F);
  Epetra_Vector& getF(){return *rhs;}
  void computeJac(RCP<Epetra_Vector> u);
  void computeJacobian(Epetra_Vector& x,Epetra_CrsMatrix& A)
  {
    computeJac(Teuchos::rcpFromRef(x));
    A=*jac;
  }
  RCP<Epetra_CrsMatrix> getJacobian(){ return jac;}
  void ThetaStepper(Teuchos::RCP<Epetra_Vector> rhs_u0);
  int LinSolve(Epetra_Vector* LHS, Epetra_Vector* RHS);
  bool NewtonSolver();
  void RunBackTracking(Epetra_Vector &rhs_u0);
  RCP<Epetra_Vector> get_Xdim(){return x_;}
  RCP<Epetra_Vector> getSolution(){return u_;}
  void setExpVyVy(RCP<Epetra_Vector>ExpVyVy){ExpVyVy_=ExpVyVy;}
  RCP<Epetra_CrsMatrix> getMassMatrix();
  void WriteSolution(std::string filename, double param,
                         const Epetra_Vector& soln);
  void createW(Epetra_Vector& diagmass);
  int get_dim_W(){return NumStchFrcVec_;}
  RCP<Epetra_MultiVector> get_W(){return W_;}
  double theta;
  double Tol;
  double MaxIter;
  int n, ny, nx;
  double mu,rynldsNum;
  RCP<Epetra_Vector> u0;
  RCP<Epetra_Vector> RHS;
  RCP<Epetra_CrsMatrix> ThetaJac;
  RCP<Epetra_CrsMatrix> jac;
  RCP<Epetra_Vector> rhs, uu, tmp2, ThetaRHS;
  RCP<Epetra_LinearProblem> Prblm;
  Teuchos::RCP<Amesos_BaseSolver> v_solve;
  Teuchos::RCP<AztecOO> v_solve_iter;
  Teuchos::RCP<QG::QG> qg;
  //virtual ~Mean();
  void printnormMV(Epetra_MultiVector &mv, int normType, string str);
  bool newtonLineSearchSolve(Epetra_Vector &x0);
  void setSolution(Epetra_Vector &x);
private:
  RCP<Epetra_Comm> Comm_;
  RCP<Epetra_CrsMatrix> mass_;
  RCP<Epetra_Vector> u_,dx_,x_; /* dx_ should be du_ */
  RCP<Epetra_Vector> ExpVyVy_;
  RCP<Epetra_MultiVector> W_;
  bool isConverged_,backTracking_,test_,debug_;
  int NumStchFrcVec_;
  std::string solver_type;

  //Teuchos::ParameterList MeanParam_;
  double *t_, *dt_;
  double NormRHS_,NormRHStest_,toleranceRHS_;
  int iter_,maxNumIterations_,backTrack_,numBackTrackingSteps_;  int MyPID;
  
};
#endif

