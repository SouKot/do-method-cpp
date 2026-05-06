/*
 * =====================================================================================
 *
 *       Filename:  Mean.hpp
 *
 *    Description:  Mean-field solver for the stochastic Quasi-Geostrophic (QG)
 *                  equation.  Integrates the mean PDE with a theta-method time
 *                  discretisation and Newton-Raphson nonlinear solve.  The
 *                  stochastic forcing W is time-invariant and set once in the
 *                  constructor; refreshForcing() is a no-op for this model.
 *
 *        Version:  1.0
 *        Created:  04/15/2018 03:38:12 PM
 *       Revision:  none
 *       Compiler:  gcc / clang (C++17)
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
#include "../DO/DOUtils.hpp"
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
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
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

/**
 * @brief Mean-field solver for the stochastic Quasi-Geostrophic equation.
 *
 * @details Integrates the deterministic (mean) part of the stochastic QG
 * equation using a theta-method time discretisation and Newton-Raphson
 * nonlinear solve.  The stochastic forcing W is computed once in the
 * constructor (time-invariant); refreshForcing() is therefore a no-op.
 * Used exclusively by the QGDO build target.
 *
 * @see Burger::Mean  Equivalent solver for the Burgers model.
 */
class Mean
{
  public:

  /**
   * @brief Construct and initialise the QG mean solver.
   *
   * Builds the mass matrix, Jacobian structure, and initial stochastic forcing W
   * from @p PrmLst.  The diagonal mass vector required by createW() is derived
   * internally from the QG operator.
   *
   * @param[in] PrmLst  Parameter list (nx, ny, Reynolds number, theta, solver).
   * @param[in] t       Pointer to the current simulation time (advanced externally).
   * @param[in] dt      Pointer to the current time-step size (updated externally).
   * @param[in] Comm    Epetra communicator (serial or MPI).
   */
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
  /**
   * @brief Advance the mean solution by one time step using Newton-Raphson.
   * @return @c true on convergence, @c false if the solver fails.
   */
  bool NewtonSolver();
  void RunBackTracking(Epetra_Vector &rhs_u0);
  RCP<Epetra_Vector> get_Xdim(){return x_;}

  /// @brief Return an RCP to the current mean solution vector.
  RCP<Epetra_Vector> getSolution(){return u_;}

  /**
   * @brief Inject the E[V y V y^T] term computed by the DO stochastic solver.
   *
   * Called once per time step by the DO driver after the stochastic sub-step
   * so that the mean RHS can include the second-moment correction.
   *
   * @param[in] ExpVyVy  E[Vy(Vy)^T] vector computed by Y_Stoch.
   */
  void setExpVyVy(RCP<Epetra_Vector>ExpVyVy){ExpVyVy_=ExpVyVy;}

  /// @brief Return an RCP to the (diagonal, scaled by -1) mass matrix.
  RCP<Epetra_CrsMatrix> getMassMatrix();

  /**
   * @brief Write the solution vector to a text file.
   * @param[in] filename  Output file path.
   * @param[in] param     Parameter value printed in the header (e.g. current time).
   * @param[in] soln      Solution vector to write.
   */
  void WriteSolution(std::string filename, double param,
                         const Epetra_Vector& soln);

  /**
   * @brief Compute the time-invariant stochastic forcing vector W.
   *
   * Fills W_ with a Gaussian-weighted divergence pattern.  Entries corresponding
   * to zero mass-matrix rows (boundary / pressure dofs) are zeroed via @p diagmass.
   * Called once in the constructor; should not be called again during time stepping
   * as it reallocates W_ and invalidates shared RCPs held by Y_Stoch / Problem_Interface.
   *
   * @param[in] diagmass  Diagonal of the QG mass matrix; used to mask zero-mass dofs.
   */
  void createW(Epetra_Vector& diagmass);

  /**
   * @brief No-op forcing-refresh adapter for the DO time loop.
   *
   * QG's forcing W is time-invariant and already constructed in the constructor,
   * so this method intentionally does nothing.  Provides the same call site as
   * Burger::Mean::refreshForcing() so the driver contains no model-specific guards.
   */
  void refreshForcing(double) {}

  /// @brief Return the number of stochastic forcing vectors (columns of W).
  int get_dim_W(){return NumStchFrcVec_;}

  /// @brief Return an RCP to the stochastic forcing multi-vector W.
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
  typedef Epetra_MultiVector MV;
  typedef Epetra_CrsMatrix MAT;
  Teuchos::RCP<Amesos2::Solver<MAT,MV> > amesos2_solve;
  Teuchos::RCP<QG::QG> qg;
  //virtual ~Mean();
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

