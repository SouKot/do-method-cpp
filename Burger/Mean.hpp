/*
 * =====================================================================================
 *
 *       Filename:  Mean.hpp
 *
 *    Description:  Mean-field solver for the stochastic Burgers equation.
 *                  Integrates the mean PDE with a theta-method time
 *                  discretisation and Newton-Raphson nonlinear solve with
 *                  optional back-tracking line search.  The time-dependent
 *                  stochastic forcing W is provided via createW() and the
 *                  unified refreshForcing() adapter.
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
//#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
//#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Ifpack.h"
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
#include "BelosConfigDefs.hpp"
#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.
#include "AnasaziEpetraAdapter.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
using namespace Teuchos;

/**
 * @brief Mean-field solver for the stochastic Burgers equation.
 *
 * @details Integrates the deterministic (mean) part of the stochastic Burgers
 * equation in time using a theta-method.  At each step Newton-Raphson is used
 * to solve the nonlinear system; optional back-tracking provides robustness.
 * Used exclusively by the BRGRDO build target.
 *
 * @see QG::Mean   Equivalent solver for the quasi-geostrophic model.
 */
class Mean {
public:
    /**
     * @brief Construct and initialise the Burgers mean solver.
     *
     * Builds the linear operator, mass matrix, and initial stochastic forcing W.
     *
     * @param[in] PrmLst  Parameter list (mesh size, Reynolds number, theta, solver).
     * @param[in] t       Pointer to the current simulation time (advanced externally).
     * @param[in] dt      Pointer to the current time-step size (updated externally).
     * @param[in] Comm    Epetra communicator (serial or MPI).
     */
    Mean(
        RCP<Teuchos::ParameterList> PrmLst,
        double *t, double *dt,
        RCP<Epetra_Comm> Comm);

    void createLinOp(const Teuchos::RCP<Epetra_Comm>& Comm);

    void BilinearTerm(
        RCP<Epetra_Vector> u1,
        RCP<Epetra_Vector> u2,
        RCP<Epetra_Vector> u3);

    void createRHS(RCP<Epetra_Vector> x);

    void computeF(Epetra_Vector &u, Epetra_Vector &F) {
        createRHS(Teuchos::rcpFromRef(u));
        F = *rhs;
    }

    Epetra_Vector &getF() { return *rhs; }

    void computeJac(RCP<Epetra_Vector> u);

    void computeJacobian(Epetra_Vector &x, Epetra_CrsMatrix &A) {
        computeJac(Teuchos::rcpFromRef(x));
        A = *jac;
    }

    RCP<Epetra_CrsMatrix> getJacobian() { return jac; }

    void ThetaStepper();

    int LinSolve(Epetra_Vector &LHS, Epetra_Vector &RHS);

    /**
     * @brief Advance the mean solution by one time step using Newton-Raphson.
     * @return @c true on convergence, @c false if the solver fails.
     */
    bool NewtonSolver();

    void RunBackTracking();

    RCP<Epetra_Vector> get_Xdim() { return x_; }

    /// @brief Return an RCP to the current mean solution vector.
    RCP<Epetra_Vector> getSolution() { return u_; }

    /**
     * @brief Inject the E[V y V y^T] term computed by the DO stochastic solver.
     *
     * Called once per time step by the DO driver after the stochastic sub-step
     * so that the mean RHS can include the second-moment correction.
     *
     * @param[in] ExpVyVy  E[Vy(Vy)^T] vector computed by Y_Stoch.
     */
    void setExpVyVy(RCP<Epetra_Vector> ExpVyVy) { ExpVyVy_ = ExpVyVy; }

    /// @brief Return an RCP to the mass matrix (identity for Burgers).
    RCP<Epetra_CrsMatrix> getMassMatrix() {
        return eye_;
    }

    /**
     * @brief Write the solution vector to a text file.
     * @param[in] filename  Output file path.
     * @param[in] param     Parameter value printed in the header (e.g. current time).
     * @param[in] soln      Solution vector to write.
     */
    void WriteSolution(std::string filename, double param,
                       const Epetra_Vector &soln);

    /**
     * @brief Compute the time-dependent stochastic forcing vector W at time @p t.
     *
     * Fills W_ with the harmonic forcing pattern @f$ 0.5[\cos(4\pi x)e^{-10t} +
     * \cos(2\pi x)e^{5(t-1)}] @f$.  W_ is allocated on first call and reused
     * in-place on subsequent calls, so shared RCPs obtained via get_W() remain valid.
     *
     * @param[in] t  Current simulation time.
     */
    void createW(double t = 0);

    /**
     * @brief Unified forcing-refresh adapter called by the DO time loop.
     *
     * Delegates to createW(t).  Provides a model-agnostic call site in the driver
     * so that no @c #if @c brgr / @c quasi_geo guard is needed in the time loop.
     *
     * @param[in] t  Current simulation time passed to createW().
     */
    void refreshForcing(double t) { createW(t); }

    /// @brief Return the number of stochastic forcing vectors (columns of W).
    int get_dim_W() { return NumStchFrcVec_; }

    /// @brief Return an RCP to the stochastic forcing multi-vector W.
    RCP<Epetra_MultiVector> get_W() { return W_; }
    double theta;
    double Tol;
    double MaxIter;
    int m;
    double mu;
    RCP<Epetra_Vector> u0;
    RCP<Epetra_Vector> RHS;
    RCP<Epetra_CrsMatrix> LinJac;
    RCP<Epetra_CrsMatrix> NlinJac;
    RCP<Epetra_CrsMatrix> ThetaJac;
    RCP<Epetra_CrsMatrix> jac;
    RCP<Epetra_CrsMatrix> SpDiag;
    RCP<Epetra_Vector> Op1x_u;
    RCP<Epetra_Vector> Op2x_u;
    RCP<Epetra_CrsMatrix> TmpMat1;
    RCP<Epetra_CrsMatrix> TmpMat2;
    RCP<Epetra_Vector> rhs, uu, tmp2, ThetaRHS;
    RCP<Epetra_LinearProblem> Prblm;
    Teuchos::RCP<Amesos_BaseSolver> v_solve;
    typedef Epetra_MultiVector MV;
    typedef Epetra_CrsMatrix MAT;
    Teuchos::RCP<Amesos2::Solver<MAT, MV> > amesos2_solve;
    typedef double ST;
    typedef Teuchos::ScalarTraits<ST> SCT;
    typedef SCT::magnitudeType MT;
    typedef Epetra_Operator OP;
    typedef Belos::MultiVecTraits<ST, MV> MVT;
    typedef Belos::OperatorTraits<ST, MV, OP> OPT;
    Teuchos::RCP<Ifpack_Preconditioner> prec;
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec;
    Teuchos::RCP<Belos::LinearProblem<double, MV, OP> > problem;
    Teuchos::RCP<Belos::SolverManager<double, MV, OP> > v_solve_iter;

    //virtual ~Mean();
private:
    RCP<Epetra_Comm> Comm_;
    RCP<Epetra_CrsMatrix> LinOp_;
    RCP<Epetra_CrsMatrix> Op1x_;
    RCP<Epetra_CrsMatrix> Op2x_;
    RCP<Epetra_CrsMatrix> eye_;
    RCP<Epetra_Vector> u_, dx_, x_; /* dx_ should be du_ */
    RCP<Epetra_Vector> ExpVyVy_;
    RCP<Epetra_MultiVector> W_;
    bool isConverged_, backTracking_, test_, debug_;
    int NumStchFrcVec_;
    std::string solver_type;

    //Teuchos::ParameterList MeanParam_;
    double *t_, *dt_;
    double NormRHS_, NormRHStest_, toleranceRHS_;
    int iter_, maxNumIterations_, backTrack_, numBackTrackingSteps_;
    int MyPID;
};


#endif
