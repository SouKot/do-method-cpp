/*
 * =====================================================================================
 *
 *       Filename:  BasisSolver.hpp
 *
 *    Description:  DO stochastic basis (V) solver.
 *
 *                  BasisSolver manages the stochastic basis matrix V (size N×m),
 *                  its M-orthogonalisation, and the blocked linear system that
 *                  advances V one time step.  It replaces the former
 *                  Problem_Interface class; a backward-compatible alias
 *                  "using Problem_Interface = BasisSolver" is provided in
 *                  Problem_Interface.hpp.
 *
 *                  Key design changes vs. the old Problem_Interface:
 *                    - Solver back-end selection (Amesos / Amesos2 / Belos) is
 *                      delegated to LinearSolverWrapper, keeping BasisSolver free
 *                      of solver-specific ifdef guards.
 *                    - solver_ is held as std::optional<LinearSolverWrapper> so that
 *                      it can be emplaced in the constructor body after LHS_block_1_
 *                      is built — avoiding the initialiser-list ordering problem.
 *                    - syncJacobian() exposes an explicit API for updating A_ and
 *                      rebuilding LHS_block_1_ when the deterministic Jacobian
 *                      changes (e.g. after a mean-field Newton step).
 *
 *        Version:  1.0
 *        Created:  2026-05-13
 *       Revision:  none
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef DO_BASIS_SOLVER_HPP
#define DO_BASIS_SOLVER_HPP

#include "LinearSolverWrapper.hpp"
#include "StochasticState.hpp"
#include "DOUtils.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "BelosIMGSOrthoManager.hpp"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <optional>
#include <string>

extern "C" void dsytri_(char*, int*, double*, int*, int*, double*, int*);

// =====================================================================================
/// @name DO basis solver
// =====================================================================================
/// @{

/**
 * @brief Stochastic basis (V) solver for the Dynamical Orthogonality method.
 *
 * Manages the N×m basis matrix V and advances it one DO time step via the
 * bordered linear system:
 * @f[
 *   (M - \Delta t J)\, X = F_v
 * @f]
 * followed by a bordered solve for the updating direction.
 *
 * The linear solver back-end (Amesos / Amesos2 / Belos) is selected at
 * construction time from the "Solver Package" entry of @p SolverParams and
 * encapsulated in a LinearSolverWrapper member.
 *
 * @see LinearSolverWrapper
 */
class BasisSolver {
public:
    /**
     * @brief Construct the basis solver.
     *
     * Reads solver configuration, allocates workspace, builds the initial
     * @c M - dt·A matrix, and performs the symbolic factorisation.
     *
     * @param A           Deterministic Jacobian (shared; updated externally).
     * @param SolverParams Parameter list; must contain "Solver Package",
     *                     "StochBasisFile", "Class Testing", etc.
     * @param udet        Mean-field solution vector.
     * @param m           Number of stochastic basis modes.
     * @param t           Pointer to current simulation time (read only).
     * @param dt          Initial time-step size.
     * @param W           Stochastic forcing basis (N × nW).
     * @param mass        Mass matrix M.
     * @param iter        Initial stochastic iteration counter.
     */
    BasisSolver(const Teuchos::RCP<Epetra_CrsMatrix>& A,
                const Teuchos::RCP<Teuchos::ParameterList>& SolverParams,
                const Teuchos::RCP<Epetra_Vector>& udet,
                int& m,
                double* t,
                double& dt,
                const Teuchos::RCP<Epetra_MultiVector>& W,
                const Teuchos::RCP<Epetra_CrsMatrix>& mass,
                int iter,
                const Teuchos::RCP<StochasticState>& sharedState);

    ~BasisSolver();

    // -----------------------------------------------------------------------
    /// @name NOX/LOCA compatibility stubs
    // -----------------------------------------------------------------------
    /// @{
    bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);
    bool computeShiftedMatrix(double alpha, double beta, const Epetra_Vector& x, Epetra_Operator& A);
    /// @}

    // -----------------------------------------------------------------------
    /// @name Basis management
    // -----------------------------------------------------------------------
    /// @{

    /// M-orthogonalise V in place (Belos::IMGSOrthoManager).
    void MOrth();

    /**
     * @brief Numerically factorise and solve the bordered V sub-system.
     *
     * @param Vold Basis at the previous time step.
     * @return EXIT_SUCCESS.
     */
    int v_stoch_init(Epetra_MultiVector* Vold);

    /**
     * @brief Rebuild @c M-dt·A and the RHS forcing @c F_v.
     *
     * Must be called once per time step, before v_stoch_init().
     *
     * @param dt Current time-step size.
     */
    void computeBlocks(double dt);

    /// Transfer the QR factor R into the stochastic coefficients y.
    void TransferNorm();

    /**
     * @brief Update the stored Jacobian and rebuild @c M-dt·A in place.
     *
     * Call after the deterministic solver updates A (e.g. after a Newton step)
     * so the next computeBlocks() / v_stoch_init() cycle uses the new operator.
     *
     * @param newA New deterministic Jacobian; values are copied into A_ in place.
     */
    void syncJacobian(const Epetra_CrsMatrix& newA);

    /// Solve A·LHS = RHS via the LinearSolverWrapper.
    int SolveV(Epetra_MultiVector& LHS, Epetra_MultiVector& RHS, const std::string& label);

    /// Initialise V from a file or random vectors.
    void init_v(double t);

    // Legacy stubs (not yet implemented)
    void computeBlocksold(double dt);
    void computeExpectations(double dt);
    void PostProcess();

    /// @}

    // -----------------------------------------------------------------------
    /// @name Accessors
    // -----------------------------------------------------------------------
    /// @{

    Epetra_Vector& getSolution()    { return *(sharedState_->V->operator()(0)); }
    Epetra_Vector& getMean()        { return getSolution(); }
    Epetra_CrsMatrix& getJacobian() { return *LHS_block_1_; }
    Teuchos::RCP<Epetra_MultiVector> getEyy() { return Eyy; }
    Teuchos::RCP<Epetra_MultiVector> get_y()     { return sharedState_->y; }

    /// Return the stochastic basis matrix V (RCP; caller shares ownership).
    Teuchos::RCP<Epetra_MultiVector> getBasisV() const { return sharedState_->V; }

    void   set_y(const Epetra_MultiVector& yy) { *(sharedState_->y) = yy; }
    double& get_frcStrenth()                   { return stchFrcStren_; }
    void   set_frcStrength(double v)           { stchFrcStren_ = v; }

    /// @}

    // -----------------------------------------------------------------------
    /// @name Public data (shared with Y_Stoch / driver)
    // -----------------------------------------------------------------------

    double dt_, stchFrcStren_;
    int    m_;
    Teuchos::RCP<Teuchos::ParameterList> SolverParams_;
    std::string                          InitFile_;

    Teuchos::RCP<Epetra_CrsMatrix>  A_;
    Teuchos::RCP<Epetra_CrsMatrix>  LHS_block_1_;
    Teuchos::RCP<Epetra_CrsMatrix>  mass_;
    Teuchos::RCP<Epetra_Vector>     udet_;

    int    stochiter, MyPID;
    double *t_, Tol_;

    Teuchos::RCP<Epetra_MultiVector> W_;
    Teuchos::RCP<Epetra_MultiVector> V1;
    Teuchos::RCP<Epetra_MultiVector> eye;
    Teuchos::RCP<Epetra_MultiVector> Ezy;
    Teuchos::RCP<Epetra_MultiVector> Eyy, Rvec;
    Teuchos::RCP<Epetra_MultiVector> Eyyy;
    Teuchos::RCP<Epetra_MultiVector> EVyVyBasis;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>> R;
    Teuchos::RCP<Epetra_MultiVector> EyyInv;
    Teuchos::RCP<Epetra_Map>         map_, y_map;
    Teuchos::RCP<Epetra_MultiVector> z, RHS_block_1;
    Teuchos::RCP<Epetra_LocalMap>    locmap_;

    // Scalar typedefs kept for callers that reference them via BasisSolver::ST etc.
    typedef double                  ST;
    typedef Teuchos::ScalarTraits<ST> SCT;
    typedef SCT::magnitudeType      MT;
    typedef Epetra_MultiVector      MV;
    typedef Epetra_Operator         OP;
    typedef Epetra_CrsMatrix        MAT;

    int  iteration, MaxIter_, DbgLvl_;
    bool test_, debug_;

private:
    /// Shared state with Y_Stoch — replaces former sync methods.
    Teuchos::RCP<StochasticState> sharedState_;

    // Emplaced in the constructor body after LHS_block_1_ is built.
    std::optional<LinearSolverWrapper> solver_;
    std::string                        solver_type_;
};

/// @}

#endif // DO_BASIS_SOLVER_HPP
