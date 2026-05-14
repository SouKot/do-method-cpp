/*
 * =====================================================================================
 *
 *       Filename:  LinearSolverWrapper.hpp
 *
 *    Description:  Solver-backend abstraction for the DO basis solver.
 *                  Encapsulates the Amesos / Amesos2 / Belos solver selection that
 *                  previously lived inline inside Problem_Interface / BasisSolver.
 *
 *                  Supported back-ends (selected at runtime via "Solver Package"
 *                  in StochasticParams.xml):
 *                    - "Amesos"  — direct, serial/parallel via Amesos factory
 *                    - "Amesos2" — direct, via Amesos2 factory
 *                    - anything else — iterative via Belos with Ifpack preconditioner
 *
 *                  Life-cycle:
 *                    1. Construct: symbolic factorisation only.
 *                    2. factorize(): numeric factorisation (call once per time step
 *                       after the matrix values change).
 *                    3. solve(LHS, RHS, label): solve A*LHS = RHS.
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
#ifndef DO_LINEAR_SOLVER_WRAPPER_HPP
#define DO_LINEAR_SOLVER_WRAPPER_HPP

#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Amesos2.hpp"
#include "Amesos2_Solver_decl.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVec.hpp"
#include "BelosSolverManager.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <string>

// =====================================================================================
/// @name Linear solver abstraction
// =====================================================================================
/// @{

/**
 * @brief Solver-backend wrapper used by the DO basis (V) solver.
 *
 * Hides Amesos / Amesos2 / Belos selection behind a uniform three-step API:
 *   - Constructor — symbolic factorisation (sparsity pattern analysis).
 *   - factorize() — numeric factorisation / preconditioner compute.
 *   - solve()     — back-substitution / iterative solve.
 *
 * The matrix is passed by RCP and modified in place by the caller (e.g.,
 * BasisSolver::computeBlocks rebuilds @c M - dt·A without allocating a new
 * object), so the solver always sees the current values when factorize() is
 * called.
 */
class LinearSolverWrapper {
public:
    typedef double                      ST;
    typedef Epetra_MultiVector          MV;
    typedef Epetra_Operator             OP;
    typedef Epetra_CrsMatrix            MAT;

    /**
     * @brief Construct the wrapper and perform symbolic factorisation.
     *
     * @param solverType  Back-end tag: "Amesos", "Amesos2", or Belos solver name.
     * @param solverName  Specific solver within the back-end (e.g. "Amesos_Klu",
     *                    "KLU2", "GMRES").
     * @param SolParams   Full parameter list for the chosen back-end, loaded from
     *                    @c <solverName>.xml.
     * @param matrix      The system matrix (held by RCP; modified in place by caller).
     * @param test        Enable timing output when true.
     * @param debug       Enable residual/diagnostic output when true.
     * @param dbgLvl      Bitmask controlling which diagnostics are printed.
     */
    LinearSolverWrapper(const std::string& solverType,
                        const std::string& solverName,
                        const Teuchos::RCP<Teuchos::ParameterList>& SolParams,
                        const Teuchos::RCP<Epetra_CrsMatrix>& matrix,
                        bool test,
                        bool debug,
                        int  dbgLvl);

    /**
     * @brief Numeric factorisation (Amesos/Amesos2) or preconditioner compute (Belos).
     *
     * Must be called once per time step — after the matrix values have been updated
     * by BasisSolver::computeBlocks() — and before any solve() calls.
     */
    void factorize();

    /**
     * @brief Solve the linear system @p matrix · LHS = RHS.
     *
     * Timing is reported when @p test_ is set; the residual norm is printed when
     * @p debug_ is set and DbgLvl_ % 2 == 0.
     *
     * @param LHS   On entry: initial guess (Belos) or ignored (direct).
     *              On exit: solution.
     * @param RHS   Right-hand side (unmodified).
     * @param label String printed in timing/residual output.
     * @return 0 on success.
     */
    int solve(Epetra_MultiVector& LHS,
              Epetra_MultiVector& RHS,
              const std::string&  label);

private:
    std::string                     solver_type_;
    bool                            test_;
    bool                            debug_;
    int                             dbgLvl_;
    Teuchos::RCP<Epetra_CrsMatrix>  matrix_;

    // Amesos (direct, serial)
    Teuchos::RCP<Epetra_LinearProblem>  v_prob_;
    Teuchos::RCP<Amesos_BaseSolver>     v_solve_;

    // Amesos2 (direct, parallel-capable)
    Teuchos::RCP<Amesos2::Solver<MAT, MV>> amesos2_solver_;

    // Belos (iterative) + Ifpack preconditioner
    Teuchos::RCP<Ifpack_Preconditioner>                  prec_;
    Teuchos::RCP<Belos::EpetraPrecOp>                    belosPrec_;
    Teuchos::RCP<Belos::LinearProblem<ST, MV, OP>>       problem_;
    Teuchos::RCP<Belos::SolverManager<ST, MV, OP>>       v_solve_iter_;
};

/// @}

#endif // DO_LINEAR_SOLVER_WRAPPER_HPP
