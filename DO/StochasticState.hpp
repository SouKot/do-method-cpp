/*
 * =====================================================================================
 *
 *       Filename:  StochasticState.hpp
 *
 *    Description:  Shared state struct for coupling BasisSolver (V-solver) and
 *                  CoeffSolver (coefficient solver) without circular references.
 *
 *                  Both classes are constructed with an RCP<StochasticState> and
 *                  access shared data exclusively through it, making the coupling
 *                  explicit and the initialisation order enforceable.
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
#ifndef DO_STOCHASTIC_STATE_HPP
#define DO_STOCHASTIC_STATE_HPP

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include <functional>

// =====================================================================================
/// @name Shared stochastic state
// =====================================================================================
/// @{

/**
 * @brief Shared data contract between CoeffSolver and BasisSolver.
 *
 * Instead of exchanging raw pointers or calling sync methods, both classes
 * hold an @c RCP<StochasticState> and read/write the fields they own:
 *
 * | Field      | Written by     | Read by        |
 * |------------|----------------|----------------|
 * | V          | BasisSolver    | CoeffSolver        |
 * | y          | CoeffSolver        | BasisSolver    |
 * | EDEyy  | CoeffSolver        | BasisSolver    |
 * | uMean       | Mean solver    | both           |
 * | A          | Mean solver    | both           |
 * | W          | Mean solver    | both           |
 */
struct StochasticState {
    /// Stochastic basis matrix (N × m). Owned by BasisSolver.
    Teuchos::RCP<Epetra_MultiVector> V;

    /// Stochastic coefficient matrix (m × S). Owned by CoeffSolver.
    Teuchos::RCP<Epetra_MultiVector> y;

    /// E[ d(ExpDExp) / yy ] coupling term (N × m). Computed by CoeffSolver.
    Teuchos::RCP<Epetra_MultiVector> EDEyy;

    /// Mean-field solution vector. Owned by the Mean solver.
    Teuchos::RCP<Epetra_Vector> uMean;

    /// Deterministic Jacobian. Owned by the Mean solver.
    Teuchos::RCP<Epetra_CrsMatrix> A;

    /// Stochastic forcing basis (N × nW). Owned by the Mean solver.
    Teuchos::RCP<Epetra_MultiVector> W;

    // -----------------------------------------------------------------------
    /// @name Bilinear-term callback
    // -----------------------------------------------------------------------
    /// @{

    /**
     * @brief Callback for the model-specific bilinear form u⊗v → uv.
     *
     * Replaces the former @c model_->BilinearTerm() coupling between CoeffSolver
     * and Mean.  The driver sets this once after constructing the Mean object.
     *
     * Signature:  void(RCP<Epetra_Vector> u, RCP<Epetra_Vector> v,
     *                   RCP<Epetra_Vector> uv)
     *
     * For the LOCA path (need_locaInterface == 1), this is unused — the
     * Fortran @c model_bil() call is used directly.
     */
    using BilinearFn = std::function<void(const Teuchos::RCP<Epetra_Vector>&,
                                          const Teuchos::RCP<Epetra_Vector>&,
                                          const Teuchos::RCP<Epetra_Vector>&)>;
    BilinearFn bilinearTerm;

    /// @}
};

/// @}

#endif // DO_STOCHASTIC_STATE_HPP
