/*
 * =====================================================================================
 *
 *       Filename:  DOTimeLoop.hpp
 *
 *    Description:  Inline time-stepping helpers shared by the Burgers/QG driver
 *                  (timedependent_do.cpp) and the SWE driver (timedependent_swe.cpp).
 *                  Contains no model-type knowledge; all helpers are compiled into
 *                  every DO target via inclusion.
 *
 *        Version:  1.0
 *        Created:  2026-05-06
 *       Revision:  none
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef DO_TIME_LOOP_HPP
#define DO_TIME_LOOP_HPP

#include "DOUtils.hpp"
#include "StochasticState.hpp"
#include "StochIO.hpp"
#include "Problem_Interface.hpp"
#include "CoeffSolver.hpp"
#include "EpetraExt_MatrixMatrix.h"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StrUtils.hpp"
#include <string>

// =====================================================================================
/// @name Data structures
// =====================================================================================
/// @{

/**
 * @brief Teuchos timer handles for the three stages of a stochastic update step.
 *
 * Passed as a bundle into runStochStep() so that timing instrumentation can be
 * enabled or disabled in one place without changing the helper's signature.
 */
struct StochTimers {
    Teuchos::RCP<Teuchos::Time> coeff; ///< Timer for the Y (coefficient) solve.
    Teuchos::RCP<Teuchos::Time> bilin; ///< Timer for bilinear-form / expectation computation.
    Teuchos::RCP<Teuchos::Time> basis; ///< Timer for the V (basis) update.
};

/// @}

// =====================================================================================
/// @name Initialisation helpers
// =====================================================================================
/// @{

/**
 * @brief Write the initial mass and Jacobian matrices to disk, then sync detA <- A.
 *
 * Writes @c mass.mm and @c jac.mm for offline diagnostics in MatrixMarket format.
 * Also copies @p A into @p detA so that the V-solver's Jacobian is consistent with
 * the assembled operator at the very start of time stepping.
 *
 * @param[in]  A       Current assembled Jacobian (RCP to avoid an extra copy).
 * @param[out] detA    V-solver's local Jacobian copy; overwritten with @p A.
 * @param[in]  massmat Mass matrix written to @c mass.mm.
 */
inline void dumpInitialJacobian(const Teuchos::RCP<Epetra_CrsMatrix>& A,
                                Epetra_CrsMatrix &detA,
                                Epetra_CrsMatrix &massmat) {
    StochIO::writeMatrix("mass.mm", massmat);
    StochIO::writeMatrix("jac.mm", detA);
    EpetraExt::MatrixMatrix::Add(*A, false, 1.0, detA, 0.0);
}

/// @}

// =====================================================================================
/// @name Per-step stochastic update
// =====================================================================================
/// @{

/**
 * @brief Advance the stochastic (Y, V) unknowns for one time step.
 *
 * Performs the three-stage stochastic sub-step:
 *   -# **Coefficient update** — CoeffSolver::StochasticIterations() updates the coefficient matrix Y.
 *   -# **Basis update** — Problem_Interface::computeBlocks(), v_stoch_init(), TransferNorm()
 *      update the stochastic basis V.
 *   -# **Statistics update** — HBilinV(), computeEyyTyT(), computeEVyVy() update the
 *      second-moment statistics needed by the mean solver.
 *
 * The caller is responsible for refreshing the stochastic forcing (via
 * @c model->refreshForcing(t)) **before** this call; CoeffSolver and Problem_Interface
 * already hold a shared pointer to Wbase and see the update automatically.
 *
 * @param[in]     isStochOn    If false the function returns immediately (deterministic run).
 * @param[in]     timeProf     Enable Teuchos timer instrumentation when true.
 * @param[in]     dt           Current time-step size.
 * @param[in,out] y_interface  Stochastic coefficient solver.
 * @param[in,out] Vstoch       Stochastic basis / V-solver.
 * @param[in,out] Vn           Current stochastic basis vectors (updated in place).
 * @param[in]     tmr          Timer bundle; fields used only when @p timeProf is true.
 */
inline void runStochStep(bool isStochOn,
                         bool timeProf,
                         double dt,
                         const Teuchos::RCP<CoeffSolver>& y_interface,
                         const Teuchos::RCP<Problem_Interface>& Vstoch,
                         const Teuchos::RCP<Epetra_MultiVector>& Vn,
                         const StochTimers &tmr) {
    if (!isStochOn) return;

    if (timeProf) tmr.coeff->start();
    y_interface->StochasticIterations();
    if (timeProf) {
        tmr.coeff->stop();
        tmr.coeff->incrementNumCalls();
    }

    Epetra_MultiVector Vtemp(*Vn);
    if (timeProf) tmr.basis->start();
    Vstoch->computeBlocks(dt);
    Vstoch->v_stoch_init(&Vtemp);
    Vstoch->TransferNorm();
    if (timeProf) {
        tmr.basis->stop();
        tmr.basis->incrementNumCalls();
    }

    if (timeProf) tmr.bilin->start();
    y_interface->HBilinV();
    y_interface->computeEyyTyT();
    y_interface->computeEVyVy();
    if (timeProf) {
        tmr.bilin->stop();
        tmr.bilin->incrementNumCalls();
    }
}

/// @}

// =====================================================================================
/// @name Periodic output
// =====================================================================================
/// @{

/**
 * @brief Write periodic output snapshots during time integration.
 *
 * Two independent output actions are performed at each call:
 *   - **Eyy ring buffer**: the @p expyy snapshot is appended to @p Eyy.  Once the
 *     buffer holds @p nv snapshots it is flushed to @c tsEyy_<t>.mm and reset.
 *   - **Snapshot files**: when the simulation time crosses the next print interval
 *     the current V, yT, and mean-solution vectors are written to MatrixMarket files.
 *
 * File format is controlled at compile time by the @c need_locaInterface macro:
 *   - @c 1 (SWEDO): HYMLS::MatrixUtils::mmwrite
 *   - @c 0 (BRGRDO/QGDO): EpetraExt::MultiVectorToMatrixMarketFile
 *
 * @param[in]     t            Current simulation time.
 * @param[in,out] eyyBufferIdx Current write position in the Eyy ring buffer (0-based).
 * @param[in]     nv           Eyy ring-buffer capacity (flush threshold).
 * @param[in,out] Eyy          Ring-buffer multi-vector accumulating Eyy snapshots.
 * @param[in]     expyy        Current E[yy^T] snapshot to append.
 * @param[in]     isStochOn    Suppress V/yT output when false.
 * @param[in]     Vn           Current stochastic basis.
 * @param[in]     yTrans       Current stochastic coefficient matrix.
 * @param[in]     soln         Current mean solution.
 * @param[in]     prntintvl    Snapshot print interval in simulation-time units.
 * @param[in,out] count        Time of the most recently written snapshot; updated on write.
 */
inline void saveTimestepOutputs(double t,
                                int &eyyBufferIdx,
                                int nv,
                                Epetra_MultiVector &Eyy,
                                Epetra_MultiVector &expyy,
                                bool isStochOn,
                                Epetra_MultiVector &Vn,
                                Epetra_MultiVector &yTrans,
                                Epetra_Vector &soln,
                                double prntintvl,
                                double &count) {
    for (int i = 0; i < expyy.NumVectors(); i++)
        for (int j = 0; j < expyy.MyLength(); j++)
            Eyy.ReplaceGlobalValue(i * expyy.MyLength() + j, eyyBufferIdx, expyy[i][j]);
    eyyBufferIdx++;
    if (eyyBufferIdx == nv && isStochOn) {
        std::string fname = "tsEyy_" + Teuchos::toString(float(t)) + ".mm";
        HYMLS::MatrixUtils::mmwrite(fname, Eyy);
        eyyBufferIdx = 0;
    }

    if ((t - count - prntintvl) >= 0.0) {
        std::string fname;
        if (isStochOn) {
            fname = "v_" + Teuchos::toString(float(t)) + ".mm";
#if need_locaInterface == 1
            HYMLS::MatrixUtils::mmwrite(fname, Vn);
#else
            StochIO::writeMV(fname, Vn);
#endif
            fname = "yT_" + Teuchos::toString(float(t)) + ".mm";
#if need_locaInterface == 1
            HYMLS::MatrixUtils::mmwrite(fname, yTrans);
#else
            StochIO::writeMV(fname, yTrans);
#endif
        }
        fname = "mean_" + Teuchos::toString(float(t)) + ".mm";
#if need_locaInterface == 1
        HYMLS::MatrixUtils::mmwrite(fname, soln);
#else
        StochIO::writeMV(fname, soln);
#endif
        count += prntintvl;
    }
}

/// @}

#endif // DO_TIME_LOOP_HPP
