/*
 * =====================================================================================
 *
 *       Filename:  PDEAssembler.hpp
 *
 *    Description:  PDE operator assembly for the barotropic Quasi-Geostrophic equation.
 *                  Delegates to the QG::QG Fortran/C++ wrapper for RHS, Jacobian,
 *                  and bilinear-form evaluation.
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
#ifndef QG_PDE_ASSEMBLER_HPP
#define QG_PDE_ASSEMBLER_HPP

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "QG.hpp"

namespace QG {

/**
 * @brief PDE operators for the barotropic Quasi-Geostrophic vorticity equation.
 *
 * Wraps the QG::QG solver (Fortran/C++ core) to provide a uniform assembly
 * interface for the RHS @f$ F(\psi) @f$, the Jacobian @f$ J(\psi) @f$, and
 * the bilinear advection form.  Also builds the theta-method operator
 * @f$ M - \Delta t\,\theta\,J @f$ used by QG::Mean's Newton solver.
 */
class PDEAssembler {
public:
    /**
     * @brief Construct the QG PDE assembler.
     * @param nx         Grid points in x.
     * @param ny         Grid points in y.
     * @param reynoldsNum Reynolds number.
     * @param topography  Topography parameter.
     * @param theta       Theta-method parameter.
     * @param Comm        Epetra communicator.
     */
    PDEAssembler(int nx, int ny, double reynoldsNum, double topography,
                 double theta, const Teuchos::RCP<Epetra_Comm>& Comm);

    /**
     * @brief Bilinear form via QG::QG::bilin().
     */
    void BilinearTerm(Teuchos::RCP<Epetra_Vector> u1,
                      Teuchos::RCP<Epetra_Vector> u2,
                      Teuchos::RCP<Epetra_Vector> u3);

    /**
     * @brief Assemble the PDE right-hand side F(u) scaled by dt.
     */
    void assembleRHS(const Teuchos::RCP<Epetra_Vector>& u,
                     const Teuchos::RCP<Epetra_Vector>& EVyVy,
                     double dt);

    /**
     * @brief Assemble the Jacobian and theta-method operator M - dt*θ*J.
     */
    void assembleJacobian(const Teuchos::RCP<Epetra_Vector>& u,
                          double dt, double theta);

    /**
     * @brief Compute the bare PDE right-hand side (without stochastic correction).
     */
    void getProblemRHS(Epetra_Vector& u, Epetra_Vector& F);

    Teuchos::RCP<Epetra_Vector> getRHS() { return rhs_; }
    Epetra_Vector& getRHSRef() { return *rhs_; }
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return jac_; }
    Teuchos::RCP<Epetra_CrsMatrix> getThetaJacobian() { return ThetaJac_; }
    Teuchos::RCP<Epetra_CrsMatrix> getMassMatrix() { return mass_; }

    /// @brief Return the diagonal mass vector (needed for ForcingProvider).
    Epetra_Vector& getDiagMass() { return *diagMass_; }

    /// @brief Return the QG object (needed for legacy access).
    Teuchos::RCP<::QG::QG> getQG() { return qg_; }

    /// @brief Total number of unknowns.
    int getN() const { return n_; }

private:
    int n_, nx_, ny_;
    Teuchos::RCP<::QG::QG> qg_;
    Teuchos::RCP<Epetra_CrsMatrix> mass_;
    Teuchos::RCP<Epetra_CrsMatrix> jac_;
    Teuchos::RCP<Epetra_CrsMatrix> ThetaJac_;
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> uu_;
    Teuchos::RCP<Epetra_Vector> diagMass_;
    bool debug_ = false;
};

} // namespace QG

#endif // QG_PDE_ASSEMBLER_HPP
