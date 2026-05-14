/*
 * =====================================================================================
 *
 *       Filename:  PDEAssembler.hpp
 *
 *    Description:  PDE operator assembly for the Burgers equation.
 *                  Extracted from Burger::Mean as part of Phase 5 decomposition.
 *                  Owns the FD stencil matrices (Laplacian, gradient, identity) and
 *                  assembles the RHS F(u), Jacobian J(u), and theta-method operator.
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
#ifndef BURGER_PDE_ASSEMBLER_HPP
#define BURGER_PDE_ASSEMBLER_HPP

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <string>

namespace Burger {

/**
 * @brief PDE operator assembly for the Burgers equation.
 *
 * Builds and manages the finite-difference stencil operators (Laplacian, gradient,
 * identity) and assembles the nonlinear RHS and Jacobian for the mean-field system.
 */
class PDEAssembler {
public:
    /**
     * @brief Construct and build FD operators on the given mesh.
     * @param m     Number of grid points.
     * @param mu    Viscosity parameter.
     * @param theta Theta-method parameter (0=explicit, 0.5=CN, 1=implicit).
     * @param Comm  Epetra communicator.
     */
    PDEAssembler(int m, double mu, double theta,
                 const Teuchos::RCP<Epetra_Comm>& Comm);

    /**
     * @brief Bilinear form: u3 = (Op1x * u1) .* (Op2x * u2).
     */
    void BilinearTerm(Teuchos::RCP<Epetra_Vector> u1,
                      Teuchos::RCP<Epetra_Vector> u2,
                      Teuchos::RCP<Epetra_Vector> u3);

    /**
     * @brief Assemble the PDE right-hand side F(u) scaled by dt.
     * @param u        Current solution.
     * @param ExpVyVy  Stochastic correction term.
     * @param dt       Time-step size.
     */
    void assembleRHS(const Teuchos::RCP<Epetra_Vector>& u,
                     const Teuchos::RCP<Epetra_Vector>& ExpVyVy,
                     double dt);

    /**
     * @brief Assemble the Jacobian J(u) and theta-method operator M - dt*θ*J.
     * @param u     Current solution.
     * @param dt    Time-step size.
     * @param theta Theta-method parameter.
     */
    void assembleJacobian(const Teuchos::RCP<Epetra_Vector>& u,
                          double dt, double theta);

    /// @brief Return an RCP to the current RHS vector.
    Teuchos::RCP<Epetra_Vector> getRHS() { return rhs_; }
    Epetra_Vector& getRHSRef() { return *rhs_; }

    /// @brief Return the assembled Jacobian J(u).
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return jac_; }

    /// @brief Return the theta-method operator (I - dt*θ*J).
    Teuchos::RCP<Epetra_CrsMatrix> getThetaJacobian() { return ThetaJac_; }

    /// @brief Return the mass matrix (identity for Burgers).
    Teuchos::RCP<Epetra_CrsMatrix> getMassMatrix() { return eye_; }

    /// @brief Return the grid coordinate vector.
    Teuchos::RCP<Epetra_Vector> getGridCoords() { return x_; }

    /// @brief Return the initial condition vector.
    Teuchos::RCP<Epetra_Vector> getInitialCondition() { return ICom_; }

private:
    void createLinOp(int m, double mu, const Teuchos::RCP<Epetra_Comm>& Comm);

    // FD stencil operators
    Teuchos::RCP<Epetra_CrsMatrix> LinOp_;   ///< Laplacian (diffusion)
    Teuchos::RCP<Epetra_CrsMatrix> Op1x_;    ///< Identity-like (for bilinear)
    Teuchos::RCP<Epetra_CrsMatrix> Op2x_;    ///< Gradient (advection)
    Teuchos::RCP<Epetra_CrsMatrix> eye_;     ///< Mass matrix (identity)

    // Jacobian workspace
    Teuchos::RCP<Epetra_CrsMatrix> LinJac_;
    Teuchos::RCP<Epetra_CrsMatrix> NlinJac_;
    Teuchos::RCP<Epetra_CrsMatrix> ThetaJac_;
    Teuchos::RCP<Epetra_CrsMatrix> jac_;
    Teuchos::RCP<Epetra_CrsMatrix> SpDiag_;
    Teuchos::RCP<Epetra_CrsMatrix> TmpMat1_;
    Teuchos::RCP<Epetra_CrsMatrix> TmpMat2_;

    // Work vectors
    Teuchos::RCP<Epetra_Vector> Op1x_u_;
    Teuchos::RCP<Epetra_Vector> Op2x_u_;
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> uu_;
    Teuchos::RCP<Epetra_Vector> tmp2_;
    Teuchos::RCP<Epetra_Vector> x_;     ///< Grid coordinates
    Teuchos::RCP<Epetra_Vector> ICom_;  ///< Initial condition

    bool debug_ = false;
};

} // namespace Burger

#endif // BURGER_PDE_ASSEMBLER_HPP
