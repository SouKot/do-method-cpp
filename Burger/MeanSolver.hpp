/*
 * =====================================================================================
 *
 *       Filename:  Mean.hpp
 *
 *    Description:  Mean-field solver for the stochastic 1-D viscous Burgers equation
 *                  u_t + u*u_x = mu*u_xx + sigma*dW.  Coordinates PDE assembly
 *                  (PDEAssembler), time-dependent stochastic forcing
 *                  (ForcingProvider), and linear solves (LinearSolverWrapper)
 *                  behind an implicit theta-method Newton iteration.
 *
 *        Version:  2.0
 *        Created:  04/15/2018 03:38:12 PM
 *       Revision:  2026-05-13
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef MEAN_SOLVER_HPP
#define MEAN_SOLVER_HPP

#include "PDEAssembler.hpp"
#include "ForcingProvider.hpp"
#include "../DO/LinearSolverWrapper.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <cmath>
#include <optional>
#include <string>

/**
 * @brief Mean-field solver for the stochastic 1-D viscous Burgers equation.
 *
 * Solves the deterministic component @f$ \bar{u}(x,t) @f$ of the DO expansion
 * @f$ u = \bar{u} + V Y @f$ using an implicit theta-method with Newton iteration.
 * PDE assembly is delegated to Burger::PDEAssembler, stochastic forcing to
 * Burger::ForcingProvider, and the linear solver back-end to LinearSolverWrapper.
 */
class MeanSolver {
public:
    MeanSolver(Teuchos::RCP<Teuchos::ParameterList> PrmLst,
         double* t, double* dt,
         Teuchos::RCP<Epetra_Comm> Comm);

    // -- PDE assembly (delegated) --
    void BilinearTerm(Teuchos::RCP<Epetra_Vector> u1,
                      Teuchos::RCP<Epetra_Vector> u2,
                      Teuchos::RCP<Epetra_Vector> u3)
    { pde_.BilinearTerm(u1, u2, u3); }

    /// @brief Evaluate the PDE right-hand side F(u) including stochastic correction.
    void computeF(Epetra_Vector& u, Epetra_Vector& F) {
        pde_.assembleRHS(Teuchos::rcpFromRef(u), EVyVy_, *dt_);
        F = pde_.getRHSRef();
    }

    Epetra_Vector& getF() { return pde_.getRHSRef(); }

    /// @brief Assemble the Jacobian J(u) and the theta-method operator I − dt·θ·J.
    void computeJacobian(Epetra_Vector& x, Epetra_CrsMatrix& A) {
        pde_.assembleJacobian(Teuchos::rcpFromRef(x), *dt_, theta_);
        A = *pde_.getJacobian();
    }

    Teuchos::RCP<Epetra_CrsMatrix> getJacobian() { return pde_.getJacobian(); }
    Teuchos::RCP<Epetra_CrsMatrix> getMassMatrix() { return pde_.getMassMatrix(); }

    // -- Forcing (delegated) --
    void refreshForcing(double t) { forcing_.refreshForcing(t); }
    int getDimW() { return forcing_.getDimW(); }
    Teuchos::RCP<Epetra_MultiVector> getW() { return forcing_.getW(); }

    // -- Time stepping --
    /// @brief Run Newton iteration to advance the mean field one implicit time step.
    bool NewtonSolver();

    /// @brief Write the mean solution to a text file.
    void WriteSolution(std::string filename, double param,
                       const Epetra_Vector& soln);

    // -- Accessors --
    Teuchos::RCP<Epetra_Vector> getSolution() { return u_; }
    void setEVyVy(Teuchos::RCP<Epetra_Vector> EVyVy) { EVyVy_ = EVyVy; }

private:
    // -----------------------------------------------------------------------
    /// @brief Selects which convergence criterion governs the Newton loop.
    ///
    /// Set via "Newton Convergence Criterion" in params.xml:
    ///  - "Absolute"          : ||F_k|| < tol_desired
    ///  - "Relative Fixed"    : ||F_k||/||F_0|| < tol_relative
    ///  - "Relative Adaptive" : ||F_k||/||F_0|| < alpha * dt^(p-1)
    ///  - "Eisenstat-Walker"  : eta_k = gamma*(||F_k||/||F_{k-1}||)^ew_alpha
    // -----------------------------------------------------------------------
    enum class NewtonCriterion { Absolute, RelativeFixed, RelativeAdaptive, EisenstatWalker };
    static NewtonCriterion parseCriterion(const std::string& s);
    NewtonCriterion criterion_;

    void ThetaStepper();
    void runBackTracking();
    int  LinSolve(Epetra_Vector& LHS, Epetra_Vector& RHS);
    /// Evaluate the active convergence criterion; returns true when converged.
    bool checkConvergence(double normFk, double normF0,
                          double normFprev, double etaEW) const;

    // Composed subsystems
    Burger::PDEAssembler     pde_;
    Burger::ForcingProvider  forcing_;
    std::optional<LinearSolverWrapper> solver_;

    // Solution state
    Teuchos::RCP<Epetra_Vector> u_;
    Teuchos::RCP<Epetra_Vector> u0_;
    Teuchos::RCP<Epetra_Vector> dx_;
    Teuchos::RCP<Epetra_Vector> ThetaRHS_;
    Teuchos::RCP<Epetra_Vector> EVyVy_;

    // Parameters
    double* t_;
    double* dt_;
    double  theta_;

    // Newton parameters
    double toleranceRHS_;   ///< Absolute desired tolerance (silent success if met)
    double tolAcceptable_;  ///< Absolute acceptable tolerance (warn but continue)
    double tolRelative_;    ///< Fixed relative tolerance  [Relative Fixed criterion]
    double ewGamma_;        ///< Eisenstat-Walker gamma constant (default 0.9)
    double ewAlphaEW_;      ///< Eisenstat-Walker exponent  (default 1.618)
    double alphaAdapt_;     ///< Safety factor for Relative Adaptive: tol = alpha*dt^(p-1)
    int    schemeOrder_;    ///< Temporal order p of the theta-scheme (1=BE, 2=CN)
    double NormRHS_, NormRHStest_;
    int    iter_, maxNumIterations_;
    int    backTrack_, numBackTrackingSteps_;
    bool   isConverged_, backTracking_;
    bool   test_, debug_;
    int    MyPID_;
};

#endif
