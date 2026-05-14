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
    int get_dim_W() { return forcing_.get_dim_W(); }
    Teuchos::RCP<Epetra_MultiVector> get_W() { return forcing_.get_W(); }

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
    void ThetaStepper();
    void RunBackTracking();
    int  LinSolve(Epetra_Vector& LHS, Epetra_Vector& RHS);

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
    double toleranceRHS_, NormRHS_, NormRHStest_;
    int    iter_, maxNumIterations_;
    int    backTrack_, numBackTrackingSteps_;
    bool   isConverged_, backTracking_;
    bool   test_, debug_;
    int    MyPID_;
};

#endif
