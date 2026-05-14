/*
 * =====================================================================================
 *
 *       Filename:  Mean.hpp
 *
 *    Description:  Mean-field solver for the stochastic Burgers equation.
 *                  Thin coordinator that delegates to:
 *                    - PDEAssembler:       RHS and Jacobian assembly
 *                    - ForcingProvider:    time-dependent stochastic forcing W
 *                    - LinearSolverWrapper: linear solve (Amesos/Amesos2/Belos)
 *
 *        Version:  2.0  (Phase 5 decomposition)
 *        Created:  04/15/2018 03:38:12 PM
 *       Revision:  2026-05-13  Phase 5 SRP decomposition
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#ifndef mean_solve_h
#define mean_solve_h

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
 * @brief Mean-field solver for the stochastic Burgers equation.
 *
 * Thin coordinator: PDE assembly lives in Burger::PDEAssembler,
 * stochastic forcing in Burger::ForcingProvider, solver back-end in
 * LinearSolverWrapper.
 */
class Mean {
public:
    Mean(Teuchos::RCP<Teuchos::ParameterList> PrmLst,
         double* t, double* dt,
         Teuchos::RCP<Epetra_Comm> Comm);

    // -- PDE assembly (delegated) --
    void BilinearTerm(Teuchos::RCP<Epetra_Vector> u1,
                      Teuchos::RCP<Epetra_Vector> u2,
                      Teuchos::RCP<Epetra_Vector> u3)
    { pde_.BilinearTerm(u1, u2, u3); }

    void computeF(Epetra_Vector& u, Epetra_Vector& F) {
        pde_.assembleRHS(Teuchos::rcpFromRef(u), ExpVyVy_, *dt_);
        F = pde_.getRHSRef();
    }

    Epetra_Vector& getF() { return pde_.getRHSRef(); }

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
    bool NewtonSolver();
    void WriteSolution(std::string filename, double param,
                       const Epetra_Vector& soln);

    // -- Accessors --
    Teuchos::RCP<Epetra_Vector> getSolution() { return u_; }
    void setExpVyVy(Teuchos::RCP<Epetra_Vector> ExpVyVy) { ExpVyVy_ = ExpVyVy; }

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
    Teuchos::RCP<Epetra_Vector> ExpVyVy_;

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
