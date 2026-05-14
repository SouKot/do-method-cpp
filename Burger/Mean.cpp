/*
 * =====================================================================================
 *
 *       Filename:  Mean.cpp
 *
 *    Description:  Implementation of the Burger::Mean thin coordinator.
 *                  Constructor, Newton solver, theta-method, and I/O.
 *                  PDE assembly delegated to PDEAssembler, solver to
 *                  LinearSolverWrapper.
 *
 *        Version:  2.0  (Phase 5 decomposition)
 *        Created:  04/14/2018 07:55:56 AM
 *       Revision:  2026-05-13  Phase 5 SRP decomposition
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#include "Mean.hpp"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Teuchos;

// =====================================================================================
Mean::Mean(RCP<Teuchos::ParameterList> PrmLst, double* t, double* dt,
           RCP<Epetra_Comm> Comm)
    : pde_(PrmLst->get("nx", 200),
           PrmLst->get("mu", 0.1),
           PrmLst->get("theta", 0.5),
           Comm)
    , forcing_(PrmLst->get("No. of vectors in stoch. forcing", 2),
               pde_.getGridCoords())
    , t_(t)
    , dt_(dt)
    , theta_(PrmLst->get("theta", 0.5))
    , toleranceRHS_(10e-12)
    , maxNumIterations_(30)
    , numBackTrackingSteps_(5)
    , backTracking_(true)
    , isConverged_(false)
    , test_(PrmLst->get("Testing", true))
    , debug_(PrmLst->get("Debugging", true))
{
    // Allocate solution vectors
    Epetra_Map Map(PrmLst->get("nx", 200), 0, *Comm);
    u_        = rcp(new Epetra_Vector(Map));
    u0_       = rcp(new Epetra_Vector(Map));
    dx_       = rcp(new Epetra_Vector(Map));
    ThetaRHS_ = rcp(new Epetra_Vector(Map));
    ExpVyVy_  = rcp(new Epetra_Vector(Map));
    MyPID_    = Comm->MyPID();

    // Set initial condition from PDEAssembler
    *u_  = *pde_.getInitialCondition();
    *u0_ = *u_;

    // Build initial forcing
    forcing_.createW();

    // Setup solver
    std::string solver_type = PrmLst->get("Solver Package", "Amesos2");
    Teuchos::ParameterList& SolverSublist = PrmLst->sublist(solver_type);
    std::string solver_name = SolverSublist.get("Solver name", "Amesos_Scalapack");
    int numProc = Comm->NumProc();

    if (numProc == 1 && solver_type == "Amesos")
        solver_name = "Amesos_Klu";
    else if (numProc == 1 && solver_type == "Amesos2")
        solver_name = "KLU2";

    const RCP<Teuchos::ParameterList> SolParams =
        Teuchos::getParametersFromXmlFile(solver_name + ".xml");
    if (MyPID_ == 0)
        std::cout << "================================================\n"
                  << solver_type << " solver(" << solver_name
                  << ") has been chosen as mean solver\n"
                  << "================================================\n";

    // Build initial ThetaJac so LinearSolverWrapper has a filled matrix
    pde_.assembleJacobian(u_, *dt_, theta_);

    solver_.emplace(solver_type, solver_name, SolParams,
                    pde_.getThetaJacobian(), test_, debug_, 0);
}

// =====================================================================================
void Mean::ThetaStepper()
{
    pde_.assembleRHS(u_, ExpVyVy_, *dt_);
    ThetaRHS_->Update(theta_, pde_.getRHSRef(), 0.0);

    pde_.assembleRHS(u0_, ExpVyVy_, *dt_);
    ThetaRHS_->Update(1 - theta_, pde_.getRHSRef(), 1.0);
    ThetaRHS_->Update(1.0, *u_, -1.0, *u0_, -1.0);
}

// =====================================================================================
int Mean::LinSolve(Epetra_Vector& LHS, Epetra_Vector& RHS)
{
    solver_->factorize();
    return solver_->solve(LHS, RHS, "Mean solve");
}

// =====================================================================================
bool Mean::NewtonSolver()
{
    isConverged_ = false;
    dx_->PutScalar(0.0);

    ThetaStepper();
    ThetaRHS_->Scale(-1.0);
    ThetaRHS_->Norm2(&NormRHS_);

    if (debug_)
        std::cout << "\nMean:: Newton:      initial norm: " << NormRHS_;

    for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
        pde_.assembleJacobian(u_, *dt_, theta_);
        LinSolve(*dx_, *ThetaRHS_);
        u_->Update(1.0, *dx_, 1.0);

        ThetaStepper();
        ThetaRHS_->Scale(-1.0);
        ThetaRHS_->Norm2(&NormRHStest_);

        if (debug_) {
            std::cout << "\nNewton:      iter: " << iter_;
            std::cout << "\nNewton:      norm: " << NormRHStest_;
        }

        if (NormRHStest_ < toleranceRHS_) {
            if (debug_) std::cout << "\nSuccess...";
            break;
        }
        if (backTracking_ && (NormRHS_ < NormRHStest_))
            RunBackTracking();
        NormRHS_ = NormRHStest_;
    }

    if (iter_ == maxNumIterations_) {
        std::cout << "\nNewton: ---> TROUBLE" << __FILE__ << __LINE__;
        std::cout << "\nNewton not converged after Max number of iter!!!!!!\n";
        std::cout << "\nrhs norm = " << NormRHStest_;
    } else {
        *u0_ = *u_;
        isConverged_ = true;
    }
    return isConverged_;
}

// =====================================================================================
void Mean::RunBackTracking()
{
    double reduction = -1.0 / 2;
    for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_) {
        if (NormRHStest_ < NormRHS_) {
            if (debug_) std::cout << "\nSuccess...";
            break;
        }
        u_->Update(reduction, *dx_, 1.0);
        ThetaStepper();
        ThetaRHS_->Norm2(&NormRHStest_);
        reduction /= 2.0;
    }
    if (backTrack_ == numBackTrackingSteps_)
        std::cout << "\nNewton: --> BACKTRACKING FAILED" << __FILE__ << __LINE__;
}

// =====================================================================================
void Mean::WriteSolution(std::string filename, double param,
                         const Epetra_Vector& soln)
{
    Teuchos::RCP<std::ostream> out;
    if (soln.Comm().MyPID() == 0)
        out = Teuchos::rcp(new std::ofstream(filename.c_str()));
    else
        out = Teuchos::rcp(new Teuchos::oblackholestream());
    (*out) << std::setw(15) << std::setprecision(15);
    out->setf(std::ios::scientific);
    (*out) << param;
    (*out) << *(HYMLS::MatrixUtils::Gather(soln, 0));
}
