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
 *        Version:  2.0
 *        Created:  04/14/2018 07:55:56 AM
 *       Revision:  2026-05-13
 *       Compiler:  gcc / clang (C++17)
 *
 *         Author:  Sourabh Kotnala (), sauravkotnala@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */
#include "MeanSolver.hpp"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Teuchos;

// =====================================================================================
MeanSolver::MeanSolver(RCP<Teuchos::ParameterList> PrmLst, double* t, double* dt,
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
    , toleranceRHS_  (PrmLst->get("Newton Tolerance",             1.0e-10))
    , tolAcceptable_ (PrmLst->get("Newton Acceptable Tolerance",   1.0e-5))
    , tolRelative_   (PrmLst->get("Newton Relative Tolerance",     1.0e-4))
    , ewGamma_       (PrmLst->get("Newton EW Gamma",               0.9))
    , ewAlphaEW_     (PrmLst->get("Newton EW Alpha",               1.618))
    , alphaAdapt_    (PrmLst->get("Newton Relative Tolerance Factor", 0.01))
    , schemeOrder_   (PrmLst->get("Newton Scheme Order",           2))
    , maxNumIterations_(PrmLst->get("Newton Max Iter",             50))
    , numBackTrackingSteps_(PrmLst->get("Newton Backtrack Steps",  20))
    , backTracking_(true)
    , isConverged_(false)
    , criterion_(parseCriterion(
          PrmLst->get("Newton Convergence Criterion", std::string("Relative Fixed"))))
    , test_(PrmLst->get("Testing", true))
    , debug_(PrmLst->get("Debugging", true))
{
    // Allocate solution vectors
    Epetra_Map Map(PrmLst->get("nx", 200), 0, *Comm);
    u_        = rcp(new Epetra_Vector(Map));
    u0_       = rcp(new Epetra_Vector(Map));
    dx_       = rcp(new Epetra_Vector(Map));
    ThetaRHS_ = rcp(new Epetra_Vector(Map));
    EVyVy_  = rcp(new Epetra_Vector(Map));
    MyPID_    = Comm->MyPID();

    // Set initial condition from PDEAssembler
    *u_  = *pde_.getInitialCondition();
    *u0_ = *u_;

    // Build initial forcing
    forcing_.createW();

    // Report chosen criterion
    if (MyPID_ == 0) {
        const std::string cnames[] = {
            "Absolute", "Relative Fixed", "Relative Adaptive", "Eisenstat-Walker"};
        std::cout << "[Mean Newton] Convergence criterion: "
                  << cnames[static_cast<int>(criterion_)] << "\n";
    }

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

/**
 * @brief Assemble the theta-method residual for the mean Burgers equation.
 *
 * Computes  R = θ·F(u) + (1−θ)·F(u₀) + (u − u₀)  where F is the
 * spatially discretised RHS including diffusion, advection, and stochastic
 * correction.  The residual is stored in ThetaRHS_.
 */
void MeanSolver::ThetaStepper()
{
    pde_.assembleRHS(u_, EVyVy_, *dt_);
    ThetaRHS_->Update(theta_, pde_.getRHSRef(), 0.0);

    pde_.assembleRHS(u0_, EVyVy_, *dt_);
    ThetaRHS_->Update(1 - theta_, pde_.getRHSRef(), 1.0);
    ThetaRHS_->Update(1.0, *u_, -1.0, *u0_, -1.0);
}

// =====================================================================================
// -- Criterion helpers ----------------------------------------------------------------

MeanSolver::NewtonCriterion
MeanSolver::parseCriterion(const std::string& s)
{
    if (s == "Absolute")          return NewtonCriterion::Absolute;
    if (s == "Relative Fixed")    return NewtonCriterion::RelativeFixed;
    if (s == "Relative Adaptive") return NewtonCriterion::RelativeAdaptive;
    if (s == "Eisenstat-Walker")  return NewtonCriterion::EisenstatWalker;
    // Unknown string: warn and fall back to Relative Fixed.
    std::cout << "[Mean Newton] WARNING: unknown criterion '" << s
              << "' — falling back to 'Relative Fixed'.\n";
    return NewtonCriterion::RelativeFixed;
}

bool MeanSolver::checkConvergence(double normFk, double normF0,
                                  double normFprev, double etaEW) const
{
    switch (criterion_) {
    case NewtonCriterion::Absolute:
        return normFk < toleranceRHS_;

    case NewtonCriterion::RelativeFixed:
        return (normFk < toleranceRHS_) || (normFk / normF0 < tolRelative_);

    case NewtonCriterion::RelativeAdaptive: {
        // tol = alpha * dt^(p-1);  p=schemeOrder_, alpha=alphaAdapt_
        const double adaptTol = alphaAdapt_ * std::pow(*dt_, schemeOrder_ - 1);
        return (normFk < toleranceRHS_) || (normFk / normF0 < adaptTol);
    }
    case NewtonCriterion::EisenstatWalker:
        // etaEW is updated by caller; convergence when ||F_k|| <= eta * ||F_0||
        return (normFk < toleranceRHS_) || (normFk <= etaEW * normF0);
    }
    return false; // unreachable
}

// =====================================================================================
int MeanSolver::LinSolve(Epetra_Vector& LHS, Epetra_Vector& RHS)
{
    solver_->factorize();
    return solver_->solve(LHS, RHS, "Mean solve");
}

/**
 * @brief Full Newton iteration for the implicit mean-field Burgers step.
 *
 * Starting from the previous solution u₀, iterates
 * @f$ u^{k+1} = u^k + \Delta u @f$ with backtracking until the theta-method
 * residual drops below toleranceRHS_ or maxNumIterations_ is reached.
 *
 * @return @c true if the Newton loop converged.
 */
bool MeanSolver::NewtonSolver()
{
    isConverged_ = false;
    dx_->PutScalar(0.0);

    // Evaluate initial residual — used as normalisation base for relative criteria.
    ThetaStepper();
    ThetaRHS_->Scale(-1.0);
    ThetaRHS_->Norm2(&NormRHS_);
    const double normF0   = (NormRHS_ > 0.0) ? NormRHS_ : 1.0;
    double normFprev      = normF0;  // for Eisenstat-Walker tracking
    double etaEW          = 1.0;     // Eisenstat-Walker force tolerance (init to 1)

    if (debug_) {
        // Print which effective tolerance will be used this step
        std::string critLabel;
        switch (criterion_) {
            case NewtonCriterion::Absolute:
                critLabel = "abs<" + std::to_string(toleranceRHS_); break;
            case NewtonCriterion::RelativeFixed:
                critLabel = "rel<" + std::to_string(tolRelative_); break;
            case NewtonCriterion::RelativeAdaptive: {
                double atol = alphaAdapt_ * std::pow(*dt_, schemeOrder_ - 1);
                critLabel = "rel-adaptive<" + std::to_string(atol); break;
            }
            case NewtonCriterion::EisenstatWalker:
                critLabel = "EW(γ=" + std::to_string(ewGamma_) + ")"; break;
        }
        std::cout << "\nMean:: Newton: initial ||F|| = " << normF0
                  << "  criterion: " << critLabel;
    }

    for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
        pde_.assembleJacobian(u_, *dt_, theta_);
        LinSolve(*dx_, *ThetaRHS_);
        u_->Update(1.0, *dx_, 1.0);

        ThetaStepper();
        ThetaRHS_->Scale(-1.0);
        ThetaRHS_->Norm2(&NormRHStest_);

        // Update Eisenstat-Walker eta for next iteration
        if (criterion_ == NewtonCriterion::EisenstatWalker && iter_ > 0) {
            const double etaCandidate =
                ewGamma_ * std::pow(NormRHStest_ / normFprev, ewAlphaEW_);
            // Guard: prevent eta from growing (Eisenstat-Walker safeguard)
            etaEW = std::min(etaCandidate, ewGamma_ * etaEW * etaEW);
            etaEW = std::max(etaEW, 1.0e-14); // floor to avoid underflow
        }

        if (debug_)
            std::cout << "\nMean:: Newton: iter " << iter_
                      << "  ||F|| = " << NormRHStest_
                      << "  ||F||/||F0|| = " << NormRHStest_ / normF0
                      << (criterion_ == NewtonCriterion::EisenstatWalker
                              ? "  eta_EW = " + std::to_string(etaEW) : "");

        if (checkConvergence(NormRHStest_, normF0, normFprev, etaEW)) {
            if (debug_) std::cout << "\nMean:: Newton: converged.";
            isConverged_ = true;
            break;
        }

        if (backTracking_ && (NormRHS_ < NormRHStest_)) {
            runBackTracking();
            ThetaRHS_->Scale(-1.0);  // restore sign after backtracking
        }
        normFprev = NormRHStest_;
        NormRHS_  = NormRHStest_;
    }

    // ---- Post-loop: always accept the step; only warn if not converged -----------
    *u0_ = *u_;

    if (!isConverged_) {
        const double relRes = NormRHStest_ / normF0;
        if (NormRHStest_ < tolAcceptable_) {
            if (MyPID_ == 0)
                std::cout << "\n[Mean Newton] WARNING t=" << *t_
                          << ": primary criterion not met."
                          << "  ||F|| = " << NormRHStest_
                          << "  ||F||/||F0|| = " << relRes
                          << "  (within acceptable tol " << tolAcceptable_
                          << ") — continuing.";
            isConverged_ = true;
        } else {
            if (MyPID_ == 0)
                std::cout << "\n[Mean Newton] WARNING t=" << *t_
                          << ": acceptable tol (" << tolAcceptable_
                          << ") not met.  ||F|| = " << NormRHStest_
                          << "  ||F||/||F0|| = " << relRes
                          << " — step accepted with degraded accuracy.";
        }
    }
    return true; // simulation always continues
}

// =====================================================================================
void MeanSolver::runBackTracking()
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
void MeanSolver::WriteSolution(std::string filename, double param,
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
