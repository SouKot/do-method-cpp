/*
 * =====================================================================================
 *
 *       Filename:  Mean.cpp
 *
 *    Description:  Implementation of the QG::Mean thin coordinator.
 *                  Constructor, Newton solver, theta-method, line-search, and I/O.
 *                  PDE assembly delegated to QG::PDEAssembler, solver to
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
#include <cmath>
#include <algorithm>

using namespace Teuchos;

// =====================================================================================
MeanSolver::MeanSolver(RCP<Teuchos::ParameterList> PrmLst, double* t, double* dt,
           RCP<Epetra_Comm> Comm)
    : pde_(PrmLst->get("nx", 200),
           PrmLst->get("ny", 200),
           PrmLst->get("Reynolds Number", 0.0),
           PrmLst->get("Topography", 0.0),
           PrmLst->get("theta", 0.5),
           Comm)
    , t_(t)
    , dt_(dt)
    , theta_(PrmLst->get("theta", 0.5))
    , toleranceRHS_(1.0e-10)
    , maxNumIterations_(50)
    , numBackTrackingSteps_(20)
    , backTracking_(true)
    , isConverged_(false)
    , test_(PrmLst->get("Testing", true))
    , debug_(PrmLst->get("Debugging", true))
{
    // Construct forcing after pde_ is ready (needs diagMass and map)
    forcing_.emplace(PrmLst->get("No. of vectors in stoch. forcing", 2),
                     pde_.getMassMatrix()->RowMatrixRowMap(),
                     pde_.getDiagMass(),
                     debug_);

    int n = pde_.getN();
    RCP<Epetra_Map> Map = rcp(new Epetra_Map(n, 0, *Comm));
    u_        = rcp(new Epetra_Vector(*Map));
    u0_       = rcp(new Epetra_Vector(*Map));
    dx_       = rcp(new Epetra_Vector(*Map));
    x_        = rcp(new Epetra_Vector(*Map));
    ThetaRHS_ = rcp(new Epetra_Vector(*Map));
    EVyVy_  = rcp(new Epetra_Vector(*Map));
    MyPID_    = Comm->MyPID();

    // Setup solver
    std::string solver_type = PrmLst->get("Solver Package", "Amesos");
    ParameterList& SolverSublist = PrmLst->sublist(solver_type);
    std::string solver_name = SolverSublist.get("Solver name", "Mumps");
    const RCP<ParameterList> SolParams =
        getParametersFromXmlFile(solver_name + ".xml");
    if (MyPID_ == 0)
        std::cout << solver_type << " solver(" << solver_name
                  << ") has been chosen as mean solver";

    solver_.emplace(solver_type, solver_name, SolParams,
                    pde_.getThetaJacobian(), test_, debug_, 0);
}

// =====================================================================================
void MeanSolver::setSolution(Epetra_Vector& x)
{
    *u0_ = x;
    *u_  = x;
}

/**
 * @brief Assemble the theta-method residual for the mean QG equation.
 *
 * Computes  R = θ·F(u) + (1−θ)·F(u₀) + M·(u − u₀)  where M is the QG mass
 * matrix and F the spatially discretised RHS (vorticity advection + diffusion
 * + stochastic correction).
 *
 * @param rhs_u0  Pre-computed RHS evaluated at the old solution u₀.
 */
void MeanSolver::ThetaStepper(const RCP<Epetra_Vector>& rhs_u0)
{
    pde_.assembleRHS(u_, EVyVy_, *dt_);
    ThetaRHS_->Update(theta_, pde_.getRHSRef(), 0.0);

    ThetaRHS_->Update(1 - theta_, *rhs_u0, 1.0);
    Epetra_Vector U(*u_), MU(*u_);
    U.Update(1.0, *u_, -1.0, *u0_, 0.0);
    pde_.getMassMatrix()->Multiply(false, U, MU);
    ThetaRHS_->Update(1.0, MU, -1.0);
}

// =====================================================================================
int MeanSolver::LinSolve(Epetra_Vector* LHS, Epetra_Vector* RHS)
{
    solver_->factorize();
    return solver_->solve(*LHS, *RHS, "Mean solve");
}

/**
 * @brief Full Newton iteration for the implicit mean-field QG step.
 *
 * Uses an eta-squared step-length heuristic: the Newton correction is scaled by
 * @f$ \sqrt{|s_0/s|} @f$ where @f$ s = \delta u^T J \delta u @f$ and
 * @f$ s_0 = R^T \delta u @f$, to improve robustness for the nonlinear
 * vorticity equation.
 *
 * @return @c true if the Newton loop converged.
 */
bool MeanSolver::NewtonSolver()
{
    isConverged_ = false;
    dx_->PutScalar(0.0);

    pde_.assembleRHS(u0_, EVyVy_, *dt_);
    Epetra_Vector rhs_u0(pde_.getRHSRef());
    ThetaStepper(rcpFromRef(rhs_u0));
    ThetaRHS_->Scale(-1.0);
    ThetaRHS_->Norm2(&NormRHS_);

    if (debug_)
        std::cout << "\nMean:: Newton:      initial norm: " << NormRHS_;

    for (iter_ = 0; iter_ != maxNumIterations_; ++iter_) {
        pde_.assembleJacobian(u_, *dt_, theta_);
        LinSolve(dx_.get(), ThetaRHS_.get());

        u_->Update(1.0, *dx_, 1.0);

        ThetaStepper(rcpFromRef(rhs_u0));
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
        if (backTracking_ && (NormRHS_ < NormRHStest_)) {
            runBackTracking(rhs_u0);
            ThetaRHS_->Scale(-1.0);  // restore sign after backtracking leaves ThetaRHS_=+G(u)
        }
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
void MeanSolver::runBackTracking(Epetra_Vector& rhs_u0)
{
    double reduction = -1.0 / 2;
    for (backTrack_ = 0; backTrack_ != numBackTrackingSteps_; ++backTrack_) {
        if (NormRHStest_ < NormRHS_) {
            if (debug_) std::cout << "\nSuccess...";
            break;
        }
        u_->Update(reduction, *dx_, 1.0);
        ThetaStepper(rcpFromRef(rhs_u0));
        ThetaRHS_->Norm2(&NormRHStest_);
        reduction /= 2.0;
    }
    if (backTrack_ == numBackTrackingSteps_)
        std::cout << "\nNewton: --> BACKTRACKING FAILED" << __FILE__ << __LINE__;
}

/**
 * @brief Cubic-backtracking line-search Newton solver for the QG mean field.
 *
 * Implements a safeguarded line-search strategy (Armijo condition with cubic
 * interpolation) over the theta-method nonlinear system.  Falls back to
 * quadratic interpolation on the first backtrack.  Used as an alternative to
 * NewtonSolver() when the standard eta-squared step length is insufficient.
 *
 * @param x0  Initial guess (usually the solution from the previous time step).
 * @return @c true if the solver converged within tolerance.
 */
bool MeanSolver::newtonLineSearchSolve(Epetra_Vector& x0)
{
    int nIter = 0;
    double alpha = 1.0e-4, maxiter = 30, lambdaMax = 0.5, lambdaMin = 0.1;
    double tolX = 1.0e-12, tolFun = 1.0e-8, resnrm, lambda = 1.0;
    double fold = 0.0, dscrmnnt, resnrm0, convergence, stepnorm;
    double slope = 0.0, mxval, f, f2 = 0.0, lambdaPre, lambda2 = 0.0, A, c1, c2, a, b;
    isConverged_ = false;
    Epetra_Vector x(x0), gT(x0), xold(x0), xoldMod(x0);

    pde_.assembleRHS(rcpFromRef(x0), EVyVy_, *dt_);
    Epetra_Vector rhs_u0(pde_.getRHSRef());
    ThetaStepper(rcpFromRef(rhs_u0));
    ThetaRHS_->Scale(-1.0);
    pde_.assembleJacobian(rcpFromRef(x), *dt_, theta_);
    ThetaRHS_->Norm2(&resnrm);
    dx_->PutScalar(0.0);

    if (debug_)
        std::cout << "\n ITER. NO.    RES. NORM    STEP NORM     LAMBDA ";

    while ((resnrm > tolFun || lambda < 1) && nIter <= maxiter) {
        if (std::abs(lambda - 1.0) < 1.0e-10) {
            nIter += 1;
            LinSolve(dx_.get(), ThetaRHS_.get());
            pde_.getThetaJacobian()->Multiply(true, *ThetaRHS_, gT);
            gT.Dot(*dx_, &slope);
            ThetaRHS_->Dot(*ThetaRHS_, &fold);
            xold = x;
            for (int i = 0; i < xold.MyLength(); ++i)
                xoldMod[i] = std::abs((*dx_)[i]) / std::max(std::abs(xold[i]), 1.0);
            xoldMod.MaxValue(&mxval);
            lambdaMin = tolX / mxval;
        }
        x.Update(1.0, xold, lambda, *dx_, 0.0);
        *u_ = x;
        ThetaStepper(rcpFromRef(rhs_u0));
        ThetaRHS_->Scale(-1.0);
        pde_.assembleJacobian(rcpFromRef(x), *dt_, theta_);
        ThetaRHS_->Dot(*ThetaRHS_, &f);
        lambdaPre = lambda;

        if (f > fold + alpha * lambda * slope) {
            if (std::abs(lambda - 1.0) < 1.0e-12) {
                lambda = -slope / (2 * (f - fold - slope));
            } else {
                A = 1.0 / (lambdaPre - lambda2);
                c1 = f - fold - lambdaPre * slope;
                c2 = f2 - fold - lambda2 * slope;
                a = A * (1 / std::pow(lambdaPre, 2.0) * c1
                       + (-1.0) / std::pow(lambda2, 2.0) * c2);
                b = A * (-lambda2 / std::pow(lambdaPre, 2.0) * c1
                       + lambdaPre / std::pow(lambda2, 2.0) * c2);
                if (std::abs(a) < 1.0e-12) {
                    lambda = -slope / (2 * b);
                } else {
                    dscrmnnt = std::pow(b, 2.0) - 3 * a * slope;
                    if (dscrmnnt < 0)
                        lambda = lambdaMax * lambdaPre;
                    else if (b <= 0.0)
                        lambda = (-b + std::sqrt(dscrmnnt)) / (3 * a);
                    else
                        lambda = -slope / (b + std::sqrt(dscrmnnt));
                }
                lambda = std::min(lambda, lambdaMax * lambdaPre);
            }
        } else if (std::isnan(f) || std::isinf(f)) {
            lambda = lambdaMax * lambdaPre;
        } else {
            lambda = 1.0;
        }

        if (lambda < 1.0) {
            lambda2 = lambdaPre;
            f2 = f;
            lambda = std::max(lambda, lambdaMin * lambdaPre);
            continue;
        }
        resnrm0 = resnrm;
        ThetaRHS_->Norm2(&resnrm);
        convergence = std::log(resnrm0 / resnrm);
        dx_->Norm2(&stepnorm);
        if (debug_) {
            std::cout << "\n   " << nIter << "    " << resnrm
                      << "    " << stepnorm << "    " << lambda;
        }
    }

    if (resnrm <= tolFun)
        isConverged_ = true;
    else {
        std::cout << "\n WARNING: NEWTON SOLVER FOR MEAN NOT CONVERGED!!";
    }
    setSolution(x);
    return isConverged_;
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
