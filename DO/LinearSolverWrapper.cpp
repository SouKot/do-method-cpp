/*
 * =====================================================================================
 *
 *       Filename:  LinearSolverWrapper.cpp
 *
 *    Description:  Implementation of LinearSolverWrapper — see header for design notes.
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
#include "LinearSolverWrapper.hpp"

#include "Amesos_ConfigDefs.h"
#include "BelosSolverFactory.hpp"
#include "BelosTypes.hpp"
#include "Epetra_Time.h"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_ParameterList.hpp"
#include <iostream>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

// =====================================================================================
LinearSolverWrapper::LinearSolverWrapper(
        const std::string& solverType,
        const std::string& solverName,
        const RCP<Teuchos::ParameterList>& SolParams,
        const RCP<Epetra_CrsMatrix>& matrix,
        bool  test,
        bool  debug,
        int   dbgLvl)
    : solver_type_(solverType)
    , test_(test)
    , debug_(debug)
    , dbgLvl_(dbgLvl)
    , matrix_(matrix)
{
    const Epetra_Comm& comm = matrix_->Comm();

    if (solver_type_ == "Amesos")
    {
        v_prob_ = rcp(new Epetra_LinearProblem);
        v_prob_->SetOperator(matrix_.get());
        Amesos Factory;
        v_solve_ = rcp(Factory.Create(solverName.c_str(), *v_prob_), false);
        v_solve_->SetUseTranspose(false);
        v_solve_->SetParameters(*SolParams);

        Epetra_Time timer(comm);
        v_solve_->SymbolicFactorization();
        if (test_)
        {
            double LocTime = timer.ElapsedTime();
            double time    = 0.0;
            comm.SumAll(&LocTime, &time, 1);
            if (comm.MyPID() == 0)
            {
                std::cout << "\n*******************************************\n"
                          << "Average time per proc. for symbolic factorization"
                          << " = " << time / comm.NumProc() << " sec\n"
                          << "*******************************************\n";
            }
        }
    }
    else if (solver_type_ == "Amesos2")
    {
        amesos2_solver_ = Amesos2::create<MAT, MV>(solverName, matrix_);
        amesos2_solver_->setParameters(SolParams);

        Epetra_Time timer(comm);
        amesos2_solver_->symbolicFactorization();
        if (test_)
        {
            double LocTime = timer.ElapsedTime();
            double time    = 0.0;
            comm.SumAll(&LocTime, &time, 1);
            if (comm.MyPID() == 0)
            {
                std::cout << "\n*******************************************\n"
                          << "Average time per proc. for symbolic factorization"
                          << " = " << time / comm.NumProc() << " sec\n"
                          << "*******************************************\n";
            }
        }
    }
    else // Belos iterative path
    {
        std::string precType   = SolParams->get("Ifpack Preconditioner Name", "any");
        int overlapLevel       = SolParams->get("Overlap Level", 0);
        std::string solverType = SolParams->get("Solver", "any");

        Ifpack Factory;
        prec_ = rcp(Factory.Create(precType, &*matrix_, overlapLevel));
        TEUCHOS_ASSERT(prec_ != Teuchos::null);

        Teuchos::ParameterList ifpackList = SolParams->sublist(precType);
        prec_->SetParameters(ifpackList);
        prec_->Initialize();

        belosPrec_ = rcp(new Belos::EpetraPrecOp(prec_));

        Teuchos::ParameterList belosList = SolParams->sublist(solverType);
        bool leftPrec  = belosList.get("Left Preconditioner", true);
        bool verbose   = belosList.get("Verbose", false);
        int  numrhs    = matrix_->NumMyRows(); // placeholder; updated per solve

        if (verbose)
        {
            belosList.set("Verbosity",
                Belos::Errors + Belos::Warnings +
                Belos::TimingDetails + Belos::StatusTestDetails);
            int freq = belosList.get("output frequency", 0);
            if (freq > 0)
                belosList.set("Output Frequency", freq);
        }
        else
        {
            belosList.set("Verbosity", Belos::Errors + Belos::Warnings);
        }

        problem_ = rcp(new Belos::LinearProblem<ST, MV, OP>());
        problem_->setOperator(matrix_);
        if (leftPrec)
            problem_->setLeftPrec(belosPrec_);
        else
            problem_->setRightPrec(belosPrec_);

        Belos::SolverFactory<ST, MV, OP> belosFactory;
        v_solve_iter_ = belosFactory.create(solverType, rcp(&belosList, false));
    }
}

// =====================================================================================
void LinearSolverWrapper::factorize()
{
    if (solver_type_ == "Amesos")
    {
        AMESOS_CHK_ERRV(v_solve_->NumericFactorization());
    }
    else if (solver_type_ == "Amesos2")
    {
        amesos2_solver_->numericFactorization();
    }
    else // Belos
    {
        IFPACK_CHK_ERRV(prec_->Compute());
    }
}

// =====================================================================================
int LinearSolverWrapper::solve(Epetra_MultiVector& LHS,
                               Epetra_MultiVector& RHS,
                               const std::string&  label)
{
    const Epetra_Comm& comm = matrix_->Comm();
    Epetra_Time timer(comm);

    if (solver_type_ == "Amesos")
    {
        v_prob_->SetLHS(&LHS);
        v_prob_->SetRHS(&RHS);
        int err = v_solve_->Solve();
        if (debug_)
            v_solve_->PrintStatus();
        (void)err;
    }
    else if (solver_type_ == "Amesos2")
    {
        amesos2_solver_->solve(&LHS, &RHS);
    }
    else // Belos
    {
        problem_->setOperator(matrix_);
        problem_->setLHS(rcpFromRef(LHS));
        problem_->setRHS(rcpFromRef(RHS));
        problem_->setProblem();
        v_solve_iter_->setProblem(problem_);
        v_solve_iter_->solve();

        if (v_solve_iter_->isLOADetected())
            std::cout << "\n Loss of accuracy detected!!! \n"
                      << "achieved tol = " << v_solve_iter_->achievedTol();

        if (debug_ && comm.MyPID() == 0)
            std::cout << "\nIterative solution\n"
                      << "Solver performed " << v_solve_iter_->getNumIters()
                      << " iterations.\n"
                      << "Achieved tolerance = " << v_solve_iter_->achievedTol()
                      << std::endl;
    }

    // Optional residual check
    if (debug_ && dbgLvl_ % 2 == 0)
    {
        Epetra_MultiVector rh(LHS);
        matrix_->Multiply(false, LHS, rh);
        rh.Update(-1.0, RHS, 1.0);
        std::vector<double> nrm(rh.NumVectors());
        rh.Norm2(nrm.data());
        std::cout << "\nNorm of residual (" << label << ")";
        for (int i = 0; i < rh.NumVectors(); i++)
            std::cout << "  " << nrm[i];
        std::cout << std::endl;
    }

    if (test_)
    {
        double LocTime = timer.ElapsedTime();
        double time    = 0.0;
        comm.SumAll(&LocTime, &time, 1);
        if (comm.MyPID() == 0)
        {
            std::cout << "\n*******************************************\n"
                      << label << " = " << time / comm.NumProc() << " sec\n"
                      << "*******************************************\n";
        }
    }
    return 0;
}
