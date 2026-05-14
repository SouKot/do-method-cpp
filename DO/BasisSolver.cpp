/*
 * =====================================================================================
 *
 *       Filename:  BasisSolver.cpp
 *
 *    Description:  Implementation of the DO stochastic basis (V) solver.
 *                  Replaces Problem_Interface.cpp.  The class is now called
 *                  BasisSolver; the old name remains usable via the alias in
 *                  Problem_Interface.hpp.
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
#include "BasisSolver.hpp"
#include "StochIO.hpp"

#include <EpetraExt_MatrixMatrix.h>
#include "Amesos_Lapack.h"
#include "AnasaziBasicOrthoManager.hpp"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "HYMLS_MatrixUtils.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#if need_locaInterface == 1
#include "FVM_model_interface.h"
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

// =====================================================================================
BasisSolver::BasisSolver(const Teuchos::RCP<Epetra_CrsMatrix>& A,
                         const Teuchos::RCP<Teuchos::ParameterList>& SolverParams,
                         const Teuchos::RCP<Epetra_Vector>& uMean,
                         int& m, double* t, double& dt,
                         const Teuchos::RCP<Epetra_MultiVector>& W,
                         const Teuchos::RCP<Epetra_CrsMatrix>& mass,
                         int iter,
                         const Teuchos::RCP<StochasticState>& sharedState)
    : A_(A)
    , SolverParams_(SolverParams)
    , udet_(uMean)
    , m_(m)
    , t_(t)
    , dt_(dt)
    , W_(W)
    , stochIter(iter)
    , mass_(mass)
    , sharedState_(sharedState)
    , solver_type_(SolverParams->get("Solver Package", "any"))
{
    MyPID = A_->Comm().MyPID();
    InitFile_    = SolverParams_->get("StochBasisFile", "v2.mm");
    test_        = SolverParams_->get("Class Testing", false);
    debug_       = SolverParams_->get("Class Debugging", false);
    DbgLvl_      = SolverParams_->get("Debug Level", 0);

    Teuchos::ParameterList& SolverSublist = SolverParams_->sublist(solver_type_);
    const std::string solver_name = SolverSublist.get("Solver name", "any");
    const RCP<Teuchos::ParameterList> SolParams =
        Teuchos::getParametersFromXmlFile(solver_name + ".xml");

    if (MyPID == 0)
        cout << solver_type_ << " solver(" << solver_name
             << ") has been chosen as stoch basis solver";

    Epetra_Map rowmap(A_->RowMap());
    sharedState_->V       = rcp(new Epetra_MultiVector(rowmap, m_));
    map_    = rcp(new Epetra_Map(m_, 0, A_->Comm()));
    locmap_ = rcp(new Epetra_LocalMap(m_, 0, A_->Comm()));
    y_map   = rcp(new Epetra_Map(W_->NumVectors(), 0, A_->Comm()));
    iteration = 1;

    EVyVyBasis       = rcp(new Epetra_MultiVector(*sharedState_->V));
    identity         = rcp(new Epetra_MultiVector(*map_, m_));
    Ezy      = rcp(new Epetra_MultiVector(*y_map, m_));
    Eyy      = rcp(new Epetra_MultiVector(*map_, m_));
    Eyyy     = rcp(new Epetra_MultiVector(*map_, m_));
    EyyInv  = rcp(new Epetra_MultiVector(*map_, m_));
    RHS_block_1 = rcp(new Epetra_MultiVector(A_->RowMap(), m_));
    sharedState_->EDEyy   = rcp(new Epetra_MultiVector(A->RowMap(), m_));
    Rfactor        = rcp(new Epetra_MultiVector(*locmap_, m_));

    {
        double* r_val;
        Rfactor->ExtractView(&r_val, &m_);
        R = rcp(new Teuchos::SerialDenseMatrix<int, double>(
            Teuchos::DataAccess::View, r_val, m_, m_, m_));
    }

    // Build LHS_block_1_ = M - dt*A
    LHS_block_1_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, A_->RowMap(), 3));
    { Epetra_CrsMatrix* lhs = LHS_block_1_.get();
      EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *A_, false, -(dt_), lhs); }
    LHS_block_1_->FillComplete();

    // Emplace LinearSolverWrapper now that LHS_block_1_ is ready.
    solver_.emplace(solver_type_, solver_name, SolParams, LHS_block_1_,
                    test_, debug_, DbgLvl_);

    double one = 1.0;
    for (int i = 0; i < m_; i++)
        identity->ReplaceGlobalValue(i, i, one);
}

// =====================================================================================
BasisSolver::~BasisSolver() {}

// =====================================================================================
bool BasisSolver::computeJacobian(const Epetra_Vector& /*x*/, Epetra_Operator& Jac)
{
    Jac = Teuchos::dyn_cast<Epetra_Operator>(getJacobian());
    return true;
}

// =====================================================================================
bool BasisSolver::computeShiftedMatrix(double /*alpha*/, double /*beta*/,
                                       const Epetra_Vector& /*x*/, Epetra_Operator& /*A*/)
{
    return false;
}

// =====================================================================================
int BasisSolver::SolveV(Epetra_MultiVector& LHS, Epetra_MultiVector& RHS,
                        const std::string& label)
{
    return solver_->solve(LHS, RHS, label);
}

// =====================================================================================
void BasisSolver::syncJacobian(const Epetra_CrsMatrix& newA)
{
    // Replace A_ values in place; sparsity pattern must match.
    *A_ = newA;
}

// =====================================================================================
void BasisSolver::init_v(double /*t*/)
{
    if (InitFile_ == "None")
    {
        sharedState_->V->Random();
#if need_locaInterface == 1
        Epetra_Vector massDiag(sharedState_->V->Map());
        mass_->ExtractDiagonalCopy(massDiag);
        for (int i = 0; i < massDiag.MyLength(); i++)
        {
            if (abs(massDiag[i]) < 10e-10)
                for (int j = 0; j < sharedState_->V->NumVectors(); j++)
                    (*sharedState_->V)[j][i] = 0.0;
        }
#endif
        MOrth();
        if (debug_ && DbgLvl_ % 3 == 0)
            printnormMV(*sharedState_->V, 2, "norm of V after m-orthogonalization:");
    }
    else
    {
        if (A_->Comm().MyPID() == 0)
            std::cout << "\nInitializing basis from " << InitFile_ << "\n";

        auto vt = StochIO::readMV(InitFile_, sharedState_->V->Map());
        sharedState_->V = vt;
        if (debug_ && A_->Comm().MyPID() == 0)
            printnormMV(*sharedState_->V, 2, "initial norm of V without m-orthogonalization:");
        MOrth();
    }
}

// =====================================================================================
void BasisSolver::MOrth()
{
    Epetra_MultiVector Mmv(*sharedState_->V);
    double time = 0.0, LocTime = 0.0;
    Epetra_Time timer2(A_->Comm());
    R->putScalar(0.0);

    typedef Epetra_MultiVector mv;
    typedef Epetra_Operator    OP;
    mass_->Multiply(false, *sharedState_->V, Mmv);
    A_->Comm().Barrier();

    Belos::IMGSOrthoManager<double, mv, OP> orthman;
    orthman.setOp(mass_);
    int rank = orthman.normalize(*sharedState_->V, rcpFromRef(Mmv), R);
    (void)rank;

    if (test_)
    {
        LocTime = timer2.ElapsedTime();
        A_->Comm().SumAll(&LocTime, &time, 1);
        if (A_->Comm().MyPID() == 0)
            cout << "\n*******************************************\n"
                 << " Avergage time per proc. for Orthogonalization = "
                 << time / A_->Comm().NumProc() << " sec\n"
                 << "*******************************************\n";
    }
    if (debug_ && DbgLvl_ % 5 == 0)
    {
        Teuchos::RCP<mv> temv = rcp(new mv(*sharedState_->V));
        std::vector<double> nrm1(sharedState_->V->NumVectors());
        mass_->Multiply(false, *sharedState_->V, *temv);
        temv->Dot(*sharedState_->V, nrm1.data());
        cout << "\n Norm of V'MV\n";
        for (int i = 0; i < sharedState_->V->NumVectors(); i++)
            std::cout << "  " << nrm1[i];
    }
}

// =====================================================================================
int BasisSolver::v_stoch_init(Epetra_MultiVector* Vold)
{
    Epetra_MultiVector X1(*RHS_block_1);
    Epetra_MultiVector X2(*sharedState_->V);
    std::string label;

    // Numeric factorisation (Amesos/Amesos2) or preconditioner compute (Belos)
    solver_->factorize();

    /**** Solve X1 = J^-1 * RHS_block_1 ****/
    label = debug_ ? "J*X1-RHS_block_1"
                   : "average time per proc. for solving X1=J^-1*RHS_block_1";
    SolveV(X1, *RHS_block_1, label);

    /**** Solve X2 = J^-1 * (M*sharedState_->V) ****/
    Epetra_MultiVector MV(*sharedState_->V);
    label = debug_ ? "J*X2-MV"
                   : "average time per proc. for solving X2=J^-1*MV";
    mass_->Multiply(false, *sharedState_->V, MV);
    SolveV(X2, MV, label);

    /**** Solve Y = (MV'*X2)^-1 * (MV'*X1) ****/
    int numvecV = sharedState_->V->NumVectors();
    double one = 1.0;
    Epetra_LocalMap MVtX_map(numvecV, 0, sharedState_->V->Comm());
    Epetra_CrsMatrix MVtX2(Epetra_DataAccess::Copy, MVtX_map, numvecV);
    for (int i = 0; i < numvecV; i++)
        for (int j = 0; j < numvecV; j++)
            MVtX2.InsertGlobalValues(i, 1, &one, &j);
    MVtX2.FillComplete();

    Epetra_MultiVector MVtX2_temp(MVtX_map, numvecV);
    Epetra_MultiVector MVtX1(MVtX2_temp), Y_sol(MVtX2_temp);
    MVtX2_temp.Multiply('T', 'N', 1.0, MV, X2, 0.0);
    MVtX1.Multiply('T', 'N', 1.0, MV, X1, 0.0);

    for (int i = 0; i < numvecV; i++)
        for (int j = 0; j < numvecV; j++)
        {
            double jacval = MVtX2_temp.operator()(j)->operator[](i);
            MVtX2.ReplaceGlobalValues(i, 1, &jacval, &j);
        }

    Epetra_LinearProblem v_prob2;
    Amesos_Lapack v_solve2(v_prob2);
    v_prob2.SetOperator(&MVtX2);
    v_prob2.SetLHS(&Y_sol);
    v_prob2.SetRHS(&MVtX1);
    v_solve2.SetUseTranspose(false);
    v_solve2.SymbolicFactorization();
    v_solve2.NumericFactorization();
    v_solve2.Solve();

    if (debug_ && DbgLvl_ % 2 == 0)
    {
        Epetra_MultiVector rh(MVtX1);
        HYMLS::MatrixUtils::Dump(MVtX1, "MVtX1");
        HYMLS::MatrixUtils::Dump(MVtX2, "MVtX2");
        HYMLS::MatrixUtils::Dump(Y_sol, "Y_sol");
        MVtX2.Multiply(false, Y_sol, rh);
        rh.Update(-1.0, MVtX1, 1.0);
        std::vector<double> nrm1(rh.NumVectors());
        rh.Norm2(nrm1.data());
        std::cout << "\nNorm of residual (MVtX2*Y_sol-MVtX1)";
        for (int i = 0; i < rh.NumVectors(); i++)
            cout << "  " << nrm1[i];
    }

    /**** X = X1 - X2*Y, then V = Vold + X ****/
    Epetra_MultiVector X_sol(X1);
    X_sol.Multiply('N', 'N', -1.0, X2, Y_sol, 1.0);
    sharedState_->V->Update(1.0, *Vold, 0.0);
    sharedState_->V->Update(1.0, X_sol, 1.0);
    MOrth();

    if (debug_ && DbgLvl_ % 3 == 0)
        printnormMV(*sharedState_->V, 2, "norm of V after m-orthogonalization:");

    return EXIT_SUCCESS;
}

// =====================================================================================
void BasisSolver::computeBlocks(double dt)
{
    dt_ = dt;
    // Rebuild LHS_block_1_ = M - dt*J in place
    { Epetra_CrsMatrix* lhs = LHS_block_1_.get();
      EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *A_, false, -(dt_), lhs); }

    if (debug_ && DbgLvl_ % 7 == 0)
        std::cout << "\nNorm of deterministic Jacobian A = "
                  << A_->NormFrobenius() << std::endl;

    Epetra_Map zeromap_(m_, 0, A_->Comm());
    RCP<Epetra_MultiVector> zeros = rcp(new Epetra_MultiVector(zeromap_, m_));
    zeros->PutScalar(0.0);

    // AV = A * V
    RCP<Epetra_CrsMatrix> Acopy = rcp(new Epetra_CrsMatrix(*A_));
    *Acopy = *A_;
    Epetra_MultiVector AV(*sharedState_->V);
    Acopy->Multiply(false, *sharedState_->V, AV);

    // Fv = dt*AV + sharedState_->EDEyy
    Epetra_MultiVector Fv(AV);
    Fv.Update(dt_, AV, 1.0, *sharedState_->EDEyy, 0.0);

    if (debug_)
    {
        HYMLS::MatrixUtils::Dump(AV, "AV");
        HYMLS::MatrixUtils::Dump(Fv, "FV");
        HYMLS::MatrixUtils::Dump(*sharedState_->EDEyy, "Expect");
        double* FVarray;
        int Ldim;
        Fv.ExtractView(&FVarray, &Ldim);
        double nrm = 0;
        for (int i = 1; i < Fv.MyLength(); i += 2)
            nrm += abs(FVarray[i]) / dt_;
        std::ofstream NormConstrFile("Constraint.txt", std::ios::out | std::ios::app);
        NormConstrFile << std::scientific;
        NormConstrFile.precision(16);
        NormConstrFile << nrm << "\n";
    }

    if (debug_ && DbgLvl_ % 11 == 0)
    {
        std::vector<double> nrm1(m_);
        AV.Norm2(nrm1.data());
        std::cout << "\n nrm of AV";
        for (int i = 0; i < sharedState_->V->NumVectors(); i++) cout << "  " << nrm1[i];
        sharedState_->EDEyy->Norm2(nrm1.data());
        std::cout << "\n nrm of (dt*E[<Vy,Vy>y']+numsubtimestep*E[zy'])/E[yy']";
        for (int i = 0; i < sharedState_->V->NumVectors(); i++) cout << "  " << nrm1[i];
    }

    RHS_block_1->Update(1.0, Fv, 0.0);

    if (debug_ && DbgLvl_ % 11 == 0)
    {
        std::vector<double> nrm1(m_);
        RHS_block_1->Norm2(nrm1.data());
        std::cout << "\n  nrm of RHS_block_1";
        for (int i = 0; i < sharedState_->V->NumVectors(); i++) cout << "  " << nrm1[i];
    }
}

// =====================================================================================
void BasisSolver::TransferNorm()
{
    auto& y = sharedState_->y;
    Epetra_MultiVector ycpy(*y);
    ycpy = *y;
    y->Multiply('N', 'N', 1.0, *Rfactor, ycpy, 0.0);
}

// =====================================================================================
// Legacy stubs (declared in header for link compatibility; not currently used)
void BasisSolver::computeBlocksold(double /*dt*/) {}
void BasisSolver::computeExpectations(double /*dt*/)    {}
void BasisSolver::PostProcess()                   {}
