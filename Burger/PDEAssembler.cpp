/*
 * =====================================================================================
 *
 *       Filename:  PDEAssembler.cpp
 *
 *    Description:  Finite-difference PDE operators for the 1-D viscous Burgers equation
 *                  @f$ u_t + u\,u_x = \mu\,u_{xx} @f$ with periodic boundary conditions.
 *                  Builds the Laplacian, gradient, and identity stencil matrices and
 *                  assembles the nonlinear RHS and Jacobian at each time step.
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
#include "PDEAssembler.hpp"
#include <EpetraExt_MatrixMatrix.h>
#include "EpetraExt_RowMatrixOut.h"
#include <cmath>
#include <iostream>
#include <vector>

using Teuchos::rcp;
using Teuchos::RCP;

static const double PI = 3.14159265358979323846264338328;

namespace Burger {

// =====================================================================================
PDEAssembler::PDEAssembler(int m, double mu, double theta,
                           const Teuchos::RCP<Epetra_Comm>& Comm)
{
    RCP<Epetra_Map> Map = rcp(new Epetra_Map(m, 0, *Comm));
    LinOp_  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    Op2x_   = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 2));
    Op1x_   = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
    eye_    = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
    x_      = rcp(new Epetra_Vector(*Map));
    ICom_   = rcp(new Epetra_Vector(*Map));

    createLinOp(m, mu, Comm);

    LinOp_->FillComplete();
    Op2x_->FillComplete();
    Op1x_->FillComplete();
    eye_->FillComplete();

    // Allocate Jacobian workspace
    LinJac_   = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    NlinJac_  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    ThetaJac_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    jac_      = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    SpDiag_   = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
    *SpDiag_ = *eye_;
    SpDiag_->FillComplete();

    Op1x_u_ = rcp(new Epetra_Vector(*Map));
    Op2x_u_ = rcp(new Epetra_Vector(*Map));
    rhs_    = rcp(new Epetra_Vector(*Map));
    uu_     = rcp(new Epetra_Vector(*Map));
    tmp2_   = rcp(new Epetra_Vector(*Map));

    TmpMat1_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));
    TmpMat2_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));

    // Build initial LinJac = LinOp
    EpetraExt::MatrixMatrix::Add(*LinOp_, false, 1.0, *LinJac_, 0.0);
    LinJac_->FillComplete();

    // Build initial NlinJac (zero, then fill structure)
    Epetra_Vector rndmvct(*x_);
    rndmvct = *ICom_;
    Op1x_->Multiply(false, rndmvct, *Op1x_u_);
    Op2x_->Multiply(false, rndmvct, *Op2x_u_);
    SpDiag_->ReplaceDiagonalValues(*Op1x_u_);
    EpetraExt::MatrixMatrix::Multiply(*SpDiag_, false, *Op2x_, false, *TmpMat1_, false);
    SpDiag_->ReplaceDiagonalValues(*ICom_);
    EpetraExt::MatrixMatrix::Multiply(*Op2x_, false, *SpDiag_, false, *TmpMat2_);
    EpetraExt::MatrixMatrix::Add(*TmpMat2_, false, 1.0, *NlinJac_, 0.0);
    if (!TmpMat1_->Filled()) TmpMat1_->FillComplete();
    NlinJac_->FillComplete();

    // Build initial jac = LinJac - NlinJac
    EpetraExt::MatrixMatrix::Add(*LinJac_, false, 1.0, *jac_, 0.0);
    EpetraExt::MatrixMatrix::Add(*NlinJac_, false, -1.0, *jac_, 1.0);
    jac_->FillComplete();

    // Build initial ThetaJac = I - dt*theta*jac (with dt=0 initially => ThetaJac = I)
    EpetraExt::MatrixMatrix::Add(*jac_, false, 0.0, *ThetaJac_, 0.0);
    EpetraExt::MatrixMatrix::Add(*eye_, false, 1.0, *ThetaJac_, -1.0);
    ThetaJac_->FillComplete();
}

// =====================================================================================
void PDEAssembler::createLinOp(int m, double mu,
                               const Teuchos::RCP<Epetra_Comm>& Comm)
{
    RCP<Epetra_Map> Map = rcp(new Epetra_Map(m, 0, *Comm));
    double x_end = 1.0, x_in = 0.0;
    double dx = (x_end - x_in) / m;

    Epetra_Vector& x = *x_;
    Epetra_Vector& ICom = *ICom_;
    int GlblRw;
    for (int i = 0; i < x.MyLength(); ++i) {
        GlblRw = x.Map().GID(i);
        x[i] = GlblRw * dx;
    }
    for (int i = 0; i < ICom.MyLength(); ++i) {
        ICom[i] = 0.5 * (std::exp(std::cos(2 * PI * x[i])) - 1.5)
                  * std::sin(2 * PI * (x[i] + 0.37));
    }

    int NumGlobalElements = Map->NumGlobalElements();
    int NumMyElements = Map->NumMyElements();
    int* MyGlobalElements = nullptr;
    Map->MyGlobalElementsPtr(MyGlobalElements);

    std::vector<double> Values(2), GradVal(2);
    std::vector<int> Indices(2);
    int NumEntries;
    double val = 1.0;

    for (int i = 0; i < NumMyElements; ++i) {
        if (MyGlobalElements[i] == 0) {
            Indices[0] = 1;
            Values[0] = 1.0;
            Indices[1] = NumGlobalElements - 1;
            Values[1] = 1.0;
            NumEntries = 2;
            GradVal[0] = 1.0;
            GradVal[1] = -1.0;
        } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
            Indices[0] = NumGlobalElements - 2;
            Values[0] = 1.0;
            Indices[1] = 0;
            Values[1] = 1.0;
            NumEntries = 2;
            GradVal[0] = -1.0;
            GradVal[1] = 1.0;
        } else {
            Indices[0] = MyGlobalElements[i] - 1;
            Values[1] = 1.0;
            Indices[1] = MyGlobalElements[i] + 1;
            Values[0] = 1.0;
            NumEntries = 2;
            GradVal[0] = -1.0;
            GradVal[1] = 1.0;
        }
        LinOp_->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
        Op2x_->InsertGlobalValues(MyGlobalElements[i], NumEntries, &GradVal[0], &Indices[0]);

        Values[0] = -2.0;
        LinOp_->InsertGlobalValues(MyGlobalElements[i], 1, &Values[0], MyGlobalElements + i);
        eye_->InsertGlobalValues(MyGlobalElements[i], 1, &val, MyGlobalElements + i);
    }

    LinOp_->Scale(mu / std::pow(dx, 2));
    Op2x_->Scale(1.0 / (2 * dx));

    // Op1x_ is identity (diagonal)
    *Op1x_ = *eye_;
}

// =====================================================================================
void PDEAssembler::BilinearTerm(Teuchos::RCP<Epetra_Vector> u1,
                                Teuchos::RCP<Epetra_Vector> u2,
                                Teuchos::RCP<Epetra_Vector> u3)
{
    Op1x_->Multiply(false, *u1, *Op1x_u_);
    Op2x_->Multiply(false, *u2, *Op2x_u_);
    u3->Multiply(1.0, *Op1x_u_, *Op2x_u_, 0.0);
}

// =====================================================================================
void PDEAssembler::assembleRHS(const Teuchos::RCP<Epetra_Vector>& u,
                               const Teuchos::RCP<Epetra_Vector>& ExpVyVy,
                               double dt)
{
    LinOp_->Multiply(false, *u, *rhs_);
    uu_->Multiply(1.0, *u, *u, 0.0);
    Op2x_->Multiply(false, *uu_, *tmp2_);
    tmp2_->Scale(0.5);
    rhs_->Update(-1.0, *tmp2_, 1.0);
    rhs_->Update(-1.0, *ExpVyVy, 1.0);
    rhs_->Scale(dt);
}

// =====================================================================================
void PDEAssembler::assembleJacobian(const Teuchos::RCP<Epetra_Vector>& u,
                                    double dt, double theta)
{
    Op1x_->Multiply(false, *u, *Op1x_u_);
    Op2x_->Multiply(false, *u, *Op2x_u_);

    SpDiag_->ReplaceDiagonalValues(*Op1x_u_);
    EpetraExt::MatrixMatrix::Multiply(*SpDiag_, false, *Op2x_, false, *TmpMat1_);
    SpDiag_->ReplaceDiagonalValues(*u);
    EpetraExt::MatrixMatrix::Multiply(*Op2x_, false, *SpDiag_, false, *TmpMat2_);
    EpetraExt::MatrixMatrix::Add(*TmpMat2_, false, 1.0, *NlinJac_, 0.0);

    EpetraExt::MatrixMatrix::Add(*LinOp_, false, 1.0, *jac_, 0.0);
    EpetraExt::MatrixMatrix::Add(*NlinJac_, false, -1.0, *jac_, 1.0);

    // ThetaJac = I - dt*theta*jac
    EpetraExt::MatrixMatrix::Add(*jac_, false, dt * theta, *ThetaJac_, 0.0);
    EpetraExt::MatrixMatrix::Add(*eye_, false, 1.0, *ThetaJac_, -1.0);
}

} // namespace Burger
