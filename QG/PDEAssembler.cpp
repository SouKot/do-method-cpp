/*
 * =====================================================================================
 *
 *       Filename:  PDEAssembler.cpp
 *
 *    Description:  PDE operator assembly for the barotropic quasi-geostrophic (QG)
 *                  vorticity equation on a rectangular domain.  Delegates spatial
 *                  discretisation to QG::QG (Fortran/C++ wrapper) and provides the
 *                  theta-method operator M − dt·θ·J for the implicit time step.
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
#include "../DO/DOUtils.hpp"
#include <iostream>

using Teuchos::rcp;
using Teuchos::RCP;

namespace QG {

// =====================================================================================
PDEAssembler::PDEAssembler(int nx, int ny, double reynoldsNum, double topography,
                           double theta, const Teuchos::RCP<Epetra_Comm>& Comm)
    : nx_(nx), ny_(ny), n_(nx * ny * 2)
{
    qg_ = Teuchos::rcp(new ::QG::QG(nx, ny));
    double zeta = 1.0 / n_;
    qg_->set_par(20, zeta);
    qg_->set_par(21, topography);
    qg_->set_par(5, reynoldsNum);
    qg_->set_par(11, 1.0);

    RCP<Epetra_Map> Map = rcp(new Epetra_Map(n_, 0, *Comm));
    mass_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1));

    double diagmassval[n_];
    qg_->mass(diagmassval);
    diagMass_ = rcp(new Epetra_Vector(Epetra_DataAccess::Copy, *Map, diagmassval));
    for (int i = 0; i < n_; i++)
        mass_->InsertGlobalValues(i, 1, &(*diagMass_)[i], &i);
    mass_->FillComplete();
    mass_->Scale(-1.0); // NOTE: mass matrix multiplied by -1.0 !! see QG.m

    // Allocate workspace
    jac_      = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    ThetaJac_ = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 3));
    rhs_      = rcp(new Epetra_Vector(*Map));
    uu_       = rcp(new Epetra_Vector(*Map));

    // Build initial jac and ThetaJac with a random vector to establish sparsity
    Epetra_Vector tmp(*Map);
    tmp.Random();
    assembleJacobian(Teuchos::rcpFromRef(tmp), 0.0, theta);
    ThetaJac_->FillComplete();
}

// =====================================================================================
void PDEAssembler::BilinearTerm(Teuchos::RCP<Epetra_Vector> u1,
                                Teuchos::RCP<Epetra_Vector> u2,
                                Teuchos::RCP<Epetra_Vector> u3)
{
    double *viewU1, *viewU2, *viewU3;
    u1->ExtractView(&viewU1);
    u2->ExtractView(&viewU2);
    u3->ExtractView(&viewU3);
    qg_->bilin(viewU1, viewU2, viewU3);
}

// =====================================================================================
void PDEAssembler::assembleRHS(const Teuchos::RCP<Epetra_Vector>& u,
                               const Teuchos::RCP<Epetra_Vector>& ExpVyVy,
                               double dt)
{
    double* viewRHS;
    double* viewX;
    u->ExtractView(&viewX);
    rhs_->ExtractView(&viewRHS);
    qg_->rhs(viewX, viewRHS);
    rhs_->Update(1.0, *ExpVyVy, -1.0); // f1=ExpVyVy-rhs!! check stochastic.m
    if (debug_) {
        printnormMV(*rhs_, 2, "norm of qg.rhs");
        printnormMV(*ExpVyVy, 2, "PDEAssembler:: norm of ExpVyVy");
    }
    rhs_->Scale(dt);
}

// =====================================================================================
void PDEAssembler::getProblemRHS(Epetra_Vector& u, Epetra_Vector& F)
{
    double* viewRHS;
    double* viewX;
    u.ExtractView(&viewX);
    F.ExtractView(&viewRHS);
    qg_->rhs(viewX, viewRHS);
}

// =====================================================================================
void PDEAssembler::assembleJacobian(const Teuchos::RCP<Epetra_Vector>& u,
                                    double dt, double theta)
{
    int mxcolsz = n_ * 20;
    int numentries;
    double* viewU;
    int* rowmrkr = new int[n_ + 1];
    int* col = new int[mxcolsz];
    double* val = new double[mxcolsz];
    u->ExtractView(&viewU);
    qg_->jacobian(viewU, 0.0, rowmrkr, col, val);

    if (jac_->Filled()) {
        jac_->PutScalar(0.0);
        int* indexOffSet;
        int* indices;
        double* values;
        jac_->ExtractCrsDataPointers(indexOffSet, indices, values);
        for (int i = 0; i < n_ + 1; ++i)
            indexOffSet[i] = rowmrkr[i];
        for (int i = 0; i < indexOffSet[n_]; ++i) {
            indices[i] = col[i];
            values[i] = val[i];
        }
    } else {
        for (int i = 0; i < n_ + 1; ++i) {
            numentries = rowmrkr[i + 1] - rowmrkr[i];
            if (numentries == 0) continue;
            for (int j = rowmrkr[i]; j < rowmrkr[i + 1]; ++j)
                jac_->InsertGlobalValues(i, 1, &val[j], &col[j]);
        }
    }
    if (!jac_->Filled()) jac_->FillComplete();
    jac_->Scale(-1.0); // NOTE: jac is multiplied by -1.0 !!! see QG.m

    delete[] rowmrkr;
    delete[] col;
    delete[] val;

    // ThetaJac = mass - dt*theta*jac
    EpetraExt::MatrixMatrix::Add(*jac_, false, dt * theta, *ThetaJac_, 0.0);
    EpetraExt::MatrixMatrix::Add(*mass_, false, 1.0, *ThetaJac_, -1.0);
}

} // namespace QG
