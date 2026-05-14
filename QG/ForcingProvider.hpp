/*
 * =====================================================================================
 *
 *       Filename:  ForcingProvider.hpp
 *
 *    Description:  Time-invariant stochastic forcing provider for the QG equation.
 *                  W is computed once from a Gaussian-weighted divergence pattern
 *                  and never updated during time stepping (refreshForcing is a no-op).
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
#ifndef QG_FORCING_PROVIDER_HPP
#define QG_FORCING_PROVIDER_HPP

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MultiVectorOut.h"
#include "Teuchos_RCP.hpp"
#include "../DO/DOUtils.hpp"
#include <cmath>

namespace QG {

/**
 * @brief Time-invariant stochastic forcing for the Quasi-Geostrophic equation.
 *
 * Fills W with a Gaussian-weighted divergence pattern.  Entries corresponding
 * to zero mass-matrix rows (boundary / pressure dofs) are zeroed via the
 * diagonal mass vector provided at construction.
 */
class ForcingProvider {
public:
    /**
     * @brief Construct and compute the forcing.
     * @param numVec   Number of stochastic forcing vectors.
     * @param map      Epetra map for the solution vector.
     * @param diagmass Diagonal of the QG mass matrix; used to mask zero-mass dofs.
     * @param debug    Enable diagnostic output.
     */
    ForcingProvider(int numVec, const Epetra_Map& map,
                    Epetra_Vector& diagmass, bool debug)
        : NumStchFrcVec_(numVec)
    {
        double l = 0.125;
        double C = 1.0;
        int n = map.NumGlobalElements();
        double NX = std::sqrt(n / 2.0);
        W_ = Teuchos::rcp(new Epetra_MultiVector(map, NumStchFrcVec_));
        for (int col = 0; col < W_->NumVectors(); ++col) {
            double etax = C * col;
            double etay = C * (1 - col);
            for (int i = 1; i <= W_->MyLength(); i = i + 2) {
                double j   = i / 2.0;
                double rem = std::fmod(j, NX) / (NX - 1);
                double q   = std::floor(j / NX) / (NX - 1);
                (*(*W_)(col))[i - 1] =
                    std::exp(-2.0 * ((rem - 0.5) * (rem - 0.5) + (q - 0.5) * (q - 0.5))
                             / (4.0 * l * l))
                    * (etay - 2.0 * etay * rem + 2.0 * etax * q - etax) / (2.0 * l * l);
                if (std::abs(diagmass[i - 1]) < 1.0e-14)
                    (*(*W_)(col))[i - 1] = 0.0;
                (*(*W_)(col))[i] = 0.0;
            }
        }
        const char* flnm = "forcing.mm";
        EpetraExt::MultiVectorToMatrixMarketFile(flnm, *W_);
        if (debug)
            printnormMV(*W_, 2, "norm of Original W :");
    }

    /// @brief No-op — QG forcing is time-invariant.
    void refreshForcing(double) {}

    /// @brief Return the number of stochastic forcing vectors.
    int get_dim_W() const { return NumStchFrcVec_; }

    /// @brief Return an RCP to the stochastic forcing multi-vector W.
    Teuchos::RCP<Epetra_MultiVector> get_W() { return W_; }

private:
    int NumStchFrcVec_;
    Teuchos::RCP<Epetra_MultiVector> W_;
};

} // namespace QG

#endif // QG_FORCING_PROVIDER_HPP
