/*
 * =====================================================================================
 *
 *       Filename:  ForcingProvider.hpp
 *
 *    Description:  Time-dependent stochastic forcing provider for Burgers.
 *                  Extracted from Burger::Mean as part of Phase 5 decomposition.
 *                  Computes W(t) = 0.5[cos(4πx)e^{-10t} + cos(2πx)e^{5(t-1)}].
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
#ifndef BURGER_FORCING_PROVIDER_HPP
#define BURGER_FORCING_PROVIDER_HPP

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <cmath>

namespace Burger {

/**
 * @brief Time-dependent stochastic forcing for the Burgers equation.
 *
 * Fills W with the harmonic forcing pattern:
 * @f[
 *   W_i(x,t) = 0.5[\cos(4\pi x)\,e^{-10t} + \cos(2\pi x)\,e^{5(t-1)}]
 * @f]
 *
 * W_ is allocated on first call to createW() and updated in-place on subsequent
 * calls, so shared RCPs obtained via get_W() remain valid.
 */
class ForcingProvider {
public:
    /**
     * @brief Construct the forcing provider.
     * @param numVec  Number of stochastic forcing vectors (columns of W).
     * @param x       Grid coordinate vector (used for spatial pattern).
     */
    ForcingProvider(int numVec, const Teuchos::RCP<Epetra_Vector>& x)
        : NumStchFrcVec_(numVec), x_(x) {}

    /**
     * @brief Compute the time-dependent stochastic forcing at time @p t.
     * @param t  Current simulation time.
     */
    void createW(double t = 0) {
        static const double PI = 3.14159265358979323846264338328;
        if (W_ == Teuchos::null)
            W_ = Teuchos::rcp(new Epetra_MultiVector(x_->Map(), NumStchFrcVec_));
        for (int i = 0; i < W_->NumVectors(); ++i)
            for (int j = 0; j < W_->MyLength(); ++j)
                (*(*W_)(i))[j] = 0.5 * (std::cos(4 * PI * (*x_)[j]) * std::exp(-10 * t)
                                       + std::cos(2 * PI * (*x_)[j]) * std::exp(5 * (t - 1)));
    }

    /// @brief Unified adapter called by the DO time loop.
    void refreshForcing(double t) { createW(t); }

    /// @brief Return the number of stochastic forcing vectors.
    int get_dim_W() const { return NumStchFrcVec_; }

    /// @brief Return an RCP to the stochastic forcing multi-vector W.
    Teuchos::RCP<Epetra_MultiVector> get_W() { return W_; }

private:
    int NumStchFrcVec_;
    Teuchos::RCP<Epetra_Vector> x_;
    Teuchos::RCP<Epetra_MultiVector> W_;
};

} // namespace Burger

#endif // BURGER_FORCING_PROVIDER_HPP
