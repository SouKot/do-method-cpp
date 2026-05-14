/*
 * =====================================================================================
 *
 *       Filename:  NoiseGenerator.hpp
 *
 *    Description:  Compile-time-selected Gaussian noise generator for DO stochastic
 *                  solvers.  Hides the TRNG vs. Boost.Random #if guard behind a
 *                  single sample() call so that Y_Stoch never needs to know which
 *                  library is active.
 *
 *                  Build with -Duse_trng=1 to activate the TRNG back-end; omit or
 *                  set to 0 for Boost.Random (default).
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
#ifndef DO_NOISE_GENERATOR_HPP
#define DO_NOISE_GENERATOR_HPP

#if use_trng==1
#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>
#else
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <random>
#endif

#include "Teuchos_RCP.hpp"

// =====================================================================================
/// @name Gaussian noise generator
// =====================================================================================
/// @{

/**
 * @brief Rank-aware Gaussian noise generator for the DO stochastic solver.
 *
 * Wraps either the TRNG back-end (@c use_trng==1) or Boost.Random under a
 * uniform @c sample() interface.
 *
 * - **TRNG**: calls @c eng_.split(numProc, myPID) for reproducible, independent
 *   parallel streams — one sub-sequence per MPI rank.
 * - **Boost**: seeds a Mersenne-Twister with @c myPID*31 + std::random_device()
 *   for per-rank entropy.
 *
 * Typical usage:
 * @code
 *   NoiseGenerator rng(comm.NumProc(), comm.MyPID());
 *   double v = rng.sample();  // draw from N(0,1)
 * @endcode
 */
class NoiseGenerator {
public:
    /**
     * @brief Seed the engine for the given MPI topology.
     *
     * @param numProc Total number of MPI ranks (TRNG split parameter).
     * @param myPID   Rank of the calling process.
     */
    NoiseGenerator(int numProc, int myPID);

    /// Draw one sample from N(0, 1).
    double sample();

private:
#if use_trng==1
    trng::yarn2 eng_;
    typedef trng::normal_dist<> GEN;
    Teuchos::RCP<GEN> gen_;
#else
    typedef boost::mt19937                          ENG;
    typedef boost::normal_distribution<double>      DIST;
    typedef boost::variate_generator<ENG&, DIST>    GEN;
    ENG                 eng_;
    DIST                dist_;
    Teuchos::RCP<GEN>   gen_;
#endif
};

/// @}

// =====================================================================================
// Inline implementation — header-only because the code path is selected at
// compile time by the use_trng preprocessor macro.
// =====================================================================================

inline NoiseGenerator::NoiseGenerator(int numProc, int myPID)
{
#if use_trng==1
    eng_.split(numProc, myPID);
    gen_ = Teuchos::rcp(new GEN(0.0, 1.0));
#else
    std::random_device rd;
    eng_  = ENG(myPID * 31 + rd());
    dist_ = DIST(0.0, 1.0);
    gen_  = Teuchos::rcp(new GEN(eng_, dist_));
#endif
}

inline double NoiseGenerator::sample()
{
#if use_trng==1
    return (*gen_)(eng_);
#else
    return (*gen_)();
#endif
}

#endif // DO_NOISE_GENERATOR_HPP
