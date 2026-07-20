#ifndef DO_NOISE_GENERATOR_HPP
#define DO_NOISE_GENERATOR_HPP

#if use_trng==1
#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>
#include "Teuchos_RCP.hpp" // Only needed for TRNG pointer now
#else
#include <random> // Replaces all 3 boost headers
#endif


class NoiseGenerator {
public:
    NoiseGenerator(int numProc, int myPID);
    double sample();

private:
#if use_trng==1
    trng::yarn2 eng_;
    typedef trng::normal_dist<> GEN;
    Teuchos::RCP<GEN> gen_;
#else
    std::mt19937 eng_;
    std::normal_distribution<double> dist_;
#endif
};

// =====================================================================================
// Inline implementation 
// =====================================================================================

inline NoiseGenerator::NoiseGenerator(int numProc, int myPID)
{
#if use_trng==1
    eng_.split(numProc, myPID);
    gen_ = Teuchos::rcp(new GEN(0.0, 1.0));
#else
    // Use standard random_device and mt19937
    std::random_device rd;
    
    // Seed the engine. 
    // myPID ensures different ranks get different sequences even if rd() falls back to a clock.
    eng_.seed(myPID * 31 + rd()); 
    
    // Initialize N(0,1) distribution
    dist_ = std::normal_distribution<double>(0.0, 1.0);
#endif
}

inline double NoiseGenerator::sample()
{
#if use_trng==1
    return (*gen_)(eng_);
#else
    // Standard library distributions just take the engine as an argument
    return dist_(eng_); 
#endif
}

#endif // DO_NOISE_GENERATOR_HPP
