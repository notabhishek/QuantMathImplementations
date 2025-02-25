#ifndef __STATISTICS_CPP
#define __STATISTICS_CPP 

#define _USE_MATH_DEFINES

#include "statistics.h"
#include <iostream>

StatisticalDistribution::StatisticalDistribution() { }
StatisticalDistribution::~StatisticalDistribution() { }

// Constructor / Destructor 
StandardNormalDistribution::StandardNormalDistribution() { }
StandardNormalDistribution::~StandardNormalDistribution() { }

// Probability density function 
// The probability density function (pdf) of the standard normal distribution is: 
// PDF(x) = (1 / sqrt(2 * pi)) * e ^ (-x^2 / 2)

double StandardNormalDistribution::pdf(const double& x) const {
    return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-x * x / 2.0);
}

// Cumulative density function 
double StandardNormalDistribution::cdf(const double& x) const {
    // TODO
}

// Inverse cumulative density function (aka probit function)
double StandardNormalDistribution::inv_cdf(const double& quantile) const {
    // TODO
}

// Expectation/mean 
double StandardNormalDistribution::mean() const { return 0.0; }

// Variance 
double StandardNormalDistribution::var() const { return 1.0; }

// Standard deviation 
double StandardNormalDistribution::stdev() const { return 1.0; }

// Obtain a sequence of random draws from this distribution 
void StandardNormalDistribution::random_draws(const std::vector& uniform_draws,
    std::vector& dist_draws) {
    // TODO
}

#endif
