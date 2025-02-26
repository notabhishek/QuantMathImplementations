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
// The cumulative density function (cdf) of the standard normal distribution is: 
// CDF(x) = integral from -infinity to x of PDF(t) dt
// There is no closed form solution to this, so we will use an approximation 
/*
    For x < 0, we use symmetry: CDF(x) = 1 - CDF(-x)
    For x >= 0, we use the approximation: 
        1. Transformation: k = 1 / (1 + 0.2316419 * x)
        2. Compute k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))))
        3. CDF(x) = 1 - (1 / sqrt(2 * pi)) * exp(-x^2 / 2) * k_sum
*/
double StandardNormalDistribution::cdf(const double& x) const {
    if(x >= 0) {
        // transformation
        double k = 1.0 / (1.0 + 0.2316419 * x); 
        // compute k_sum, we use horner's method to evaluate the polynomial efficiently
        double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

        return 1.0 - (1.0/sqrt(2*M_PI)) * exp(-0.5*x*x) * k_sum;
    }
    else {
        return 1.0 - cdf(-x);
    }
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
