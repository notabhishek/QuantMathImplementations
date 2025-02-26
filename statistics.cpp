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
/*
Purpose of inverse CDF: 
    P(Z <= x) = quantile, where Z ~ N(0, 1).
    If you have a quantile, the inverse CDF gives you the corresponding z-score.

Algorithm: Beasly-Springer-Moro algorithm
    There is no closed form solution to the inverse CDF of the standard normal 
    distribution. We will use the Beasley-Springer-Moro algorithm to approximate
    the inverse CDF.

    The algo uses piecewise approximations to accurately approxmiate the inverse as 
    no single formula is accurate across the whole range. 

    1. Central region (0.5 <= quantile <= 0.92)
        - uses a rational (fractional) approximation. we compute 
        - a numerator: sum of odd powered terms of (quantile - 0.5) scaled by coefficients from array a
        - a denominator: sum of even powered terms of (quantile - 0.5) scaled by coefficients from array b
        - inv_cdf(quantile) = numerator / denominator
    2. Uppertail (0.92 < quantile <= 1)
        - extreme tail is tricker to approximate, so a double logarithm is used 
        - x = log(-log(1 - quantile)) 
        - then uses a polynomial expansion using coefficients from array c

    3. Lowertail (0 <= quantile < 0.5)
        - By symmetry, inverse of a quantile below 0.5 is the negative of the inverse of (1-quantile)
        - inv_cdf(quantile) = -inv_cdf(1 - quantile) 
*/
double StandardNormalDistribution::inv_cdf(const double& quantile) const {
    static double a[4] = {   2.50662823884,
                            -18.61500062529,
                            41.39119773534,
                            -25.44106049637};

    static double b[4] = {  -8.47351093090,
                            23.08336743743,
                            -21.06224101826,
                            3.13082909833};

    static double c[9] = {0.3374754822726147,
                            0.9761690190917186,
                            0.1607979714918209,
                            0.0276438810333863,
                            0.0038405729373609,
                            0.0003951896511919,
                            0.0000321767881768,
                            0.0000002888167364,
                            0.0000003960315187};

    // central region approximation
    if (quantile >= 0.5 && quantile <= 0.92) {
        double num = 0.0;
        double denom = 1.0;

        for (int i=0; i<4; i++) {
            num += a[i] * pow((quantile - 0.5), 2*i + 1);
            denom += b[i] * pow((quantile - 0.5), 2*i);
        }

        return num/denom;
    } 
    // upper tail approximation
    else if (quantile > 0.92 && quantile < 1) {
        double num = 0.0;

        for (int i=0; i<9; i++) {
            num += c[i] * pow((log(-log(1-quantile))), i);
        }

        return num;
    } 
    // using symmetry to approximate lower tail
    else {
        return -1.0*inv_cdf(1-quantile);
    }
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
