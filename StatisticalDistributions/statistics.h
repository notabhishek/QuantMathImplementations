#ifndef __STATISTICS_H
#define __STATISTICS_H

#include <cmath>
#include <vector>

// Abstract base class for statistical distributions
class StatisticalDistribution {
  public: 
    StatisticalDistribution();
    virtual ~StatisticalDistribution();

    // Distribution functions 
    virtual double pdf(const double& x) const = 0;
    virtual double cdf(const double& x) const = 0;

    // Inverse cumulative distribution function (quantile function)
    virtual double inv_cdf(const double& quantile) const = 0;

    // Discriptive stats 
    virtual double mean() const = 0;
    virtual double var() const = 0;
    virtual double stdev() const = 0;

    // Obtain a sequence of random draws from this distribution 
    virtual void random_draws(const std::vector<double>& uniform_draws, 
      std::vector<double>& dist_draws) = 0;
};

class StandardNormalDistribution : public StatisticalDistribution {
  public: 
    StandardNormalDistribution();
    virtual ~StandardNormalDistribution();

    // Distribution functions 
    virtual double pdf(const double& x) const; 
    virtual double cdf(const double& x) const;

    // Inverse cumulative distribution function (quantile function)
    virtual double inv_cdf(const double& quantile) const;

    // Discriptive stats 
    virtual double mean() const;
    virtual double var() const;
    virtual double stdev() const;

    // Obtain a sequence of random draws from this distribution 
    virtual void random_draws(const std::vector<double>& uniform_draws,
      std::vector<double>& dist_draws);
};

#endif 