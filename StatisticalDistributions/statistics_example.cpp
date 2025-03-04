#include "statistics.h"
#include <iostream>
#include <vector>

int main() {
    StandardNormalDistribution snd;
    std::vector<double> uniform_draws(20, 0.0);
    std::vector<double> normal_draws(20, 0.0);

    // Simple random number generation based on RAND 
    for(auto &u : uniform_draws) {
        u = rand() / static_cast<double>(RAND_MAX);
    }

    // Generate random draws from the standard normal distribution
    snd.random_draws(uniform_draws, normal_draws);

    for(const auto &n : normal_draws) {
        std::cout << n << std::endl;
    }
    return 0;
}