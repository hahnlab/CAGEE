#ifndef DISCRETIZED_GAMMA_H
#define DISCRETIZED_GAMMA_H

#include <vector>
#include <iosfwd>

class sigma_squared;

class discretized_gamma
{
    std::vector<double> _sigma_multipliers;
    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category
    double _alpha;
public:
    discretized_gamma() = default;
    discretized_gamma(double alpha, int bins);

    sigma_squared* get_random_sigma(const sigma_squared& ss);

    std::vector<sigma_squared> get_discrete_sigmas(const sigma_squared &ss);

    double get_alpha() const { return _alpha; }
   
    void write_multipliers(std::ostream& ost, bool single_line = false) const;
    void write_probabilities(std::ostream& ost, bool single_line = false) const;
    std::vector<double> weight( std::vector<double> likelihoods) const;
};

#endif
