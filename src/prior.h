#ifndef PRIOR_H
#define PRIOR_H

#include <string>
#include <vector>

class prior
{
    std::string distribution;
    double param1 = 0.0, param2 = 0.0;
public:
    prior(std::string dist, double p1, double p2);
    prior() {}
    double pdf(double value) const;
};

class matrix_cache;
using boundaries = std::pair<double, double>;

double computational_space_prior(double val, const prior *p_prior);
double compute_prior_likelihood(const std::vector<double>& partial_likelihood, const std::vector<double>& priors);
std::vector<double> get_priors(const matrix_cache& calc, boundaries bounds, const prior *p_prior);

#endif
