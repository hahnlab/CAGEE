#include <random>
#include <ostream>
#include <algorithm>

#include "doctest.h"
#include "easylogging++.h"

#include "discretized_gamma.h"
#include "gamma.h"
#include "sigma.h"

using namespace std;

extern mt19937 randomizer_engine;

discretized_gamma::discretized_gamma(double alpha, int bins) : _alpha(alpha) 
{
    if (bins > 1)
    {
        _sigma_multipliers.resize(bins);
        _gamma_cat_probs.resize(bins);
        get_gamma(_gamma_cat_probs, _sigma_multipliers, alpha); // passing vectors by reference
    }
}

sigma_squared* discretized_gamma::get_random_sigma(const sigma_squared& ss)
{
    discrete_distribution<int> dist(_gamma_cat_probs.begin(), _gamma_cat_probs.end());
    return new sigma_squared(ss, _sigma_multipliers[dist(randomizer_engine)]);       
}

vector<sigma_squared> discretized_gamma::get_discrete_sigmas(const sigma_squared &ss)
{
    vector<sigma_squared> sigmas;
    for (auto m : _sigma_multipliers)
    {
        sigmas.emplace_back(ss, m);
    }

    return sigmas;
}

std::vector<double> discretized_gamma::weight( std::vector<double> likelihoods) const
{
    std::vector<double> result(likelihoods.size());
    std::transform(likelihoods.begin(), likelihoods.end(), _gamma_cat_probs.begin(), result.begin(), std::multiplies<double>());
    return result;
}

void discretized_gamma::write_multipliers(std::ostream& ost, bool single_line) const
{
    for (auto& lm : _sigma_multipliers)
    {
        ost << (single_line ? "" : "  ") << lm << (single_line ? "\t" : ";\n");
    }
}

void discretized_gamma::write_probabilities(std::ostream& ost, bool single_line) const
{
    for (auto& lm : _gamma_cat_probs)
    {
        ost << (single_line ? "" : "  ") << lm << (single_line ? "\t" : ";\n");
    }
}

inline void CHECK_SIGMA_VALUE(double val, const sigma_squared& sigma)
{
    CHECK_EQ(doctest::Approx(val), sigma.get_values()[0]);
}

TEST_CASE("Check sigma multipliers for a given alpha")
{
    double alpha = 0.635735;
    discretized_gamma gamma(alpha, 3);

    sigma_squared unit(1.0);
    auto multipliers = gamma.get_discrete_sigmas(unit);
    CHECK_SIGMA_VALUE(0.0976623, multipliers[0]);
    CHECK_SIGMA_VALUE(0.653525, multipliers[1]);
    CHECK_SIGMA_VALUE(2.24881, multipliers[2]);

    sigma_squared sigma((double)0.613693);
    auto sigmas = gamma.get_discrete_sigmas(sigma);
    REQUIRE_EQ(3, sigmas.size());
    CHECK_SIGMA_VALUE(0.0599347, sigmas[0]);
    CHECK_SIGMA_VALUE(0.401064, sigmas[1]);
    CHECK_SIGMA_VALUE(1.38008, sigmas[2]);
}