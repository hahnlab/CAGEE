#include <random>
#include <ostream>
#include <algorithm>

#include <boost/math/distributions/gamma.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "discretized_gamma.h"
#include "sigma.h"

using namespace std;

extern mt19937 randomizer_engine;

discretized_gamma::discretized_gamma(double alpha, int bins) : _alpha(alpha) 
{
    if (bins > 1 && alpha > 0.0)
    {
        _sigma_multipliers.resize(bins);
        _gamma_cat_probs.resize(bins);
        fill_n(_gamma_cat_probs.begin(), bins, 1.0 / bins);
        boost::math::gamma_distribution<> gamma(alpha, 1/alpha);

        double probability = _gamma_cat_probs[0]/2.0;
        for (int i = 0; i < bins; ++i)
        {
            _sigma_multipliers[i] = boost::math::quantile(gamma, probability);
            probability += _gamma_cat_probs[i];
        }
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
    CHECK_SIGMA_VALUE(0.0818406, multipliers[0]);
    CHECK_SIGMA_VALUE(0.547651, multipliers[1]);
    CHECK_SIGMA_VALUE(1.88449, multipliers[2]);

    sigma_squared sigma((double)0.613693);
    auto sigmas = gamma.get_discrete_sigmas(sigma);
    REQUIRE_EQ(3, sigmas.size());
    CHECK_SIGMA_VALUE( 0.050225, sigmas[0]);
    CHECK_SIGMA_VALUE(0.336089, sigmas[1]);
    CHECK_SIGMA_VALUE(1.1565, sigmas[2]);
 }

 
TEST_CASE("get_random_sigma uses multiplier based on category probability")
{
    discretized_gamma gamma(1.0, 3);

    vector<double> results(100);
    generate(results.begin(), results.end(), [&gamma]() {
        unique_ptr<sigma_squared> new_lam(gamma.get_random_sigma(sigma_squared(1.2)));
        return new_lam->get_values()[0];
        });

    CHECK_EQ(doctest::Approx(0.909111), accumulate(results.begin(), results.end(), 0.0) / 100.0);

}