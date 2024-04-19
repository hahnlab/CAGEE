#include <random>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/uniform.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "prior.h"
#include "matrix_cache.h"

using namespace std;

prior::prior(std::string dist, double p1, double p2) : distribution(dist), param1(p1), param2(p2)
{
    auto supported_distributions = {"gamma", "fisher", "uniform"};
    if (std::find(supported_distributions.begin(), supported_distributions.end(), dist) == supported_distributions.end())
        throw std::domain_error("Prior must be given in the form dist:param1:param2");
}


double prior::pdf(double value) const
{
    if (distribution == "gamma")
    {
        boost::math::gamma_distribution<double> d(param1, param2);
        return boost::math::pdf(d, value);
    }
    else if (distribution == "fisher")
    {
        boost::math::fisher_f_distribution<double> f(param1, param2);
        return boost::math::pdf(f, value);
    }
    else if (distribution == "uniform")
    {
        boost::math::uniform_distribution<double> f(param1, param2);
        return boost::math::pdf(f, value);
    }
    else
    {
        throw std::domain_error("Unknown probability distribution type");
    }
}

vector<double> get_priors(const matrix_cache& calc, boundaries bounds, const prior *p_prior)
{
    vector<double> priors(calc.create_vector().size());
    for (size_t j = 0; j < priors.size(); ++j) {
        double x = (double(j) + 0.5) * double(bounds.second) / (priors.size() - 1);

        priors[j] = computational_space_prior(x, p_prior);
    }

    return priors;
}

double computational_space_prior(double val, const prior *p_prior)
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    return exp(val) * p_prior->pdf(exp(val));
#else
    return p_prior->pdf(val);
#endif

}

double compute_prior_likelihood(const vector<double>& partial_likelihood, const vector<double>& priors)
{
    std::vector<double> full(partial_likelihood.size());
    std::transform(partial_likelihood.begin(), partial_likelihood.end(), priors.begin(), full.begin(), std::multiplies<double>());
    std::transform(full.begin(), full.end(), full.begin(), [](double d) {
        return isnan(d) ? -numeric_limits<double>::infinity() : d;
        });

#ifdef USE_MAX_PROBABILITY
    double likelihood = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
#else
    double likelihood = accumulate(full.begin(), full.end(), 0.0, [](double a, double b) { return isinf(b) ? a : a+b; }); // sum over all sizes (Felsenstein's approach)
#endif
    return likelihood;
}


TEST_CASE("prior returns correct pdf values for gamma")
{
    prior p("gamma", 0.375, 1600);
    CHECK_EQ(doctest::Approx(0.0133233), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.00719489), p.pdf(8));
    CHECK_EQ(doctest::Approx(0.00505239), p.pdf(14));
}

TEST_CASE("prior returns correct pdf values for fisher")
{
    prior p("fisher", 0.75, 0.75);
    CHECK_EQ(doctest::Approx(0.0388044), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.0114423), p.pdf(8));
    CHECK_EQ(doctest::Approx(0.0054983), p.pdf(14));
}

TEST_CASE("prior returns correct pdf values for uniform")
{
    prior p("uniform", 0, 10);
    CHECK_EQ(doctest::Approx(0.1), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.1), p.pdf(8));
    CHECK_EQ(doctest::Approx(0), p.pdf(14));
}

TEST_CASE("prior throws on unknown")
{
    CHECK_THROWS_AS(prior("", 0,0), std::domain_error);
    prior p;
    CHECK_THROWS_AS(p.pdf(3), std::domain_error);
}

TEST_CASE("get_priors")
{
    matrix_cache calc;
    auto p_prior = new prior("gamma", 1.0, 1600);
    auto bounds = boundaries(0, 20);
    auto priors = get_priors(calc, bounds, p_prior);
    REQUIRE_EQ(200, priors.size());
    CHECK_EQ(doctest::Approx(-7.32816), log(priors[0]));
    CHECK_EQ(doctest::Approx(-1.02372), log(priors[75]));
    CHECK_EQ(doctest::Approx(-182.551), log(priors[125]));
    CHECK_EQ(0, priors[199]);
}

TEST_CASE("compute_prior_likelihood combines prior and inference correctly")
{
    vector<double> inf{ 0.1, 0.2, 0.3};

    vector<double> priors({ 1.43078e-15,    2.5363e-23,  5.65526e-35 });
    double actual = log(compute_prior_likelihood(inf, priors));

#ifdef USE_MAX_PROBABILITY
    CHECK_EQ(doctest::Approx(-35.7683), actual);
#else
    CHECK_EQ(doctest::Approx(-36.4831), actual);
#endif
}