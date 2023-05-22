#include <random>

#include "doctest.h"
#include "easylogging++.h"

#include "prior.h"

double gammapdf(double value, const std::gamma_distribution<double>& dist) {
    double k = dist.alpha();
    double theta = dist.beta(); // see discussion in root_distribution_gamma

    double a = pow(theta, k);
    double b = pow(value, (k - 1));
    double c = exp(-1 * value / theta);
    double d = tgamma(k);

    return (b * c) / (a * d);
}

double prior::pdf(double value) const
{
    std::gamma_distribution<double> p(param1, param2);
    return gammapdf(value, p);
}

TEST_CASE("prior returns correct pdf values for gamma")
{
    prior p("gamma", 0.375, 1600);
    CHECK_EQ(doctest::Approx(0.0133233), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.00719489), p.pdf(8));
    CHECK_EQ(doctest::Approx(0.00505239), p.pdf(14));
}

