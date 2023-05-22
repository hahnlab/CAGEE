#include <random>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "prior.h"

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
    else
    {
        throw std::domain_error("Unknown probability distribution type");
    }
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

TEST_CASE("prior throws on unknown")
{
    prior p("", 0.375, 1600);
    CHECK_THROWS_AS(p.pdf(4), std::domain_error);
}
