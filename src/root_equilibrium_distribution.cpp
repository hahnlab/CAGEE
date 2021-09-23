
#include <numeric>
#include <iostream>
#include <random>
#include <algorithm>

#include "doctest.h"
#include "easylogging++.h"

#include "root_equilibrium_distribution.h"
#include "rootdist_estimator.h"
#include "io.h"
#include "optimizer.h"
#include "gene_transcript.h"
#include "DiffMat.h"
#include "probability.h"

extern std::mt19937 randomizer_engine; // seeding random number engine
using namespace std;
using namespace Eigen;

rootdist_options::rootdist_options(std::string cfg)
{
    auto tokens = tokenize_str(cfg, ':');
    if (tokens[0] == "fixed")
    {
        type = rootdist_type::fixed;
        fixed_value = stof(tokens[1]);
    }
    if (tokens[0] == "gamma")
    {
        type = rootdist_type::gamma;
        dist = gamma_distribution<double>(stof(tokens[1]), 1.0/stof(tokens[2]));
    }
    if (tokens[0] == "file")
    {
        type = rootdist_type::file;
        filename = tokens[1];
    }
}

void root_distribution_fixed::resize(size_t new_size)
{

}

float root_distribution_fixed::select_root_value(int family_number)
{
    return _fixed_root_value;
}

void root_distribution_uniform::resize(size_t new_size)
{
    _max_value = new_size;
}

float root_distribution_uniform::select_root_value(int family_number)
{
    if ((size_t)family_number >= _max_value)
        return 0;

    std::uniform_real_distribution<float> distribution(0.0, _max_value);
    return distribution(randomizer_engine);
}

root_distribution_poisson::root_distribution_poisson(const std::vector<gene_transcript>& gene_families, size_t num_values)
{
    poisson_scorer scorer(gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);
    auto poisson_lambda = result.values[0];

    LOG(INFO) << "Empirical Prior Estimation Result : (" << result.num_iterations << " iterations)";
    LOG(INFO) << "Poisson lambda: " << poisson_lambda << " &  Score: " << result.score;

    for (int i = 0; _vectorized_distribution.size() < num_values; ++i)
    {
        double pct = poisspdf(i, poisson_lambda);
        for (size_t j = 0; j < pct * num_values; ++j)
            _vectorized_distribution.push_back(i + 1);
        _frequency_percentage.push_back(pct);
    }
    // set a few extra percentages beyond the maximum size in the distribution
    for (int i = 0; i < 5; ++i)
        _frequency_percentage.push_back(poisspdf(_frequency_percentage.size(), poisson_lambda));
}

root_distribution_poisson::root_distribution_poisson(double poisson_lambda, size_t num_values)
{
    for (int i = 0; _vectorized_distribution.size() < num_values; ++i)
    {
        double pct = poisspdf(i, poisson_lambda);
        for (size_t j = 0; j < pct * num_values; ++j)
            _vectorized_distribution.push_back(i + 1);
        _frequency_percentage.push_back(pct);
    }
    // set a few extra percentages beyond the maximum size in the distribution
    for (int i = 0; i < 5; ++i)
        _frequency_percentage.push_back(poisspdf(_frequency_percentage.size(), poisson_lambda));
}

void root_distribution_poisson::resize(size_t new_size)
{
}

float root_distribution_poisson::select_root_value(int family_number)
{
    if ((size_t)family_number >= _vectorized_distribution.size())
        return 0;

    return _vectorized_distribution[family_number];
}

root_distribution_specific::root_distribution_specific(std::map<int, float> distribution)
{
    if (distribution.empty())
        throw std::runtime_error("No root distribution specified");

    for (auto it = distribution.begin(); it != distribution.end(); ++it) {
        for (int i = 0; i < it->second; ++i) {
            _vectorized_distribution.push_back(it->first);
        }
    }
    
    auto max = (size_t)*max_element(_vectorized_distribution.begin(), _vectorized_distribution.end()) + 1;
    _frequency_percentage.resize(max);
    for (size_t i = 0; i < max; ++i)
    {
        size_t c = count(_vectorized_distribution.begin(), _vectorized_distribution.end(), i);
        _frequency_percentage[i] = float(c) / float(_vectorized_distribution.size());
    }
}

void root_distribution_specific::resize(size_t new_size)
{
    auto& v = _vectorized_distribution;
    if (new_size < v.size())
    {
        // pare back the distribution randomly
        shuffle(v.begin(), v.end(), randomizer_engine);
        v.erase(v.begin() + new_size, v.end());
    }
    else
    {
        std::uniform_int_distribution<> dis(0, v.size() - 1);
        for (size_t i = v.size(); i < new_size; ++i)
        {
            v.push_back(v.at(dis(randomizer_engine)));
        }
    }
    sort(v.begin(), v.end());
}

float root_distribution_specific::select_root_value(int family_number)
{
    if ((size_t)family_number >= _vectorized_distribution.size())
        return 0;

    return _vectorized_distribution[family_number];
}

root_distribution_gamma::root_distribution_gamma(double alpha, double beta) : _dist(alpha, beta)
{
    if (alpha == 0 || beta == 0)
        throw std::runtime_error("Invalid gamma distribution");
}

float root_distribution_gamma::compute(const gene_transcript& t, double n) const
{
    return gammapdf(n, _dist.alpha(), _dist.beta());
}

void root_distribution_gamma::resize(size_t new_size)
{
}

float root_distribution_gamma::select_root_value(int family_number)
{
    return _dist(randomizer_engine);
}


TEST_CASE("Initializing with a fixed_root_value and resizing should not crash")
{
    root_distribution_fixed ef(5.3);
    ef.resize(9.2);
    CHECK(true);

}
TEST_CASE("root_equilibrium_distribution__resize")
{
    std::map<int, float> m;
    m[2] = 5;
    m[4] = 3;
    m[8] = 3;
    root_distribution_specific rd(m);
    rd.resize(15);
    CHECK_EQ(rd.select_root_value(14), 8);
    CHECK_EQ(rd.select_root_value(15), 0);
}

TEST_CASE("root_equilibrium_distribution__poisson_select_root_size")
{
    root_distribution_poisson pd(0.75, 9);

    CHECK_EQ(1, pd.select_root_value(1));
    CHECK_EQ(1, pd.select_root_value(3));
    CHECK_EQ(2, pd.select_root_value(5));
    CHECK_EQ(2, pd.select_root_value(7));
    CHECK_EQ(0, pd.select_root_value(100));
}

TEST_CASE("root_equilibrium_distribution fixed_root_value")
{
    root_distribution_fixed pd(4.3);

    gene_transcript t;
    t.set_expression_value("A", 5.8);
    t.set_expression_value("B", 19.6);
 
    CHECK_EQ(4.3f, pd.select_root_value(1));
    CHECK_EQ(4.3f, pd.select_root_value(5));
    CHECK_EQ(4.3f, pd.select_root_value(20));
}

TEST_CASE("root_distribution_specific returns correct root values")
{
    randomizer_engine.seed(10);

    std::map<int, float> m;
    m[2] = 5;
    m[4] = 3;
    m[8] = 3;
    root_distribution_specific rd(m);
    //    root_equilibrium_distribution rd(m);
    rd.resize(5);
    CHECK_EQ(rd.select_root_value(0), 2);
    CHECK_EQ(rd.select_root_value(1), 2);
    CHECK_EQ(rd.select_root_value(2), 2);
    CHECK_EQ(rd.select_root_value(3), 4);
    CHECK_EQ(rd.select_root_value(4), 8);
    CHECK_EQ(rd.select_root_value(5), 0);
}

TEST_CASE("root_distribution_gamma select_root_value")
{
    root_distribution_gamma pd(0.75, 2.5);

    CHECK_EQ(doctest::Approx(0.40546f), pd.select_root_value(1));
    CHECK_EQ(doctest::Approx(1.04763f), pd.select_root_value(3));
    CHECK_EQ(doctest::Approx(0.00273f), pd.select_root_value(5));
    CHECK_EQ(doctest::Approx(1.57372f), pd.select_root_value(7));
    CHECK_EQ(doctest::Approx(8.9958f), pd.select_root_value(100));
}

TEST_CASE("root_distribution_gamma compute")
{
    root_distribution_gamma pd(0.75, 2.5);
    gene_transcript t;

    CHECK_EQ(doctest::Approx(0.13318f).scale(1000), pd.compute(t, 1));
    CHECK_EQ(doctest::Approx(0.00068f).scale(1000), pd.compute(t, 2));
    CHECK_EQ(doctest::Approx(0.005).scale(1000), pd.compute(t, 4));
    CHECK_EQ(0.0, pd.compute(t, 100));
}

TEST_CASE("root_distribution_gamma throws on zero alpha or beta")
{
    CHECK_THROWS_WITH(root_distribution_gamma pd(0, 1), "Invalid gamma distribution");
    CHECK_THROWS_WITH(root_distribution_gamma pd(1, 0), "Invalid gamma distribution");
}
