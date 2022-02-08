
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
#include "proportional_variance.h"

extern std::mt19937 randomizer_engine; // seeding random number engine
using namespace std;
using namespace Eigen;

float root_equilibrium_distribution::select_root_value(int family_number)
{
    return proportional_variance::to_computational_space(get_raw_root_value(family_number));
}

void root_distribution_fixed::resize(size_t new_size)
{

}

float root_distribution_fixed::get_raw_root_value(int family_number)
{
    return _fixed_root_value;
}

void root_distribution_uniform::resize(size_t new_size)
{
    _max_value = new_size;
}

float root_distribution_uniform::get_raw_root_value(int family_number)
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

float root_distribution_poisson::get_raw_root_value(int family_number)
{
    if ((size_t)family_number >= _vectorized_distribution.size())
        return 0;

    return _vectorized_distribution[family_number];
}

root_distribution_specific::root_distribution_specific(std::vector<std::pair<float, int>> distribution)
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

float root_distribution_specific::get_raw_root_value(int family_number)
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

void root_distribution_gamma::resize(size_t new_size)
{
}

float root_distribution_gamma::get_raw_root_value(int family_number)
{
    return _dist(randomizer_engine);
}

TEST_CASE("select_root_value returns the raw value in computational space")
{
    root_distribution_fixed ef(5.3);
    CHECK_EQ(doctest::Approx(proportional_variance::to_computational_space(5.3)), ef.select_root_value(1));
}

TEST_CASE("Initializing with a fixed_root_value and resizing should not crash")
{
    root_distribution_fixed ef(5.3);
    ef.resize(9.2);
    CHECK(true);

}
TEST_CASE("root_equilibrium_distribution__resize")
{
    root_distribution_specific rd({ {2,5}, {4,3}, {8, 3} });
    rd.resize(15);
    CHECK_EQ(rd.get_raw_root_value(14), 8);
    CHECK_EQ(rd.get_raw_root_value(15), 0);
}

TEST_CASE("root_equilibrium_distribution__poisson_select_root_size")
{
    root_distribution_poisson pd(0.75, 9);

    CHECK_EQ(1, pd.get_raw_root_value(1));
    CHECK_EQ(1, pd.get_raw_root_value(3));
    CHECK_EQ(2, pd.get_raw_root_value(5));
    CHECK_EQ(2, pd.get_raw_root_value(7));

    CHECK_EQ(0, pd.select_root_value(100));
}

TEST_CASE("root_equilibrium_distribution fixed_root_value")
{
    root_distribution_fixed pd(4.3);

    CHECK_EQ(4.3f, pd.get_raw_root_value(1));
    CHECK_EQ(4.3f, pd.get_raw_root_value(5));
    CHECK_EQ(4.3f, pd.get_raw_root_value(20));
}

TEST_CASE("root_distribution_specific returns correct root values")
{
    randomizer_engine.seed(10);

    root_distribution_specific rd({ {2,5}, {4,3}, {8,3} });

    rd.resize(5);
    CHECK_EQ(rd.get_raw_root_value(0), 2);
    CHECK_EQ(rd.get_raw_root_value(1), 2);
    CHECK_EQ(rd.get_raw_root_value(2), 2);
    CHECK_EQ(rd.get_raw_root_value(3), 4);
    CHECK_EQ(rd.get_raw_root_value(4), 8);
    CHECK_EQ(rd.get_raw_root_value(5), 0);
}

TEST_CASE("root_distribution_gamma select_root_value")
{
    randomizer_engine.seed(10);
    root_distribution_gamma pd(0.75, 2.5);

    CHECK_EQ(doctest::Approx(2.6534f), pd.get_raw_root_value(1));
    CHECK_EQ(doctest::Approx(0.00264f), pd.get_raw_root_value(3));
    CHECK_EQ(doctest::Approx(0.10359f), pd.get_raw_root_value(5));
    CHECK_EQ(doctest::Approx(5.33302f), pd.get_raw_root_value(7));
    CHECK_EQ(doctest::Approx(0.0601f), pd.get_raw_root_value(100));
}

TEST_CASE("root_distribution_gamma throws on zero alpha or beta")
{
    CHECK_THROWS_WITH(root_distribution_gamma pd(0, 1), "Invalid gamma distribution");
    CHECK_THROWS_WITH(root_distribution_gamma pd(1, 0), "Invalid gamma distribution");
}

TEST_CASE("root_distribution_gamma has shape/scale parameters")
{
    root_distribution_gamma pd(0.75, 30.0);

    double sum = 0.0;
    for (int i = 0; i < 1000; ++i)
        sum += pd.get_raw_root_value(5);

    size_t sz = 1000;
    vector<double> v(sz);

    generate(v.begin(), v.end(), [&pd]() {
        return pd.get_raw_root_value(5);
        });

    auto mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
    auto variance = std::accumulate(v.begin(), v.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));
        });

    CHECK_EQ(doctest::Approx(23.27714), mean);
    CHECK_EQ(doctest::Approx(800.523), variance);
}
