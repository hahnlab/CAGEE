#include <numeric>
#include <iostream>
#include <random>
#include <algorithm>

#include "doctest.h"
#include "easylogging++.h"

#include "root_equilibrium_distribution.h"
#include "poisson.h"
#include "io.h"
#include "optimizer.h"
#include "gene_transcript.h"
#include "DiffMat.h"
#include "probability.h"

using namespace Eigen;
using namespace std;

extern std::mt19937 randomizer_engine; // seeding random number engine

root_equilibrium_distribution::root_equilibrium_distribution(const map<int, int>& root_distribution)
{
    if (root_distribution.empty())
        throw std::runtime_error("No root distribution specified");

    for (auto it = root_distribution.begin(); it != root_distribution.end(); ++it) {
        for (int i = 0; i < it->second; ++i) {
            _vectorized_distribution.push_back(it->first);
        }
    }

    build_percentages();

}

root_equilibrium_distribution::root_equilibrium_distribution(size_t max_size)
{
    _vectorized_distribution.resize(max_size);
    iota(_vectorized_distribution.begin(), _vectorized_distribution.end(), 0);
    build_percentages();
}

root_equilibrium_distribution::root_equilibrium_distribution(double fixed_root_value) : _fixed_root_value(fixed_root_value)
{

}

root_equilibrium_distribution::root_equilibrium_distribution(double poisson_lambda, size_t num_values)
{
    create_from_poisson(poisson_lambda, num_values);
}

root_equilibrium_distribution::root_equilibrium_distribution(const std::vector<gene_transcript>& gene_families, size_t num_values)
{
    poisson_scorer scorer(gene_families);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    LOG(INFO) << "Empirical Prior Estimation Result : (" << result.num_iterations << " iterations)";
    LOG(INFO) << "Poisson lambda: " << result.values[0] << " &  Score: " << result.score;

    create_from_poisson(result.values[0], num_values);

}

void root_equilibrium_distribution::create_from_poisson(double poisson_lambda, size_t num_values)
{
    for (int i = 0; _vectorized_distribution.size() < num_values; ++i)
    {
        double pct = poisspdf(i, poisson_lambda);
        for (size_t j = 0; j < pct * num_values; ++j)
            _vectorized_distribution.push_back(i + 1);
        _frequency_percentage.push_back(pct);
    }
    // set a few extra percentages beyond the maximum size in the distribution
    for (int i = 0; i<5; ++i)
        _frequency_percentage.push_back(poisspdf(_frequency_percentage.size(), poisson_lambda));
}

void root_equilibrium_distribution::build_percentages()
{
    auto max = (size_t)*max_element(_vectorized_distribution.begin(), _vectorized_distribution.end()) + 1;
    _frequency_percentage.resize(max);
    for (size_t i = 0; i < max; ++i)
    {
        size_t c = count(_vectorized_distribution.begin(), _vectorized_distribution.end(), i);
        _frequency_percentage[i] = float(c) / float(_vectorized_distribution.size());
    }
}

float root_equilibrium_distribution::compute(const gene_transcript& t, size_t val) const
{
    if (_fixed_root_value > 0)
    {
        VectorXd v = VectorPos_bounds(_fixed_root_value, DISCRETIZATION_RANGE, bounds(t));
        return v[val];
    }
    if (val >= _frequency_percentage.size())
        return 0;

    return _frequency_percentage[val];
}

int root_equilibrium_distribution::select_root_size(int family_number) const
{
    if ((size_t)family_number >= _vectorized_distribution.size())
        return 0;

    return _vectorized_distribution[family_number];
}

void root_equilibrium_distribution::resize(size_t new_size)
{
    if (_fixed_root_value > 0)
        return;

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

TEST_CASE("Initializing with a fixed_root_value and resizing should not crash")
{
    root_equilibrium_distribution ef(5.3);
    ef.resize(9.2);
    CHECK(true);

}
TEST_CASE("Inference: root_equilibrium_distribution__with_no_rootdist_is_uniform")
{
    root_equilibrium_distribution ef(size_t(10));
    gene_transcript t;
    CHECK_EQ(doctest::Approx(.1), ef.compute(t, 5));
    CHECK_EQ(doctest::Approx(.1), ef.compute(t, 0));
}

TEST_CASE("root_equilibrium_distribution__resize")
{
    std::map<int, int> m;
    m[2] = 5;
    m[4] = 3;
    m[8] = 3;
    root_equilibrium_distribution rd(m);
    rd.resize(15);
    CHECK_EQ(rd.select_root_size(14), 8);
    CHECK_EQ(rd.select_root_size(15), 0);
}

TEST_CASE("root_equilibrium_distribution__poisson_compute")
{
    root_equilibrium_distribution pd(0.75, 100);
    gene_transcript t;

    CHECK_EQ(doctest::Approx(0.47059).scale(1000), pd.compute(t, 0));
    CHECK_EQ(doctest::Approx(0.13725f).scale(1000), pd.compute(t, 2));
    CHECK_EQ(doctest::Approx(0.005).scale(1000), pd.compute(t, 4));
    CHECK_EQ(0.0, pd.compute(t, 100));
}

TEST_CASE("root_equilibrium_distribution__poisson_select_root_size")
{
    root_equilibrium_distribution pd(0.75, 9);

    CHECK_EQ(1, pd.select_root_size(1));
    CHECK_EQ(1, pd.select_root_size(3));
    CHECK_EQ(2, pd.select_root_size(5));
    CHECK_EQ(2, pd.select_root_size(7));
    CHECK_EQ(0, pd.select_root_size(100));
}

TEST_CASE("root_equilibrium_distribution fixed_root_value")
{
    root_equilibrium_distribution pd(4.3);

    gene_transcript t;
    t.set_expression_value("A", 5.8);
    t.set_expression_value("B", 19.6);
    CHECK_EQ(doctest::Approx(6.055), pd.compute(t, 29));
    CHECK_EQ(doctest::Approx(0.713707), pd.compute(t, 30));
    CHECK_EQ(0.0, pd.compute(t, 100));

}

TEST_CASE("Simulation: uniform_distribution__select_root_size__returns_sequential_values")
{
    root_equilibrium_distribution ud(size_t(20));
    CHECK_EQ(1, ud.select_root_size(1));
    CHECK_EQ(2, ud.select_root_size(2));
    CHECK_EQ(3, ud.select_root_size(3));
    CHECK_EQ(4, ud.select_root_size(4));
    CHECK_EQ(5, ud.select_root_size(5));
    CHECK_EQ(0, ud.select_root_size(20));
}

