#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <random>
#include <numeric>

#include "doctest.h"
#include "easylogging++.h"

#include "rootdist_estimator.h"
#include "clade.h"
#include "optimizer.h"
#include "gene_transcript.h"
#include "optimizer_scorer.h"

extern std::mt19937 randomizer_engine;

using namespace std;

double poisspdf(double x, double lambda)
{
  return exp(x*log(lambda) - lgamma(x + 1) - lambda);
}

double gammapdf(double value, const std::gamma_distribution<double>& dist) {
    double k = dist.alpha();
    double theta = dist.beta(); // see discussion in root_distribution_gamma

    double a = pow(theta, k);
    double b = pow(value, (k - 1));
    double c = exp(-1 * value / theta);
    double d = tgamma(k);

    return (b * c) / (a * d);
}

poisson_scorer::poisson_scorer(const vector<gene_transcript>& gene_families)
{
    for (auto &fam : gene_families)
    {
        for (auto species : fam.get_species())
            if (fam.get_expression_value(species) > 0)
                leaf_family_sizes.push_back(fam.get_expression_value(species) - 1);
    }
}

// Inherited via optimizer_scorer
std::vector<double> poisson_scorer::initial_guesses()
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double my_random = distribution(randomizer_engine);
    return std::vector<double>{my_random};
}

double poisson_scorer::calculate_score(const double * values)
{
    return lnLPoisson(values);

}

double poisson_scorer::lnLPoisson(const double* plambda)
{
    double lambda = plambda[0];
    if (lambda < 0)
        return -log(0);

    double score = accumulate(leaf_family_sizes.begin(), leaf_family_sizes.end(), 0.0, [lambda](double x, double sz) {
        double ll = poisspdf(sz, lambda);
        if (std::isnan(ll) || std::isinf(ll) || ll == 0) {
            return x;
        }
        return x + log(ll);
        });

    VLOG(3) << "Lambda: " << lambda << ", score " << -score << endl;
    return -score;
}

TEST_CASE("poisson_scorer_optimizes_correct_value")
{
    vector<gene_transcript> v;
    v.push_back(gene_transcript("TestFamily1", "", ""));
    v[0].set_expression_value("A", 1);
    v[0].set_expression_value("B", 2);

    poisson_scorer scorer(v);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    // DOUBLES_EQUAL(0.5, result.values[0], 0.0001)
}

TEST_CASE("poisson_scorer returns invalid for negative sigma")
{
    vector<gene_transcript> _;
    poisson_scorer scorer(_);
    double sigma = -1;
    double actual = scorer.lnLPoisson(&sigma);
    CHECK(std::isinf(actual));
}

TEST_CASE("poisson_scorer__lnlPoisson")
{
    vector<gene_transcript> v;
    v.push_back(gene_transcript("TestFamily1", "", ""));
    v[0].set_expression_value("A", 1);
    v[0].set_expression_value("B", 2);

    poisson_scorer scorer(v);
    double lambda = 0.05;
    CHECK_EQ(doctest::Approx(3.095732), scorer.lnLPoisson(&lambda));
}

TEST_CASE("poisson_scorer__lnlPoisson_skips_incalculable_family_sizes")
{
    vector<gene_transcript> v;
    v.push_back(gene_transcript("TestFamily1", "", ""));
    v[0].set_expression_value("A", 1);
    v[0].set_expression_value("B", 2);
    v.push_back(gene_transcript("TestFamily2", "", ""));
    v[1].set_expression_value("A", 3);
    v[1].set_expression_value("B", 175);

    poisson_scorer scorer(v);
    double lambda = 0.05;
    CHECK_EQ(doctest::Approx(9.830344), scorer.lnLPoisson(&lambda));
}

TEST_CASE("gammapdf function" * doctest::skip(true))
{
    const char *plotter = 
R"(#!/usr/bin/gnuplot
reset
set macros

set terminal pngcairo size 350, 262 enhanced font 'Verdana,10'
set output 'gammapdf.png'

set style line 1 lt 1 lc rgb '#FB9A99' # light red
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

set title "Gamma PDF alpha=0.75, beta=1.0/30.0"
set ylabel 'PDF Value' offset 1, 0

set nokey

plot 'gammapdf.dat' with lines ls 1)";

    std::ofstream p("gammapdf.sh");
    p << plotter;

    std::gamma_distribution<double> dist(0.75, 1.0/30.0);

    std::ofstream d("gammapdf.dat");
    for (double j = 0; j < 5; j += 0.01)
        d << j << " " << gammapdf(j, dist) << endl;
}