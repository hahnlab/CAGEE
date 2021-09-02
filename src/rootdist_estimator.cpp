#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <numeric>

#include "doctest.h"
#include "easylogging++.h"

#include "rootdist_estimator.h"
#include "clade.h"
#include "optimizer.h"
#include "probability.h"
#include "gene_transcript.h"
#include "optimizer_scorer.h"

extern std::mt19937 randomizer_engine;

using namespace std;

double poisspdf(double x, double lambda)
{
  return exp(x*log(lambda) - lgamma(x + 1) - lambda);
}

double gammapdf(double value, double alpha, double beta) {
    return (std::pow(beta, alpha) * std::pow(value, (alpha - 1)) * std::pow(M_E, (-1 * beta * value))) / tgamma(alpha);
}

vector<double> get_prior_rfsize_poisson_lambda(int min_family_size, int max_family_size, double poisson_lambda)
{
  int num_sizes = max_family_size - min_family_size;
  vector<double> prior_rfsize(num_sizes);

  for (int i = 0; i<num_sizes; i++) {

    //param->prior_rfsize[i] = poisspdf(param->pcafe->rootfamilysizes[0]+i, parameters[0]);					// poisson
    prior_rfsize[i] = poisspdf(i, poisson_lambda);					// shifted poisson
                                                   //param->prior_rfsize[i] = gampdf(param->pcafe->rootfamilysizes[0]+i, parameters[0], parameters[1]);	// gamma
  }
  return prior_rfsize;
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

gamma_scorer::gamma_scorer(const vector<gene_transcript>& gene_families)
{
    if (gene_families.empty()) return;
    if (gene_families[0].get_species().empty()) return;

    string species = gene_families[0].get_species()[0];
    for (auto& fam : gene_families)
    {
        leaf_family_sizes.push_back(fam.get_expression_value(species));
    }
}

std::vector<double> gamma_scorer::initial_guesses() 
{
    std::uniform_real_distribution<double> distribution(0.0, 10.0);
    return std::vector<double>{distribution(randomizer_engine), distribution(randomizer_engine)};
}

double gamma_scorer::calculate_score(const double* values)
{
    double alpha = values[0];
    double beta = values[1];

    double score = accumulate(leaf_family_sizes.begin(), leaf_family_sizes.end(), 0.0, [alpha, beta](double x, double sz) {
        double ll = gammapdf(sz, alpha, beta);
        if (std::isnan(ll) || std::isinf(ll) || ll == 0) {
            return x;
        }
        return x + log(ll);
        });

    VLOG(ROOT_SIZE_ESTIMATION) << "Gamma (" << alpha << "," << beta << ") score: " << -score << endl;
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

TEST_CASE("poisson_scorer returns invalid for negative lambda")
{
    vector<gene_transcript> _;
    poisson_scorer scorer(_);
    double lambda = -1;
    double actual = scorer.lnLPoisson(&lambda);
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


TEST_CASE("gamma_scorer optimizes correct value")
{
    vector<gene_transcript> v;
    v.push_back(gene_transcript("TestFamily1", "", ""));
    v[0].set_expression_value("A", 1);
    v[0].set_expression_value("B", 2);
    v.push_back(gene_transcript("TestFamily2", "", ""));
    v[1].set_expression_value("A", 3);
    v[1].set_expression_value("B", 7);

    gamma_scorer scorer(v);
    optimizer opt(&scorer);
    optimizer_parameters params;
    auto result = opt.optimize(params);

    REQUIRE_EQ(2, result.values.size());
    CHECK_EQ(doctest::Approx(3.6343), result.values[0]);
    CHECK_EQ(doctest::Approx(1.81715), result.values[1]);
}