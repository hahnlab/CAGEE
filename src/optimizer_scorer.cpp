#include <iomanip>
#include <iostream>
#include <random>
#include <numeric>
#include <algorithm>

#include "doctest.h"
#include "easylogging++.h"

#include "optimizer_scorer.h"
#include "clade.h"
#include "sigma.h"
#include "base_model.h"
#include "gamma_core.h"
#include "gamma.h"
#include "error_model.h"
#include "root_equilibrium_distribution.h"
#include "gene_transcript.h"
#include "user_data.h"

#define GAMMA_INITIAL_GUESS_EXPONENTIAL_DISTRIBUTION_LAMBDA 1.75

extern std::mt19937 randomizer_engine;

using namespace std;

double inference_optimizer_scorer::calculate_score(const double *values)
{
    prepare_calculation(values);

    if (!quiet)
    {
        report_precalculation();
    }

    double score = _p_model->infer_family_likelihoods(_user_data, _p_sigma, _prior);

    if (std::isnan(score)) score = -log(0);

    return score;
}

sigma_optimizer_scorer::sigma_optimizer_scorer(sigma* p_lambda, model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior) :
    inference_optimizer_scorer(p_lambda, p_model, user_data, prior)
{
    if (user_data.gene_families.empty()) throw runtime_error("No gene transcripts provided");
    if (!user_data.p_tree) throw runtime_error("No tree provided");

    _tree_length = user_data.p_tree->distance_from_root_to_tip();
    auto species = user_data.gene_families.at(0).get_species();
    vector<double> variances;
    for (auto& tt : user_data.gene_families)
    {
        vector<double> v(species.size());
        transform(species.begin(), species.end(), v.begin(), [&tt](string s) {return tt.get_expression_value(s);  });
        double sz = v.size();
        auto mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
        variances.push_back(std::accumulate(v.begin(), v.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
            return accumulator + ((val - mean) * (val - mean) / (sz - 1));
            }));

    }

    _species_variance = std::accumulate(variances.begin(), variances.end(), 0.0) / double(variances.size());
}

std::vector<double> sigma_optimizer_scorer::initial_guesses()
{
    double distmean = sqrt(_species_variance / _tree_length);
    std::vector<double> result(_p_sigma->count());
    std::normal_distribution<double> distribution(distmean,0.2);
    for (auto& i : result)
    {
    	i = distribution(randomizer_engine);
    }
    return result;
}

void sigma_optimizer_scorer::prepare_calculation(const double *values)
{
    _p_sigma->update(values);
}

void sigma_optimizer_scorer::report_precalculation()
{
    LOG(INFO) << "Sigma^2: " << *_p_sigma;
}

void sigma_optimizer_scorer::finalize(double *results)
{
    _p_sigma->update(results);
}

std::vector<double> lambda_epsilon_optimizer::initial_guesses()
{
    auto result = _lambda_optimizer.initial_guesses();

    current_guesses = _p_error_model->get_epsilons();
    result.insert(result.end(), current_guesses.begin(), current_guesses.end());

    return result;

}

void lambda_epsilon_optimizer::prepare_calculation(const double *values)
{
    auto lambdas = values;
    auto epsilons = values + _p_sigma->count();

    _lambda_optimizer.prepare_calculation(lambdas);

    map<double, double> replacements;
    for (size_t i = 0; i < current_guesses.size(); ++i)
    {
        replacements[current_guesses[i]] = epsilons[i];
        current_guesses[i] = epsilons[i];
    }

    _p_error_model->replace_epsilons(&replacements);
}

void lambda_epsilon_optimizer::report_precalculation()
{
    LOG(INFO) << "Calculating probability: epsilon=" << _p_error_model->get_epsilons().back()*2.0 << ", " << "lambda=" << *_p_sigma;
}

void lambda_epsilon_optimizer::finalize(double *results)
{
    _lambda_optimizer.finalize(results);
    _p_error_model->update_single_epsilon(results[_p_sigma->count()]);
}

gamma_optimizer::gamma_optimizer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior) :
    inference_optimizer_scorer(p_model->get_lambda(), p_model, user_data, prior),
    _p_gamma_model(p_model)
{

}
//In the first line that is commented out, Alpha is initiated by randomly drawing from an exponential distribution with a mean of 1.75
//It seems a gamma distribution with an alpha of 4 and a beta scaling factor of 0.25 works better. It has a mean of 1 and 70% of the density is between .5 and 1.5
std::vector<double> gamma_optimizer::initial_guesses()
{
    //std::exponential_distribution<double> distribution(GAMMA_INITIAL_GUESS_EXPONENTIAL_DISTRIBUTION_LAMBDA);
    std::gamma_distribution<double> distribution(4.0,0.25);
    return std::vector<double>({ distribution(randomizer_engine) });
}

void gamma_optimizer::prepare_calculation(const double * values)
{
    double alpha = *values;
    _p_gamma_model->set_alpha(alpha);
}

void gamma_optimizer::report_precalculation()
{
    LOG(INFO) << "Attempting alpha: " << _p_gamma_model->get_alpha();
}

void gamma_optimizer::finalize(double * result)
{
    _p_gamma_model->set_alpha(*result);
}

double gamma_optimizer::get_alpha() const
{
    return _p_gamma_model->get_alpha();
}

gamma_lambda_optimizer::gamma_lambda_optimizer(sigma*p_lambda, gamma_model * p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, double tree_length, double species_variance) :
    inference_optimizer_scorer(p_lambda, p_model, user_data, prior),
    _lambda_optimizer(p_lambda, p_model, user_data, prior, tree_length, species_variance),
    _gamma_optimizer(p_model, user_data, prior)
{
}

gamma_lambda_optimizer::gamma_lambda_optimizer(sigma* p_lambda, gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior) :
    inference_optimizer_scorer(p_lambda, p_model, user_data, prior),
    _lambda_optimizer(p_lambda, p_model, user_data, prior),
    _gamma_optimizer(p_model, user_data, prior)
{
}

std::vector<double> gamma_lambda_optimizer::initial_guesses()
{
    auto values = _lambda_optimizer.initial_guesses();
    auto alpha = _gamma_optimizer.initial_guesses();

    values.insert(values.end(), alpha.begin(), alpha.end());
    return values;

}

void gamma_lambda_optimizer::prepare_calculation(const double *values)
{
    _lambda_optimizer.prepare_calculation(values);
    _gamma_optimizer.prepare_calculation(values + _p_sigma->count());

}

void gamma_lambda_optimizer::report_precalculation()
{
    LOG(INFO) << "Attempting lambda: " << *_p_sigma << ", alpha: " << _gamma_optimizer.get_alpha();
}

/// results consists of the desired number of lambdas and one alpha value
void gamma_lambda_optimizer::finalize(double *results) {
    _lambda_optimizer.finalize(results);
    _gamma_optimizer.finalize(results + _p_sigma->count());
}

TEST_CASE("sigma_optimizer_scorer constructor calculates tree length and variance")
{
    randomizer_engine.seed(10);

    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    ud.p_tree = p_tree.get();
    sigma s(5);
    sigma_optimizer_scorer soc(&s, nullptr, ud, std::gamma_distribution<double>(1,2));

    auto guesses = soc.initial_guesses();
    REQUIRE(guesses.size() == 1);
    CHECK_EQ(doctest::Approx(0.213353), guesses[0]);

    ud.gene_families.push_back(gene_transcript("TestFamily2", "", ""));
    ud.gene_families[1].set_expression_value("A", 5);
    ud.gene_families[1].set_expression_value("B", 8);
}

TEST_CASE("sigma_optimizer_scorer constructor averages variances across all transcripts")
{
    randomizer_engine.seed(10);

    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families.push_back(gene_transcript("TestFamily2", "", ""));

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    ud.p_tree = p_tree.get();
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);
    ud.gene_families[1].set_expression_value("A", 5);
    ud.gene_families[1].set_expression_value("B", 8);

    sigma s(5);
    sigma_optimizer_scorer soc(&s, nullptr, ud, std::gamma_distribution<double>(1, 2));

    auto guesses = soc.initial_guesses();
    REQUIRE(guesses.size() == 1);
    CHECK_EQ(doctest::Approx(0.48974), guesses[0]);
}


TEST_CASE("lambda_epsilon_optimizer guesses lambda and unique epsilons")
{
    randomizer_engine.seed(10);

    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    sigma s(10);
    user_data ud;
    lambda_epsilon_optimizer leo(nullptr, &err, ud, std::gamma_distribution<double>(1, 2), &s, 10, 1);
    auto guesses = leo.initial_guesses();
    REQUIRE(guesses.size() == 3);
    CHECK_EQ(doctest::Approx(0.30597).epsilon(0.00001), guesses[0]);
    CHECK_EQ(0.3, guesses[1]);
    CHECK_EQ(0.4, guesses[2]);
}

TEST_CASE("gamma_lambda_optimizer provides two guesses")
{
    sigma sl(0.05);
    gamma_model model(NULL, NULL, 4, .25, NULL);
    user_data ud;
    gamma_lambda_optimizer glo(&sl, &model, ud, std::gamma_distribution<double>(1, 2), 5, 1);
    auto guesses = glo.initial_guesses();
    CHECK_EQ(2, guesses.size());

    double lambda = guesses[0];
    CHECK_GT(lambda, 0);
    CHECK_LT(lambda, 1);

    double alpha = guesses[1];
    CHECK_GT(alpha, 0);
    CHECK_LT(alpha, 10);
}

TEST_CASE("gamma_optimizer creates single initial guess")
{
    user_data ud;
    gamma_model m(NULL, NULL, 0, 0, NULL);
    gamma_optimizer optimizer(&m, ud, std::gamma_distribution<double>(1, 2));
    auto initial = optimizer.initial_guesses();
    CHECK_EQ(1, initial.size());
}



