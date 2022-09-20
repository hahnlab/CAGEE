#include <iomanip>
#include <iostream>
#include <random>
#include <numeric>
#include <algorithm>
#include <iterator>

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
#include "proportional_variance.h"
#include "reconstruction.h"

#define GAMMA_INITIAL_GUESS_EXPONENTIAL_DISTRIBUTION_LAMBDA 1.75

extern std::mt19937 randomizer_engine;

using namespace std;
namespace pv = proportional_variance;

double compute_distribution_mean(const user_data& user_data)
{
    if (user_data.gene_families.empty()) throw runtime_error("No gene transcripts provided");
    if (!user_data.p_tree) throw runtime_error("No tree provided");

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

    double species_variance = std::accumulate(variances.begin(), variances.end(), 0.0) / double(variances.size());

    return species_variance / user_data.p_tree->distance_from_root_to_tip();
}

// sigma only
sigma_optimizer_scorer::sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma_squared* p_lambda) :
    _p_model(p_model), _user_data(user_data), _prior(prior), _p_sigma(p_lambda),
    optimize_sigma(true), optimize_epsilon(false), optimize_gamma(false)
{
}

// sigma and epsilon
sigma_optimizer_scorer::sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma_squared* p_lambda, error_model* p_error_model) :
    _p_model(p_model), _user_data(user_data), _prior(prior), _p_sigma(p_lambda),
    _p_error_model(p_error_model),
    optimize_sigma(true), optimize_epsilon(true), optimize_gamma(false)
{
}

// alpha and sigma
sigma_optimizer_scorer::sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma_squared* p_lambda) :
    _p_model(p_model), _user_data(user_data), _prior(prior), _p_sigma(p_lambda),
    optimize_sigma(true), optimize_epsilon(false), optimize_gamma(true)
{
}

// alpha only
sigma_optimizer_scorer::sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior) :
    _p_model(p_model), _user_data(user_data), _prior(prior), _p_sigma(nullptr),
    optimize_sigma(false), optimize_epsilon(false), optimize_gamma(true)
{
}

std::vector<double> sigma_optimizer_scorer::initial_guesses()
{
    std::vector<double> guesses;
    if (optimize_sigma)
    {
        if (!_distribution_mean)
            _distribution_mean = new double(compute_distribution_mean(_user_data));

        guesses.resize(_p_sigma->count());
        for (auto& i : guesses)
        {
            i = *_distribution_mean;
        }
    }

    if (optimize_epsilon)
    {
        current_epsilon_guesses = _p_error_model->get_epsilons();
        guesses.insert(guesses.end(), current_epsilon_guesses.begin(), current_epsilon_guesses.end());
    }

    if (optimize_gamma)
    {
        //It seems a gamma distribution with an alpha of 4 and a beta scaling factor of 0.25 works better. It has a mean of 1 and 70% of the density is between .5 and 1.5
        std::gamma_distribution<double> distribution(4.0, 0.25);
        guesses.push_back(distribution(randomizer_engine));
    }
    return guesses;
}

double sigma_optimizer_scorer::calculate_score(const double* values)
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

void sigma_optimizer_scorer::prepare_calculation(const double *values)
{
    int ptr = 0;
    if (optimize_sigma)
    {
        _p_sigma->update(values);
        ptr += _p_sigma->get_values().size();
    }
    if (optimize_epsilon)
    {
        auto epsilons = values + ptr;
        map<double, double> replacements;
        for (size_t i = 0; i < current_epsilon_guesses.size(); ++i)
        {
            replacements[current_epsilon_guesses[i]] = epsilons[i];
            current_epsilon_guesses[i] = epsilons[i];
        }

        _p_error_model->replace_epsilons(&replacements);
        ptr += current_epsilon_guesses.size();
    }
    if (optimize_gamma)
    {
        dynamic_cast<gamma_model *>(_p_model)->set_alpha(*(values + ptr));
    }
}

void sigma_optimizer_scorer::report_precalculation()
{
    ostringstream ost;
    if (optimize_sigma)
    {
        ost << "Sigma^2:" << *_p_sigma;
    }
    if (optimize_epsilon)
    {
        ost << " Epsilon:" << _p_error_model->get_epsilons().back() * 2.0;
    }
    if (optimize_gamma)
    {
        ost << " Alpha:" << dynamic_cast<gamma_model*>(_p_model)->get_alpha();
    }

    LOG(INFO) << ost.str();

}

void sigma_optimizer_scorer::finalize(double *results)
{
    _p_sigma->update(results);
    if (_p_error_model)
    {
        _p_error_model->update_single_epsilon(results[_p_sigma->count()]);
    }
}

std::string sigma_optimizer_scorer::description() const
{
    vector<string> t;
    ostringstream ost;
    if (optimize_sigma)
    {
        t.push_back("Sigma");
    }
    if (optimize_epsilon)
    {
        t.push_back("Epsilon");
    }
    if (optimize_gamma)
    {
        t.push_back("Alpha");
    }
    ost << "Optimizing ";
    std::ostream_iterator<string> out_it(ost, " ");
    std::copy(t.begin(), t.end(), out_it);
    return ost.str();
}

void sigma_optimizer_scorer::force_distribution_mean(double tree_length, double species_variance) {
    _distribution_mean = new double(sqrt(species_variance / tree_length));
}

TEST_CASE("sigma_optimizer_scorer constructor calculates tree length and variance")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    ud.p_tree = p_tree.get();
    sigma_squared s(5);
    sigma_optimizer_scorer soc((model *)nullptr, ud, std::gamma_distribution<double>(1,2), &s);

    auto guesses = soc.initial_guesses();
    REQUIRE(guesses.size() == 1);
    CHECK_EQ(doctest::Approx(0.05), guesses[0]);

    ud.gene_families.push_back(gene_transcript("TestFamily2", "", ""));
    ud.gene_families[1].set_expression_value("A", 5);
    ud.gene_families[1].set_expression_value("B", 8);
}

TEST_CASE("sigma_optimizer_scorer constructor averages variances across all transcripts")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families.push_back(gene_transcript("TestFamily2", "", ""));

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    ud.p_tree = p_tree.get();
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);
    ud.gene_families[1].set_expression_value("A", 5);
    ud.gene_families[1].set_expression_value("B", 8);

    sigma_squared s(5);
    sigma_optimizer_scorer soc((model*)nullptr, ud, std::gamma_distribution<double>(1, 2), &s);

    auto guesses = soc.initial_guesses();
    REQUIRE(guesses.size() == 1);
    CHECK_EQ(doctest::Approx(0.25), guesses[0]);
}


class mock_model : public model {
    // Inherited via model
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override { return nullptr; }
    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior) override { return nullptr; }
    bool _invalid_likelihood = false;
public:
    mock_model() : model(NULL, NULL, NULL) {}
    virtual double infer_family_likelihoods(const user_data& ud, const sigma_squared* p_lambda, const std::gamma_distribution<double>& prior) override
    { 
        return _invalid_likelihood ? nan("") : 0.0;
    }
    void set_invalid_likelihood() { _invalid_likelihood = true; }
};

TEST_CASE("lambda_epsilon_optimizer guesses lambda and unique epsilons")
{
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    sigma_squared s(10);
    user_data ud;
    mock_model model;
    sigma_optimizer_scorer leo(&model, ud, std::gamma_distribution<double>(1, 2), &s, &err);
    leo.force_distribution_mean(10, 1);
    auto guesses = leo.initial_guesses();
    REQUIRE(guesses.size() == 3);
    CHECK_EQ(doctest::Approx(0.31622).epsilon(0.00001), guesses[0]);
    CHECK_EQ(0.3, guesses[1]);
    CHECK_EQ(0.4, guesses[2]);
}

TEST_CASE("gamma_lambda_optimizer provides two guesses")
{
    sigma_squared sl(0.05);
    gamma_model model(NULL, NULL, 4, .25, NULL);
    user_data ud;
    sigma_optimizer_scorer glo(&model, ud, std::gamma_distribution<double>(1, 2), &sl);
    glo.force_distribution_mean(5, 1);
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
    sigma_optimizer_scorer optimizer(&m, ud, std::gamma_distribution<double>(1, 2));
    auto initial = optimizer.initial_guesses();
    CHECK_EQ(1, initial.size());
}


TEST_CASE("lambda_epsilon_optimizer")
{
    const double initial_epsilon = 0.01;
    error_model err;
    err.set_probabilities(0, { .0, .99, initial_epsilon });
    err.set_probabilities(1, { initial_epsilon, .98, initial_epsilon });

    mock_model model;

    sigma_squared lambda(0.05);
    user_data ud;
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(1, 2), &lambda, &err);
    optimizer.force_distribution_mean(10, 3);
    optimizer.initial_guesses();
    vector<double> values = { 0.05, 0.06 };
    optimizer.calculate_score(&values[0]);
    auto actual = err.get_probs(0);
    vector<double> expected{ 0, .94, .06 };
    CHECK(expected == actual);

    values[1] = 0.04;
    optimizer.calculate_score(&values[0]);
    actual = err.get_probs(0);
    expected = { 0, .96, .04 };
    CHECK(expected == actual);
}

TEST_CASE("prepare_calculation sets sigma correctly")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.p_tree = parse_newick("(A:1,B:1);");
    mock_model model;
    sigma_squared sig(1);
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(), &sig);
    vector<double> values{ 0.01 };
    optimizer.prepare_calculation(values.data());
    CHECK_EQ(0.01, sig.get_values()[0]);
}

TEST_CASE("prepare_calculation sets sigma and epsilon correctly ")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.p_tree = parse_newick("(A:1,B:1);");
    mock_model model;
    sigma_squared sig(1);
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(), &sig, &err);
    optimizer.initial_guesses();

    vector<double> values{ 0.01, 0.025 };
    optimizer.prepare_calculation(values.data());
    CHECK_EQ(0.01, sig.get_values()[0]);
    CHECK_EQ(0.025, err.get_epsilons()[0]);
}


TEST_CASE("prepare_calculation sets sigma and gamma correctly ")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.p_tree = parse_newick("(A:1,B:1);");
    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    gamma_model model(ud.p_lambda, &ud.gene_families, gamma_categories, multipliers, NULL);
    sigma_squared sig(1);
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(), &sig);
    optimizer.initial_guesses();

    vector<double> values{ 0.01, 0.25 };
    optimizer.prepare_calculation(values.data());
    CHECK_EQ(0.01, sig.get_values()[0]);
    CHECK_EQ(doctest::Approx(0.25), model.get_alpha());
}

#ifdef MODEL_GENE_EXPRESSION_LOGS
TEST_CASE("sigma_optimizer_scorer updates model alpha and lambda" * doctest::skip(true))
#else
TEST_CASE("sigma_optimizer_scorer updates model alpha and lambda")
#endif
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    ud.p_tree = parse_newick("(A:1,B:1);");
    sigma_squared sig(1);

    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    gamma_model m(&sig, &ud.gene_families, gamma_categories, multipliers, NULL);

    sigma_optimizer_scorer optimizer(&m, ud, std::gamma_distribution<double>(1, 2), &sig);
    optimizer.force_distribution_mean(7, 1);
    vector<double> values{ 0.01, 0.25 };
    optimizer.calculate_score(values.data());
    CHECK_EQ(doctest::Approx(0.25), m.get_alpha());
    CHECK_EQ(doctest::Approx(0.01), m.get_sigma()->get_values()[0]);
}

TEST_CASE("calculate_score translates nan to inf")
{
    sigma_squared lam(0.05);
    mock_model m;
    m.set_invalid_likelihood();
    double val;
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    ud.p_tree = parse_newick("(A:1,B:1);");
    sigma_optimizer_scorer opt(&m, ud, std::gamma_distribution<double>(1, 2), &lam);
    CHECK(std::isinf(opt.calculate_score(&val)));
}

