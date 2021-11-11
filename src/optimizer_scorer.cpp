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

// sigma only
sigma_optimizer_scorer::sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma* p_lambda) :
    inference_optimizer_scorer(p_lambda, p_model, user_data, prior),
    optimize_sigma(true), optimize_epsilon(false), optimize_gamma(false)
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

// sigma and epsilon
sigma_optimizer_scorer::sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma* p_lambda, error_model* p_error_model) :
    inference_optimizer_scorer(p_lambda, p_model, user_data, prior),
    _p_error_model(p_error_model),
    optimize_sigma(true), optimize_epsilon(true), optimize_gamma(false)
{
}

// alpha and sigma
sigma_optimizer_scorer::sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma* p_lambda) :
    inference_optimizer_scorer(p_lambda, p_model, user_data, prior),
    optimize_sigma(true), optimize_epsilon(false), optimize_gamma(true)
{
}

// alpha only
sigma_optimizer_scorer::sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior) :
    inference_optimizer_scorer(nullptr, p_model, user_data, prior),
    optimize_sigma(false), optimize_epsilon(false), optimize_gamma(true)
{
}

std::vector<double> sigma_optimizer_scorer::initial_guesses()
{
    std::vector<double> guesses;
    if (optimize_sigma)
    {
        double distmean = sqrt(_species_variance / _tree_length);
        guesses.resize(_p_sigma->count());
        std::normal_distribution<double> distribution(distmean, 0.2);
        for (auto& i : guesses)
        {
            i = distribution(randomizer_engine);
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

void sigma_optimizer_scorer::prepare_calculation(const double *values)
{
    int ptr = 0;
    if (optimize_sigma)
    {
        _p_sigma->update(values);
        ptr += _p_sigma->get_lambdas().size();
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
    sigma_optimizer_scorer soc((model *)nullptr, ud, std::gamma_distribution<double>(1,2), &s);

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
    sigma_optimizer_scorer soc((model*)nullptr, ud, std::gamma_distribution<double>(1, 2), &s);

    auto guesses = soc.initial_guesses();
    REQUIRE(guesses.size() == 1);
    CHECK_EQ(doctest::Approx(0.48974), guesses[0]);
}


class mock_model : public model {
    // Inherited via model
    virtual std::string name() const override { return "mockmodel"; }
    virtual void write_family_likelihoods(std::ostream& ost) override {}
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override { return nullptr; }
    virtual inference_optimizer_scorer* get_lambda_optimizer(const user_data& data, const std::gamma_distribution<double>& prior) override { return nullptr; }
    bool _invalid_likelihood = false;
public:
    mock_model() : model(NULL, NULL, NULL) {}
    virtual double infer_family_likelihoods(const user_data& ud, const sigma* p_lambda, const std::gamma_distribution<double>& prior) override 
    { 
        return _invalid_likelihood ? nan("") : 0.0;
    }
    void set_invalid_likelihood() { _invalid_likelihood = true; }
};

TEST_CASE("lambda_epsilon_optimizer guesses lambda and unique epsilons")
{
    randomizer_engine.seed(10);

    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    sigma s(10);
    user_data ud;
    mock_model model;
    sigma_optimizer_scorer leo(&model, ud, std::gamma_distribution<double>(1, 2), &s, &err);
    leo.force_tree_length_and_variance(10, 1);
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
    sigma_optimizer_scorer glo(&model, ud, std::gamma_distribution<double>(1, 2), &sl);
    glo.force_tree_length_and_variance(5, 1);
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

    sigma lambda(0.05);
    user_data ud;
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(1, 2), &lambda, &err);
    //sigma_optimizer_scorer optimizer(&model, &err, ud, std::gamma_distribution<double>(1, 2), &lambda, 10, 3);
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
    sigma sig(1);
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(), &sig);
    vector<double> values{ 0.01 };
    optimizer.prepare_calculation(values.data());
    CHECK_EQ(0.01, sig.get_lambdas()[0]);
}

TEST_CASE("prepare_calculation sets sigma and epsilon correctly ")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.p_tree = parse_newick("(A:1,B:1);");
    mock_model model;
    sigma sig(1);
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(), &sig, &err);
    optimizer.initial_guesses();

    vector<double> values{ 0.01, 0.025 };
    optimizer.prepare_calculation(values.data());
    CHECK_EQ(0.01, sig.get_lambdas()[0]);
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
    sigma sig(1);
    sigma_optimizer_scorer optimizer(&model, ud, std::gamma_distribution<double>(), &sig);
    optimizer.initial_guesses();

    vector<double> values{ 0.01, 0.25 };
    optimizer.prepare_calculation(values.data());
    CHECK_EQ(0.01, sig.get_lambdas()[0]);
    CHECK_EQ(doctest::Approx(0.25), model.get_alpha());
}

TEST_CASE("sigma_optimizer_scorer updates model alpha and lambda")
{
    user_data ud;
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    ud.p_tree = parse_newick("(A:1,B:1);");
    sigma sig(1);

    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    gamma_model m(&sig, &ud.gene_families, gamma_categories, multipliers, NULL);

    sigma_optimizer_scorer optimizer(&m, ud, std::gamma_distribution<double>(1, 2), &sig);
    optimizer.force_tree_length_and_variance(7, 1);
    vector<double> values{ 0.01, 0.25 };
    optimizer.calculate_score(values.data());
    CHECK_EQ(doctest::Approx(0.25), m.get_alpha());
    CHECK_EQ(doctest::Approx(0.01), m.get_lambda()->get_lambdas()[0]);
}

TEST_CASE("calculate_score translates nan to inf")
{
    sigma lam(0.05);
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

