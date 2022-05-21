#include <numeric>
#include <cmath>
#include <sstream>
#include <random>
#include <algorithm>

#include <string.h>

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include "src/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

#define DOCTEST_CONFIG_IMPLEMENT
#include "src/doctest.h"

#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"
#include "src/root_equilibrium_distribution.h"
#include "src/base_model.h"
#include "src/transcript_reconstructor.h"
#include "src/matrix_cache.h"
#include "src/probability.h"
#include "src/execute.h"
#include "src/user_data.h"
#include "src/optimizer_scorer.h"
#include "src/simulator.h"
#include "src/optimizer.h"
#include "src/error_model.h"
#include "src/likelihood_ratio.h"
#include "src/sigma.h"
#include "src/DiffMat.h"
#include "src/arguments.h"
#include "src/proportional_variance.h"

using namespace std;
using namespace Eigen;
namespace pv = proportional_variance;

#define NEED_TO_CONVERT_TO_MATRIX_ALGORITHM // note tests that are failing because they rely on the CAFE transition matrix algorithm

std::mt19937 randomizer_engine(10); // seeding random number engine

class mock_model : public model {
    // Inherited via model
    virtual std::string name() const override
    {
        return "mockmodel";
    }
    virtual void write_family_likelihoods(std::ostream& ost) override
    {
    }
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override
    {
        return nullptr;
    }
    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const vector<string>& sample_groups, const std::gamma_distribution<double>& prior) override
    {
        _p_sigma = initialize_search_sigma(data.p_lambda_tree, sample_groups);
        auto result = new sigma_optimizer_scorer(this, data, prior, _p_sigma);
        result->quiet = true;
        return result;
    }
    bool _invalid_likelihood = false;
public:
    mock_model() : model(NULL, NULL, NULL)
    {

    }
    void set_lambda(sigma_squared* lambda)
    {
        _p_sigma = lambda;
    }
    void set_invalid_likelihood() { _invalid_likelihood = true;  }

    // Inherited via model
    virtual double infer_family_likelihoods(const user_data& ud, const sigma_squared* p_lambda, const std::gamma_distribution<double>& prior) override
    {
        return _invalid_likelihood ? nan("") : 0.0;
    }
};

class Inference
{
public:
    user_data _user_data;
    sigma_squared* _p_lambda;

    Inference()
    {
        _p_lambda = new sigma_squared(0.05);
        _user_data.p_tree = parse_newick("(A:1,B:1);");
        _user_data.p_lambda = _p_lambda;
        _user_data.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
        _user_data.gene_families[0].set_expression_value("A", 1);
        _user_data.gene_families[0].set_expression_value("B", 2);
        //_user_data.p_prior = new root_equilibrium_distribution(size_t(100));
        //_user_data.p_prior = new root_distribution_uniform(100);
        randomizer_engine.seed(10);
    }

    ~Inference()
    {
        delete _user_data.p_tree;
        delete _p_lambda;
    }
};

class Optimizer
{
public:
    FMinSearch fm;
    vector<candidate> c;
    Optimizer()
    {
        for (int i = 0; i < 3; ++i)
            c.push_back(candidate(3));
        fm.variable_count = 2;
        fm.variable_count_plus_one = 3;
        fm.candidates.resize(3);
        transform(c.begin(), c.end(), fm.candidates.begin(), [](candidate& c) { return &c;  });
    }

};



#define STRCMP_EQUAL(x, y) CHECK(strcmp(x,y) == 0)
#define STRCMP_CONTAINS(x, y) CHECK(strstr(y,x) != nullptr)



TEST_CASE("GeneFamilies: species_size_is_case_insensitive")
{
    gene_transcript gf;
    gf.set_expression_value("Human", 5);
    CHECK_EQ(5, gf.get_expression_value("human"));
    CHECK_EQ(5, gf.get_expression_value("HUMAN"));
    CHECK_EQ(5, gf.get_expression_value("hUmAn"));
}

TEST_CASE("GeneFamilies: species_size_differential")
{
    gene_transcript gf;
    gf.set_expression_value("Cat", 5);
    gf.set_expression_value("Horse", 3);
    gf.set_expression_value("Cow", 1);

    CHECK_EQ(4, gf.species_size_differential());

    gf.set_expression_value("Chicken", 12);
    CHECK_EQ(11, gf.species_size_differential());
}

TEST_CASE_FIXTURE(Inference, "infer_processes"
    * doctest::skip(true))
{
    vector<gene_transcript> families;
    gene_transcript fam;
    fam.set_expression_value("A", 1);
    fam.set_expression_value("B", 2);
    families.push_back(fam);
    gene_transcript fam2;
    fam2.set_expression_value("A", 2);
    fam2.set_expression_value("B", 1);
    families.push_back(fam2);
    gene_transcript fam3;
    fam3.set_expression_value("A", 3);
    fam3.set_expression_value("B", 6);
    families.push_back(fam3);
    gene_transcript fam4;
    fam4.set_expression_value("A", 6);
    fam4.set_expression_value("B", 3);
    families.push_back(fam4);

    sigma_squared lambda(0.01);

    base_model core(&lambda, &families, NULL);

    double multi = core.infer_family_likelihoods(_user_data, &lambda, std::gamma_distribution<double>(1, 2));

    CHECK_EQ(doctest::Approx(46.56632), multi);
}

TEST_CASE("Inference: gamma_set_alpha")
{
    gamma_model model(NULL, NULL,  0, 0, NULL);
    model.set_alpha(0.5);
    CHECK(true);
}

TEST_CASE( "Inference: gamma_adjust_family_gamma_membership")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);

    gamma_model model(NULL, NULL, 0, 0, NULL);
    CHECK(true);
}

#ifdef MODEL_GENE_EXPRESSION_LOGS
TEST_CASE_FIXTURE(Inference, "gamma_model_infers_processes_without_crashing" * doctest::skip(true))
#else
TEST_CASE_FIXTURE(Inference, "gamma_model_infers_processes_without_crashing")
#endif
{
    gamma_model core(_user_data.p_lambda, &_user_data.gene_families, 1, 0, NULL);

    // TODO: make this return a non-infinite value and add a check for it
    core.infer_family_likelihoods(_user_data, _user_data.p_lambda, std::gamma_distribution<double>(1,2));
    CHECK(true);

}

TEST_CASE("Inference: stash_stream")
{
    family_info_stash stash;
    stash.family_id = "F01";
    stash.lambda_multiplier = 2.5;
    stash.family_likelihood = 3.7;
    stash.posterior_probability = 4.9;

    std::ostringstream ost;
    ost << stash;
    STRCMP_EQUAL("F01\t2.5\t0\t3.7\t4.9\tN/S", ost.str().c_str());

}

bool operator==(const MatrixXd& m1, const MatrixXd& m2)
{
    if (m1.size() != m2.size())
        return false;

    for (int i = 0; i < m1.size(); ++i)
    {
        for (int j = 0; j < m1.size(); ++j)
            if (abs(m1(i, j) - m2(i, j)) > 0.00001)
                return false;
    }

    return true;
}


TEST_CASE_FIXTURE(Inference, "base_model creates lambda_epsilon_optimizer if requested")
{
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    base_model model(_user_data.p_lambda, &_user_data.gene_families, &err);

    _user_data.p_error_model = nullptr;
    _user_data.p_lambda = nullptr;

    auto opt = model.get_sigma_optimizer(_user_data, vector<string>(), std::gamma_distribution<double>(1,2));

    REQUIRE(opt);
    CHECK_EQ("Optimizing Sigma Epsilon ", opt->description());
}

TEST_CASE_FIXTURE(Inference, "gamma_model_creates__gamma_lambda_optimizer_if_nothing_provided")
{
    gamma_model model(NULL, NULL,4, -1, NULL);
    _user_data.p_lambda = nullptr;

    auto opt = model.get_sigma_optimizer(_user_data, vector<string>(), std::gamma_distribution<double>(1, 2));
    REQUIRE(opt);
    CHECK_EQ("Optimizing Sigma Alpha ", opt->description());

    delete model.get_sigma();
}

TEST_CASE_FIXTURE(Inference, "base_model_reconstruction")
{
    sigma_squared sl(0.05);

    std::vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 4);
    boundaries b(0, 20);
    base_model model(&sl, &families, NULL);

    matrix_cache calc;
    root_distribution_uniform dist(size_t(10));

    std::unique_ptr<base_model_reconstruction> rec(dynamic_cast<base_model_reconstruction*>(model.reconstruct_ancestral_states(_user_data, &calc)));

    CHECK_EQ(1, rec->_reconstructions.size());

}

TEST_CASE("Inference: branch_length_finder")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));
    auto actual = p_tree->get_branch_lengths();
    auto expected = set<double>{ 1, 3, 7, 11, 17, 23 };
    CHECK(actual == expected);
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE_FIXTURE(Inference, "gamma_model_prune" * doctest::skip(true))
{
    vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    sigma_squared lambda(0.005);
    matrix_cache cache;

    std::gamma_distribution<double> prior(1,2);

    gamma_model model(&lambda, &families, { 0.01, 0.05 }, { 0.1, 0.5 }, NULL);

    vector<double> cat_likelihoods;
    CHECK(model.prune(families[0], prior, cache, &lambda, p_tree.get(), cat_likelihoods, 20));

    CHECK_EQ(2, cat_likelihoods.size());
    CHECK_EQ(doctest::Approx(-23.04433), log(cat_likelihoods[0]));
    CHECK_EQ(doctest::Approx(-16.68005), log(cat_likelihoods[1]));
}

TEST_CASE("gamma_model_prune_returns_false_if_saturated" * doctest::skip(true))
{
    vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    sigma_squared lambda(0.9);
    matrix_cache cache;

    vector<double> cat_likelihoods;

    gamma_model model(&lambda, &families, { 1.0,1.0 }, { 0.1, 0.5 }, NULL);

    CHECK(!model.prune(families[0], std::gamma_distribution<double>(1, 2), cache, &lambda, p_tree.get(), cat_likelihoods, 20));
}


TEST_CASE("Inference: create_one_model_if_lambda_is_null")
{
    input_parameters params;
    params.input_file_path = "foo";
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<base_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("Inference: create_gamma_model_if_alpha_provided")
{
    input_parameters params;
    params.input_file_path = "foo";
    params.fixed_alpha = 0.7;
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<gamma_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("Inference: create_gamma_model_if__n_gamma_cats__provided")
{
    input_parameters params;
    params.input_file_path = "foo";
    params.n_gamma_cats = 3;
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<gamma_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("Inference: build_models__uses_error_model_if_provided")
{
    input_parameters params;
    params.use_error_model = true;
    error_model em;
    em.set_probabilities(0, { 0, 0.99, 0.01 });
    user_data data;
    data.p_error_model = &em;
    data.p_tree = new clade("A", 5);
    sigma_squared lambda(0.05);
    data.p_lambda = &lambda;
    auto model = build_models(params, data)[0];
    std::ostringstream ost;
    model->write_vital_statistics(ost, data.p_tree, 0.07);
    STRCMP_CONTAINS("Epsilon: 0.01\n", ost.str().c_str());
    delete model;
}

TEST_CASE("Inference: build_models__creates_default_error_model_if_needed")
{
    input_parameters params;
    params.use_error_model = true;
    user_data data;
    data.p_tree = new clade("A", 5);
    sigma_squared lambda(0.05);
    data.p_lambda = &lambda;
    auto model = build_models(params, data)[0];
    std::ostringstream ost;
    model->write_vital_statistics(ost, data.p_tree, 0.01);
    STRCMP_CONTAINS("Epsilon: 0.05\n", ost.str().c_str());
    delete model;
}

#if 0
void build_matrix(Matrix3d& m)
{
    m(0, 0) = 1;
    m(0, 1) = 2;
    m(0, 2) = 3;
    m(1, 0) = 4;
    m(1, 1) = 5;
    m(1, 2) = 6;
    m(2, 0) = 7;
    m(2, 1) = 8;
    m(2, 2) = 9;
}

TEST_CASE("Probability, matrix_multiply")
{
    Matrix3d m1;
    build_matrix(m1);
    vector<double> m2({ 7, 9, 11 });
    vector<double> result;
    result.resize(3);
    m1.multiply(m2, 0, 2, 0, 2, result.data());
    CHECK_EQ(3, result.size());

    CHECK_EQ(58, result[0]);
    CHECK_EQ(139, result[1]);
    CHECK_EQ(220, result[2]);

    matrix m3(8);
    m3.set(3, 3, 1);
    m3.set(3, 4, 2);
    m3.set(3, 5, 3);
    m3.set(4, 3, 4);
    m3.set(4, 4, 5);
    m3.set(4, 5, 6);
    m3.set(5, 3, 7);
    m3.set(5, 4, 8);
    m3.set(5, 5, 9);

    m3.multiply(m2, 3, 5, 3, 5, result.data());
    CHECK_EQ(3, result.size());

    CHECK_EQ(58, result[0]);
    CHECK_EQ(139, result[1]);
    CHECK_EQ(220, result[2]);
}
#endif

TEST_CASE("Probability: error_model_set_probs")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    auto vec = model.get_probs(0);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.0, vec[0]);
    CHECK_EQ(0.7, vec[1]);
    CHECK_EQ(0.3, vec[2]);

    vec = model.get_probs(1);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);
}

TEST_CASE("Probability: error_model_get_epsilon")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(2, { .1, .8, .1 });
    model.set_probabilities(3, { .2, .6, .2 });

    auto actual = model.get_epsilons();
    CHECK_EQ(3, actual.size());
    vector<double> expected{ .1, .2, .3 };
    CHECK(expected == actual);
}

TEST_CASE("Probability: error_model_get_epsilon_zero_zero_must_be_zero")
{
    error_model model;
    CHECK_THROWS_WITH_AS(model.set_probabilities(0, { 0.4, 0.3, 0.3 }), "Cannot have a non-zero probability for family size 0 for negative deviation", runtime_error);
}

TEST_CASE("Probability: error_model__set_probabilities__cannot_set_higher_values_without_setting_zero")
{
    error_model model;
    CHECK_THROWS_WITH_AS(model.set_probabilities(5, { 0.4, 0.3, 0.3 }), "Cannot have a non-zero probability for family size 0 for negative deviation", runtime_error);
    model.set_probabilities(0, { 0, 0.7, 0.3 });
    model.set_probabilities(5, { 0.4, 0.3, 0.3 });
    CHECK(model.get_probs(5) == vector<double>({ 0.4, 0.3, 0.3 }));
}

TEST_CASE("Probability: error_model__set_probabilities__can_set_higher_values_if_valid_for_zero")
{
    error_model model;
    model.set_probabilities(5, { 0, 0.7, 0.3 });
    CHECK(model.get_probs(5) == vector<double>({ 0, 0.7, 0.3 }));
}

TEST_CASE("Probability: error_model_rows_must_add_to_one")
{
    error_model model;
    model.set_probabilities(0, { 0,1,0 });
    CHECK(true);
    CHECK_THROWS_WITH_AS(model.set_probabilities(1, { 0.3, 0.3, 0.3 }), "Sum of probabilities must be equal to one", runtime_error);
}

TEST_CASE("Probability: error_model_replace_epsilons")
{
    string input = "maxcnt: 10\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "2 0.2 0.6 0.2\n";

    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    map<double, double> replacements;
    replacements[.2] = .3;
    model.replace_epsilons(&replacements);

    auto actual = model.get_probs(0);
    CHECK_EQ(.7, actual[1]);
    CHECK_EQ(.3, actual[2]);

    actual = model.get_probs(1);
    CHECK_EQ(.3, actual[0]);
    CHECK_EQ(.4, actual[1]);
    CHECK_EQ(.3, actual[2]);
}

TEST_CASE("Probability: read_error_model")
{
    string input = "maxcnt: 10\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "2 0.2 0.6 0.2\n"
        "3 0.2 0.6 0.2\n"
        "5 0.2 0.6 0.2\n";

    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);
    auto vec = model.get_probs(0);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.0, vec[0]);
    CHECK_EQ(0.8, vec[1]);
    CHECK_EQ(0.2, vec[2]);

    vec = model.get_probs(1);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);

    vec = model.get_probs(4);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);

    vec = model.get_probs(7);
    CHECK_EQ(3, vec.size());
    CHECK_EQ(0.2, vec[0]);
    CHECK_EQ(0.6, vec[1]);
    CHECK_EQ(0.2, vec[2]);
}

TEST_CASE("Probability: write_error_model")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(2, { .1, .8, .1 });
    model.set_probabilities(3, { .2, .6, .2 });

    ostringstream ost;
    write_error_model_file(ost, model);

    const char* expected = "maxcnt: 3\ncntdiff: -1 0 1\n"
        "0 0 0.7 0.3\n"
        "1 0.2 0.6 0.2\n"
        "2 0.1 0.8 0.1\n"
        "3 0.2 0.6 0.2\n";

    STRCMP_EQUAL(expected, ost.str().c_str());
}

TEST_CASE("Probability: write_error_model_skips_unnecessary_lines")
{
    error_model model;
    model.set_probabilities(0, { .0, .7, .3 });
    model.set_probabilities(1, { .2, .6, .2 });
    model.set_probabilities(10, { .1, .8, .1 });
    model.set_probabilities(15, { .05, .9, .05 });
    model.set_probabilities(18, { .05, .9, .05 });

    ostringstream ost;
    write_error_model_file(ost, model);

    const char* expected = "maxcnt: 18\ncntdiff: -1 0 1\n"
        "0 0 0.7 0.3\n"
        "1 0.2 0.6 0.2\n"
        "10 0.1 0.8 0.1\n"
        "15 0.05 0.9 0.05\n";

    STRCMP_EQUAL(expected, ost.str().c_str());
}

TEST_CASE("Inference: build_reference_list")
{
    std::string str = "Desc\tFamily ID\tA\tB\n"
        "\t (null)1\t5\t10\n"
        "\t (null)2\t5\t7\n"
        "\t (null)3\t5\t10\n"
        "\t (null)4\t5\t7\n";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    auto actual = build_reference_list(families);
    vector<int> expected({ 0, 1, 0, 1 });
    CHECK_EQ(expected.size(), actual.size());

}


TEST_CASE("Clade: get_lambda_index_throws_from_branch_length_tree")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    CHECK_EQ(7, p_tree->get_branch_length());
    CHECK_THROWS_WITH_AS(p_tree->get_sigma_index(), "Requested sigma index from branch length tree", runtime_error);

}

TEST_CASE("Clade: get_branch_length_throws_from_lambda_tree")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7", true));
    CHECK_EQ(7, p_tree->get_sigma_index());
    CHECK_THROWS_WITH_AS(p_tree->get_branch_length(), "Requested branch length from sigma tree", runtime_error);

}

TEST_CASE("Clade: lambda_tree_root_index_is_1_if_not_specified")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:2)", true));
    CHECK_EQ(1, p_tree->get_sigma_index());
}

TEST_CASE("Clade: parse_newick_throws_exception_for_invalid_lambdas_in_tree")
{
    CHECK_THROWS_WITH_AS(parse_newick("(A:1,B:0):2", true), "Invalid sigma index set for B", runtime_error);
    CHECK_THROWS_WITH_AS(parse_newick("(A:-1,B:2)", true), "Invalid sigma index set for A", runtime_error);
}

TEST_CASE("Clade: parse_newick_throws_exception_for_invalid_branch_length_in_tree")
{
    CHECK_THROWS_WITH_AS(parse_newick("(A:1,B:0):2", false), "Invalid branch length set for B", runtime_error);
    CHECK_THROWS_WITH_AS(parse_newick("(A:-1,B:2)", false), "Invalid branch length set for A", runtime_error);
}

TEST_CASE("Clade: copy_constructor")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    unique_ptr<clade> copy(new clade(*p_tree.get()));
    CHECK_EQ(1, copy->find_descendant("A")->get_branch_length());
    CHECK_EQ(7, copy->find_descendant("AB")->get_branch_length());
}

TEST_CASE("Clade: copy_constructor_modifying_branch")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    unique_ptr<clade> copy(new clade(*p_tree.get(), nullptr, [](const clade& c) { return c.get_branch_length() * 2; }));
    CHECK_EQ(2, copy->find_descendant("A")->get_branch_length());
    CHECK_EQ(14, copy->find_descendant("AB")->get_branch_length());
}


TEST_CASE("Simulation: gamma_model_get_simulation_lambda_uses_multiplier_based_on_category_probability")
{
    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    sigma_squared lam(0.05);
    gamma_model m(&lam, NULL, gamma_categories, multipliers, NULL);
    vector<double> results(100);
    generate(results.begin(), results.end(), [&m]() {
        unique_ptr<sigma_squared> new_lam(dynamic_cast<sigma_squared*>(m.get_simulation_lambda()));
        return new_lam->get_values()[0];
        });

    CHECK_EQ(doctest::Approx(0.057), accumulate(results.begin(), results.end(), 0.0) / 100.0);

}

TEST_CASE("Inference: lambda_per_family" * doctest::skip(true))
{
    randomizer_engine.seed(10);
    user_data ud;

    ud.gene_families.resize(1);
   // ud.p_prior = new root_distribution_uniform(size_t(ud.max_root_family_size));

    gene_transcript& family = ud.gene_families[0];
    family.set_expression_value("A", 3);
    family.set_expression_value("B", 6);
    input_parameters params;
    params.lambda_per_family = true;

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    ud.p_tree = p_tree.get();
    estimator v(ud, params);

    mock_model m;
    ostringstream ost;
    v.estimate_lambda_per_family(&m, ost);
    CHECK_EQ(std::string("test\t0.30597463754818\n"), ost.str());
}

class mock_scorer : public optimizer_scorer
{
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        return std::vector<double>{0.2};
    }
    virtual double calculate_score(const double* values) override
    {
        if (force_scoring_error)
            return std::numeric_limits<double>::infinity();

        return 1000;
    }
public:
    bool force_scoring_error = false;
};

TEST_CASE("Inference, optimizer_gets_initial_guesses_from_scorer")
{
    mock_scorer scorer;
    optimizer opt(&scorer);
    auto guesses = opt.get_initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(0.2, guesses[0]);
}

TEST_CASE("Inference, optimizer_disallows_bad_initializations")
{
    mock_scorer scorer;
    scorer.force_scoring_error = true;
    optimizer opt(&scorer);

    CHECK_THROWS_WITH_AS(opt.get_initial_guesses(), "Failed to initialize any reasonable values", runtime_error);
}

TEST_CASE("Inference, optimizer_result_stream")
{
    optimizer::result r;
    r.num_iterations = 10;
    r.score = 5;
    r.values = vector<double>{ .05, .03 };
    r.duration = chrono::seconds(5000);

    ostringstream ost;
    ost << r;
    STRCMP_CONTAINS("Completed 10 iterations", ost.str().c_str());
    STRCMP_CONTAINS("Time: 1H 23M 20S", ost.str().c_str());
    STRCMP_CONTAINS("Best matches are: 0.05,0.03", ost.str().c_str());
    STRCMP_CONTAINS("Final -lnL: 5", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_no_attempts")
{
    event_monitor evm;

    ostringstream ost;
    evm.log(ost);
    STRCMP_EQUAL("No attempts made\n", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_one_attempt")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    ostringstream ost;

    evm.log(ost);
    STRCMP_EQUAL("1 values were attempted (0% rejected)\n", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_rejected_attempts")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_InvalidValues();
    ostringstream ost;

    evm.log(ost);
    STRCMP_EQUAL("2 values were attempted (50% rejected)\n", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_shows_poor_performing_families")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Saturation("test");
    ostringstream ost;

    evm.log(ost);
    STRCMP_CONTAINS("2 values were attempted (0% rejected)", ost.str().c_str());
    STRCMP_CONTAINS("The following families had failure rates >20% of the time:", ost.str().c_str());
    STRCMP_CONTAINS("test had 1 failures", ost.str().c_str());
}

TEST_CASE("Inference, event_monitor_does_not_show_decent_performing_families")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Saturation("test");
    ostringstream ost;

    evm.log(ost);
    STRCMP_EQUAL("5 values were attempted (0% rejected)\n", ost.str().c_str());
}

TEST_CASE_FIXTURE(Inference, "initialization_failure_advice_shows_20_families_with_largest_differentials")
{
    std::ostringstream ost;
    _user_data.gene_families.push_back(gene_transcript("TestFamily2", "", ""));
    _user_data.gene_families[1].set_expression_value("A", 34);
    _user_data.gene_families[1].set_expression_value("B", 86);

    initialization_failure_advice(ost, _user_data.gene_families);
    STRCMP_CONTAINS("Families with largest size differentials:", ost.str().c_str());
    STRCMP_CONTAINS("\nYou may want to try removing the top few families with the largest difference\nbetween the max and min counts and then re-run the analysis.\n", ost.str().c_str());
    STRCMP_CONTAINS("TestFamily2: 52\nTestFamily1: 1", ost.str().c_str());
}

TEST_CASE_FIXTURE(Optimizer, "fminsearch_sort_sorts_scores_and_moves_values")
{
    c[0].score = 3;
    c[1].score = 5;
    c[2].score = 1;

    __fminsearch_sort(&fm);

    CHECK_EQ(&c[2], fm.candidates[0]);
    CHECK_EQ(&c[0], fm.candidates[1]);
    CHECK_EQ(&c[1], fm.candidates[2]);
}

TEST_CASE_FIXTURE(Optimizer, "fminsearch_checkV_compares_value_difference_to_lx")
{
    c[0].values[0] = 1;
    c[1].values[0] = 2;
    c[2].values[0] = 3;
    c[0].values[1] = 3;
    c[1].values[1] = 4;
    c[2].values[1] = 5;
    fm.tolx = 3;
    CHECK(__fminsearch_checkV(&fm));
    fm.tolx = .5;
    CHECK_FALSE(__fminsearch_checkV(&fm));
}

TEST_CASE_FIXTURE(Optimizer, "fminsearch_checkF_compares_score_difference_to_lf")
{
    c[0].score = 1.0;
    c[1].score = 3.0;
    c[2].score = 5.0;

    fm.tolf = 5;
    CHECK(__fminsearch_checkF(&fm));
    fm.tolf = 1;
    CHECK_FALSE(__fminsearch_checkF(&fm));
}

class multiplier_scorer : public optimizer_scorer
{
    int num_scores;
public:
    multiplier_scorer() : multiplier_scorer(2) {}
    multiplier_scorer(int sz) : num_scores(sz)
    {

    }
    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override
    {
        vector<double> result(num_scores);
        result[0] = 5;
        result[1] = 3;
        return result;
    }
    virtual double calculate_score(const double* values) override
    {
        return std::accumulate(values, values + num_scores, 1.0, std::multiplies<double>());
    }
};

TEST_CASE_FIXTURE(Optimizer, "fminsearch_min_init")
{
    fm.delta = 0.05;
    fm.zero_delta = 0.00025;
    vector<int> indices(3);
    fm.idx = &indices[0];
    c[0].values[0] = 300;
    c[1].values[0] = 200;
    c[2].values[0] = 100;

    multiplier_scorer ms;
    fm.scorer = &ms;
    auto init = ms.initial_guesses();
    __fminsearch_min_init(&fm, &init[0]);

    CHECK_EQ(15, fm.candidates[0]->score);
    CHECK_EQ(15.75, fm.candidates[1]->score);
    CHECK_EQ(doctest::Approx(15.75), fm.candidates[2]->score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_mean")
{
    fm.variable_count = 2;
    vector<double> means(2);
    fm.x_mean = &means[0];
    fm.candidates[0]->values[0] = 300;
    fm.candidates[1]->values[0] = 200;
    fm.candidates[0]->values[1] = 12;
    fm.candidates[1]->values[1] = 44;

    __fminsearch_x_mean(&fm);

    CHECK_EQ(250, fm.x_mean[0]);
    CHECK_EQ(28, fm.x_mean[1]);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_reflection")
{
    fm.rho = 1;
    vector<double> means({ 250,28 });
    vector<double> reflections(2);
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.candidates[0]->values[0] = 300;
    fm.candidates[1]->values[0] = 200;
    fm.candidates[0]->values[1] = 12;
    fm.candidates[1]->values[1] = 44;

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_reflection(&fm);

    CHECK_EQ(500, fm.x_r[0]);
    CHECK_EQ(56, fm.x_r[1]);
    CHECK_EQ(28000, score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_expansion")
{
    fm.chi = 2;
    vector<double> means({ 250,28 });
    vector<double> expansions(2);
    vector<double> reflections({ 500, 56 });
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_expansion(&fm);

    CHECK_EQ(750, fm.x_tmp[0]);
    CHECK_EQ(84, fm.x_tmp[1]);
    CHECK_EQ(63000, score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_contract_outside")
{
    fm.psi = 0.5;
    vector<double> means({ 250,28 });
    vector<double> reflections({ 500, 56 });
    vector<double> expansions(2);
    fm.x_mean = &means[0];
    fm.x_r = &reflections[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_contract_outside(&fm);

    CHECK_EQ(375, fm.x_tmp[0]);
    CHECK_EQ(42, fm.x_tmp[1]);
    CHECK_EQ(15750, score);
}


TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_contract_inside")
{
    fm.psi = 0.5;
    vector<double> means({ 250,28 });
    vector<double> expansions(2);
    fm.candidates[2]->values[0] = 26;
    fm.candidates[2]->values[1] = 12;

    fm.x_mean = &means[0];
    fm.x_tmp = &expansions[0];

    multiplier_scorer ms;
    fm.scorer = &ms;

    double score = __fminsearch_x_contract_inside(&fm);

    CHECK_EQ(362, fm.x_tmp[0]);
    CHECK_EQ(36, fm.x_tmp[1]);
    CHECK_EQ(13032, score);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_x_shrink")
{
    fm.variable_count_plus_one = 3;
    fm.sigma = 0.5;
    vector<int> indices(3);
    fm.idx = &indices[0];
    vector<double> scores(3);

    fm.candidates[0]->values[0] = 300;
    fm.candidates[0]->values[1] = 200;
    fm.candidates[1]->values[0] = 42;
    fm.candidates[1]->values[1] = 64;
    fm.candidates[2]->values[0] = 26;
    fm.candidates[2]->values[1] = 12;

    multiplier_scorer ms;
    fm.scorer = &ms;
    __fminsearch_x_shrink(&fm);

    CHECK_EQ(300, fm.candidates[0]->values[0]);
    CHECK_EQ(200, fm.candidates[0]->values[1]);
    CHECK_EQ(163, fm.candidates[1]->values[0]);
    CHECK_EQ(106, fm.candidates[1]->values[1]);
    CHECK_EQ(171, fm.candidates[2]->values[0]);
    CHECK_EQ(132, fm.candidates[2]->values[1]);
}

TEST_CASE_FIXTURE(Optimizer, "__fminsearch_set_last_element")
{
    fm.variable_count_plus_one = 3;
    vector<int> indices(3);
    fm.idx = &indices[0];

    fm.candidates[0]->values[0] = 300;
    fm.candidates[0]->values[1] = 200;
    fm.candidates[0]->score = 2;
    fm.candidates[1]->values[0] = 42;
    fm.candidates[1]->values[1] = 64;
    fm.candidates[1]->score = 4;
    fm.candidates[2]->values[0] = 26;
    fm.candidates[2]->values[1] = 12;
    fm.candidates[2]->score = 6;

    vector<double> new_vals({ 99, 14 });
    __fminsearch_set_last_element(&fm, &new_vals[0], 3);

    CHECK_EQ(2.0, fm.candidates[0]->score);
    CHECK_EQ(3.0, fm.candidates[1]->score);
    CHECK_EQ(4.0, fm.candidates[2]->score);

}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_false_on_first_nine_attempts")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < 9; ++i)
        CHECK_FALSE(strat.threshold_achieved_checking_similarity(&fm));

}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_true_if_all_attempts_match")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < OPTIMIZER_SIMILARITY_CUTOFF_SIZE - 1; ++i)
        strat.threshold_achieved_checking_similarity(&fm);

    CHECK(strat.threshold_achieved_checking_similarity(&fm));
}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_false_if_attempt_varies")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < 9; ++i)
        strat.threshold_achieved_checking_similarity(&fm);

    fm.candidates[0]->score = 100.1;
    CHECK_FALSE(strat.threshold_achieved_checking_similarity(&fm));
}

TEST_CASE_FIXTURE(Optimizer, "NelderMeadSimilarityCutoff__threshold__returns_true_if_attempt_varies_minimally")
{
    NelderMeadSimilarityCutoff strat;
    FMinSearch fm;
    fm.variable_count = 1;
    candidate c(1), c2(1);
    c.values[0] = .0001;
    c.score = 100;
    c2.values[0] = .0005;
    c2.score = 101;
    fm.candidates = { &c, &c2 };
    for (int i = 0; i < OPTIMIZER_SIMILARITY_CUTOFF_SIZE - 1; ++i)
        strat.threshold_achieved_checking_similarity(&fm);

    fm.candidates[0]->score = 100.0001;
    CHECK(strat.threshold_achieved_checking_similarity(&fm));
}

TEST_CASE("LikelihoodRatioTest, update_branchlength")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    unique_ptr<clade> actual(LikelihoodRatioTest::update_branchlength(p_tree.get(), .5, 3));

    CHECK_EQ(15.5, actual->find_descendant("AB")->get_branch_length());
    CHECK_EQ(3.5, actual->find_descendant("A")->get_branch_length());
    CHECK_EQ(7.5, actual->find_descendant("B")->get_branch_length());
}

TEST_CASE("LikelihoodRatioTest, get_likelihood_for_diff_lambdas" * doctest::skip(true))
{
    mock_scorer s;
    optimizer opt(&s);
    gene_transcript gf;
    gf.set_expression_value("A", 5);
    gf.set_expression_value("B", 9);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    std::vector<sigma_squared*> cache(100);
    CHECK_EQ(0.0, LikelihoodRatioTest::get_likelihood_for_diff_lambdas(gf, p_tree.get(), 0, 0, cache, &opt, 20));
}

TEST_CASE("LikelihoodRatioTest, compute_for_diff_lambdas" * doctest::skip(true))
{
    sigma_squared lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    user_data data;
    data.p_lambda = &lam;
    data.p_tree = p_tree.get();
    data.gene_families.resize(1);
    data.gene_families[0].set_expression_value("A", 5);
    data.gene_families[0].set_expression_value("B", 9);
    vector<int> lambda_index(data.gene_families.size(), -1);
    vector<double> pvalues(data.gene_families.size());
    vector<sigma_squared*> lambdas(100);
    mock_scorer scorer;
    optimizer opt(&scorer);
    LikelihoodRatioTest::compute_for_diff_lambdas_i(data, lambda_index, pvalues, lambdas, &opt);
    CHECK_EQ(0, lambda_index[0]);
    CHECK(isinf(pvalues[0]));
}

void init_lgamma_cache();

int main(int argc, char** argv)
{
    init_lgamma_cache();

    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Global, el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", defaultConf);

    doctest::Context context;
    context.applyCommandLine(argc, argv);

    matrix_cache::initialize(DISCRETIZATION_RANGE);
    int res;
    try {
        res = context.run(); // run
    }
    catch (runtime_error& err) {
        LOG(ERROR) << err.what() << endl;
        res = EXIT_FAILURE;
    }

    if (context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests

    return res;
}
