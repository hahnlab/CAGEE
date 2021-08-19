#include <numeric>
#include <cmath>
#include <sstream>
#include <random>
#include <algorithm>

#include <string.h>

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define ELPP_NO_CHECK_MACROS
#include "src/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

#define DOCTEST_CONFIG_IMPLEMENT
#include "src/doctest.h"

#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"
#include "src/root_equilibrium_distribution.h"
#include "src/base_model.h"
#include "src/gene_family_reconstructor.h"
#include "src/matrix_cache.h"
#include "src/probability.h"
#include "src/execute.h"
#include "src/user_data.h"
#include "src/optimizer_scorer.h"
#include "src/simulator.h"
#include "src/poisson.h"
#include "src/optimizer.h"
#include "src/error_model.h"
#include "src/likelihood_ratio.h"
#include "src/lambda.h"
#include "src/DiffMat.h"
#include "src/arguments.h"

using namespace std;
using namespace Eigen;

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
    virtual inference_optimizer_scorer* get_lambda_optimizer(const user_data& data) override
    {
        initialize_lambda(data.p_lambda_tree);
        auto result = new sigma_optimizer_scorer(_p_lambda, this, data, 10, 1);
        result->quiet = true;
        return result;
    }
    bool _invalid_likelihood = false;
public:
    mock_model() : model(NULL, NULL, NULL)
    {

    }
    void set_lambda(lambda* lambda)
    {
        _p_lambda = lambda;
    }
    void set_invalid_likelihood() { _invalid_likelihood = true;  }

    // Inherited via model
    virtual double infer_family_likelihoods(const user_data& ud, const lambda* p_lambda) override
    {
        return _invalid_likelihood ? nan("") : 0.0;
    }
};

class Inference
{
public:
    user_data _user_data;
    single_lambda* _p_lambda;

    Inference()
    {
        _p_lambda = new single_lambda(0.05);
        _user_data.p_tree = parse_newick("(A:1,B:1);");
        _user_data.p_lambda = _p_lambda;
        _user_data.max_family_size = 10;
        _user_data.max_root_family_size = 8;
        _user_data.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
        _user_data.gene_families[0].set_expression_value("A", 1);
        _user_data.gene_families[0].set_expression_value("B", 2);

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

    single_lambda lambda(0.01);

    base_model core(&lambda, &families, NULL);

    double multi = core.infer_family_likelihoods(_user_data, &lambda);

    CHECK_EQ(doctest::Approx(46.56632), multi);
}

TEST_CASE_FIXTURE(Inference, "root_equilibrium_distribution__with_rootdist_uses_rootdist")
{
    _user_data.rootdist[1] = 3;
    _user_data.rootdist[2] = 5;

    root_equilibrium_distribution ef(_user_data.rootdist);
    gene_transcript t;

    CHECK_EQ(0.375, ef.compute(t, 1));
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

TEST_CASE_FIXTURE(Inference, "gamma_model_infers_processes_without_crashing")
{
    gamma_model core(_user_data.p_lambda, &_user_data.gene_families, 1, 0, NULL);

    // TODO: make this return a non-infinite value and add a check for it
    core.infer_family_likelihoods(_user_data, _user_data.p_lambda);
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

TEST_CASE("Probability:matrices_take_fractional_branch_lengths_into_account" * doctest::skip(true))
{
    single_lambda lambda(0.006335);
    matrix_cache calc(&lambda);
    std::set<double> branch_lengths{ 68, 68.7105 };
    calc.precalculate_matrices(set<boundaries>{boundaries(0,3)}, branch_lengths);
    CHECK_EQ(doctest::Approx(0.194661).epsilon(0.0001), calc.get_matrix(68.7105, boundaries(0,0.006335))(5, 5)); // a value 
    CHECK_EQ(doctest::Approx(0.195791).epsilon(0.0001), calc.get_matrix(68, boundaries(0, 0.006335))(5, 5));
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


TEST_CASE("Probability: probability_of_matrix" * doctest::skip(true))
{
    single_lambda lambda(0.05);
    matrix_cache calc(&lambda);
    std::set<double> branch_lengths{ 5 };
    calc.precalculate_matrices(set<boundaries>(), branch_lengths);
    auto actual = calc.get_matrix(5, boundaries());
    MatrixXd expected(5,5);
    double values[5][5] = {
    {1,0,0,0,0},
    { 0.2,0.64,0.128,0.0256,0.00512 },
    { 0.04,0.256,0.4608,0.17408,0.0512 },
    { 0.008,0.0768,0.26112,0.36352,0.187392 },
    { 0.0016,0.02048,0.1024,0.249856,0.305562 } };
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            expected(i, j) = values[i][j];
    CHECK(actual == expected);

    // a second call should get the same results as the first
    actual = calc.get_matrix(5, boundaries());
    CHECK(actual == expected);
}

TEST_CASE("Probability: get_random_probabilities" * doctest::skip(true))
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    single_lambda lam(0.05);
    matrix_cache cache(&lam);

    pvalue_parameters p = { p_tree.get(),  &lam, 12, 8, cache };
    auto probs = get_random_probabilities(p, 10, 3);
    CHECK_EQ(10, probs.size());
    CHECK_EQ(doctest::Approx(0.001905924).scale(10000), probs[0]);
}

TEST_CASE_FIXTURE(Inference, "base_optimizer_guesses_lambda_only")
{
    _user_data.p_lambda = NULL;

    base_model model(_user_data.p_lambda,  NULL, NULL);

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(_user_data));
    auto guesses = opt->initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(doctest::Approx(0.696853).epsilon(0.00001), guesses[0]);
    delete model.get_lambda();
}

TEST_CASE_FIXTURE(Inference, "base_model creates lambda_epsilon_optimizer if requested")
{
    error_model err;
    err.set_probabilities(0, { .0, .7, .3 });
    err.set_probabilities(1, { .4, .2, .4 });

    base_model model(_user_data.p_lambda, &_user_data.gene_families, &err);

    _user_data.p_error_model = nullptr;
    _user_data.p_lambda = nullptr;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(_user_data));

    CHECK(opt);
    CHECK(dynamic_cast<lambda_epsilon_optimizer*>(opt.get()) != nullptr);
}

TEST_CASE_FIXTURE(Inference, "gamma_model_creates__gamma_lambda_optimizer_if_nothing_provided")
{
    gamma_model model(NULL, NULL,4, -1, NULL);
    _user_data.p_lambda = nullptr;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(_user_data));
    REQUIRE(opt);
    CHECK(dynamic_cast<gamma_lambda_optimizer*>(opt.get()));

    delete model.get_lambda();
}

TEST_CASE("Inference: gamma_model__creates__lambda_optimizer__if_alpha_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, NULL, 4, 0.25, NULL);

    user_data data;
    data.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    data.gene_families[0].set_expression_value("A", 1);
    data.gene_families[0].set_expression_value("B", 2);
    data.p_tree = p_tree.get();

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));

    CHECK(opt);
    CHECK(dynamic_cast<sigma_optimizer_scorer*>(opt.get()));
    delete model.get_lambda();
}

TEST_CASE("Inference: gamma_model__creates__gamma_optimizer__if_lambda_provided")
{
    gamma_model model(NULL, NULL, 4, -1, NULL);

    user_data data;

    single_lambda sl(0.05);
    data.p_lambda = &sl;

    unique_ptr<inference_optimizer_scorer> opt(model.get_lambda_optimizer(data));

    CHECK(opt);
    CHECK(dynamic_cast<gamma_optimizer*>(opt.get()));

    delete model.get_lambda();
}

TEST_CASE("Inference: gamma_model_creates_nothing_if_lambda_and_alpha_provided")
{
    gamma_model model(NULL, NULL, 4, .25, NULL);

    user_data data;

    single_lambda sl(0.05);
    data.p_lambda = &sl;

    CHECK(model.get_lambda_optimizer(data) == nullptr);
}

TEST_CASE("Inference: gamma_lambda_optimizer__provides_two_guesses")
{
    single_lambda sl(0.05);
    gamma_model model(NULL, NULL, 4, .25, NULL);
    user_data ud;
    gamma_lambda_optimizer glo(&sl, &model, ud, 5, 1);
    auto guesses = glo.initial_guesses();
    CHECK_EQ(2, guesses.size());

    double lambda = guesses[0];
    CHECK_GT(lambda, 0);
    CHECK_LT(lambda, 1);

    double alpha = guesses[1];
    CHECK_GT(alpha, 0);
    CHECK_LT(alpha, 10);
}

TEST_CASE("Inference: gamma_optimizer__creates_single_initial_guess")
{
    user_data ud;
    gamma_model m(NULL, NULL, 0, 0, NULL);
    gamma_optimizer optimizer(&m, ud);
    auto initial = optimizer.initial_guesses();
    CHECK_EQ(1, initial.size());
}

TEST_CASE_FIXTURE(Inference, "base_model_reconstruction")
{
    single_lambda sl(0.05);

    std::vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 4);

    base_model model(&sl, &families, NULL);

    matrix_cache calc(&sl);
    calc.precalculate_matrices(set<boundaries>(), set<double>({ 1 }));
    root_equilibrium_distribution dist(size_t(_user_data.max_root_family_size));

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

TEST_CASE("Inference: increase_decrease")
{
    base_model_reconstruction bmr;
    gene_transcript gf("myid", "", "");

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

    auto a = p_tree->find_descendant("A");
    auto b = p_tree->find_descendant("B");
    auto ab = p_tree->find_descendant("AB");
    auto abcd = p_tree->find_descendant("ABCD");

    gf.set_expression_value("A", 4);
    gf.set_expression_value("B", 2);
    bmr._reconstructions["myid"][ab] = 3;
    bmr._reconstructions["myid"][abcd] = 3;

    CHECK_EQ(1, bmr.get_difference_from_parent(gf, a));
    CHECK_EQ(-1, bmr.get_difference_from_parent(gf, b));
    CHECK_EQ(0, bmr.get_difference_from_parent(gf, ab));
}

TEST_CASE( "Inference: precalculate_matrices_calculates_all_boundaries_all_branchlengths")
{
    single_lambda s(0.05);
    matrix_cache calc(&s);
    set<boundaries> b{ boundaries(0,0.5), boundaries(0, 0.25), boundaries(0, 0.1), boundaries(0, 0.2) };
    calc.precalculate_matrices(b, set<double>({ 1,2,3 }));
    CHECK_EQ(12, calc.get_cache_size());
}

class Reconstruction
{
public:
    gene_transcript fam;
    unique_ptr<clade> p_tree;
    cladevector order;

    Reconstruction() : fam("Family5", "", "")
    {
        p_tree.reset(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

        fam.set_expression_value("A", 11);
        fam.set_expression_value("B", 2);
        fam.set_expression_value("C", 5);
        fam.set_expression_value("D", 6);

        vector<string> nodes{ "A", "B", "C", "D", "AB", "CD", "ABCD" };
        order.resize(nodes.size());
        const clade* t = p_tree.get();
        transform(nodes.begin(), nodes.end(), order.begin(), [t](string s) { return t->find_descendant(s); });

    }
};

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_leaf_node" * doctest::skip(true))
{
    single_lambda lambda(0.1);
    fam.set_expression_value("Mouse", 3);

    clade leaf("Mouse", 7);

    matrix_cache calc(&lambda);
    calc.precalculate_matrices(set<boundaries>(), set<double>({ 7 }));
    clademap<std::vector<double>> all_node_Cs;
    clademap<std::vector<double>> all_node_Ls;
    all_node_Cs[&leaf].resize(8);
    all_node_Ls[&leaf].resize(8);

    pupko_reconstructor::reconstruct_leaf_node(&leaf, fam, all_node_Cs, all_node_Ls, &calc);

    // L holds the probability of the leaf moving from size 3 to size n
    auto L = all_node_Ls[&leaf];

    CHECK_EQ(8, L.size());
    CHECK_EQ(0.0, L[0]);
    CHECK_EQ(doctest::Approx(0.0586679), L[1]);
    CHECK_EQ(doctest::Approx(0.146916), L[2]);
    CHECK_EQ(doctest::Approx(0.193072), L[3]);
}

TEST_CASE_FIXTURE(Reconstruction, "print_reconstructed_states__prints_star_for_significant_values")
{
    gene_transcript gf("Family5", "", "");
    base_model_reconstruction bmr;
    auto& values = bmr._reconstructions[gf.id()];

    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    branch_probabilities branch_probs;
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&branch_probs, &gf](const clade* c) {branch_probs.set(gf, c, branch_probabilities::branch_probability(.5)); });

    branch_probs.set(gf, p_tree->find_descendant("AB"), 0.02);
    branch_probs.set(gf, p_tree.get(), branch_probabilities::invalid());  /// root is never significant regardless of the value

    ostringstream sig;
    bmr.print_reconstructed_states(sig, order, { fam }, p_tree.get(), 0.05, branch_probs);
    STRCMP_CONTAINS("((A<0>_11.000000:1,B<1>_2.000000:3)<4>*_8.000000:7,(C<2>_5.000000:11,D<3>_6.000000:17)<5>_6.000000:23)<6>_7.000000", sig.str().c_str());

    ostringstream insig;
    bmr.print_reconstructed_states(insig, order, { fam }, p_tree.get(), 0.01, branch_probs);
    STRCMP_CONTAINS("  TREE Family5 = ((A<0>_11.000000:1,B<1>_2.000000:3)<4>_8.000000:7,(C<2>_5.000000:11,D<3>_6.000000:17)<5>_6.000000:23)<6>_7.000000;", insig.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__print_reconstructed_states__prints_value_for_each_category_and_a_summation")
{
    gamma_model_reconstruction gmr(vector<double>({ 1.0 }));

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()] = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")] = 0;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")] = 0;

    rec.reconstruction[p_tree.get()] = 7;
    rec.reconstruction[p_tree->find_descendant("AB")] = 8;
    rec.reconstruction[p_tree->find_descendant("CD")] = 6;

    ostringstream ost;
    branch_probabilities branch_probs;
    gmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);
    STRCMP_CONTAINS("  TREE Family5 = ((A<0>_11.000000:1,B<1>_2.000000:3)<4>_8.000000:7,(C<2>_5.000000:11,D<3>_6.000000:17)<5>_6.000000:23)<6>_7.000000;", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__print_additional_data__prints_likelihoods")
{
    gamma_model_reconstruction gmr(vector<double>({ 0.3, 0.9, 1.4, 2.0 }));
    gmr._reconstructions["Family5"]._category_likelihoods = { 0.01, 0.03, 0.09, 0.07 };
    ostringstream ost;
    gmr.print_category_likelihoods(ost, order, { fam });
    STRCMP_CONTAINS("Family ID\t0.3\t0.9\t1.4\t2\t\n", ost.str().c_str());
    STRCMP_CONTAINS("Family5\t0.01\t0.03\0.09\t0.07", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__prints_lambda_multipiers")
{
    vector<double> multipliers{ 0.13, 1.4 };
    gamma_model_reconstruction gmr(multipliers);

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()] = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")] = 8;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")] = 6;

    rec.reconstruction[p_tree.get()] = 7;
    rec.reconstruction[p_tree->find_descendant("AB")] = 8;
    rec.reconstruction[p_tree->find_descendant("CD")] = 6;

    branch_probabilities branch_probs;

    std::ostringstream ost;
    gmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);

    STRCMP_CONTAINS("BEGIN LAMBDA_MULTIPLIERS;", ost.str().c_str());
    STRCMP_CONTAINS("  0.13;", ost.str().c_str());
    STRCMP_CONTAINS("  1.4;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "base_model_reconstruction__print_reconstructed_states")
{
    base_model_reconstruction bmr;
    auto& values = bmr._reconstructions["Family5"];

    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    branch_probabilities branch_probs;

    ostringstream ost;

    bmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);
    STRCMP_CONTAINS("#nexus", ost.str().c_str());
    STRCMP_CONTAINS("BEGIN TREES;", ost.str().c_str());
    CHECK_MESSAGE(ost.str().find("  TREE Family5 = ((A<0>_11.000000:1,B<1>_2.000000:3)<4>_8.000000:7,(C<2>_5.000000:11,D<3>_6.000000:17)<5>_6.000000:23)<6>_7.000000;") != string::npos, ost.str());
    STRCMP_CONTAINS("END;", ost.str().c_str());

}

TEST_CASE_FIXTURE(Reconstruction, "reconstruction_process_internal_node" * doctest::skip(true))
{
    single_lambda s_lambda(0.1);
    fam.set_expression_value("A", 3);
    fam.set_expression_value("B", 6);

    matrix_cache calc(&s_lambda);
    calc.precalculate_matrices(set<boundaries>(), set<double>({ 1, 3, 7, 11, 17, 23 }));

    clademap<std::vector<double>> all_node_Cs;
    clademap<std::vector<double>> all_node_Ls;
    all_node_Cs[p_tree->find_descendant("A")].resize(25);
    all_node_Ls[p_tree->find_descendant("A")].resize(25);
    all_node_Cs[p_tree->find_descendant("B")].resize(25);
    all_node_Ls[p_tree->find_descendant("B")].resize(25);
    all_node_Cs[p_tree->find_descendant("AB")].resize(25);
    all_node_Ls[p_tree->find_descendant("AB")].resize(25);

    pupko_reconstructor::reconstruct_leaf_node(p_tree->find_descendant("A"), fam, all_node_Cs, all_node_Ls, &calc);
    pupko_reconstructor::reconstruct_leaf_node(p_tree->find_descendant("B"), fam, all_node_Cs, all_node_Ls, &calc);

    auto internal_node = p_tree->find_descendant("AB");
    pupko_reconstructor::reconstruct_internal_node(internal_node, fam, all_node_Cs, all_node_Ls, &calc);
    auto L = all_node_Ls[internal_node];

    // L holds the probability of the node moving from size 3 to size n
    CHECK_EQ(25, L.size());
    CHECK_EQ(0.0, L[0]);
    CHECK_EQ(doctest::Approx(0.00101688), L[1]);
    CHECK_EQ(doctest::Approx(0.00254648), L[2]);
    CHECK_EQ(doctest::Approx(0.0033465), L[3]);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruction_process_internal_node with 0s at the leafs" * doctest::skip(true))
{
    single_lambda s_lambda(0.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc(&s_lambda);
    calc.precalculate_matrices(set<boundaries>(), set<double>({ 1, 3, 7, 11, 17, 23 }));

    clademap<std::vector<double>> all_node_Cs;
    clademap<std::vector<double>> all_node_Ls;
    all_node_Cs[p_tree->find_descendant("A")].resize(25);
    all_node_Ls[p_tree->find_descendant("A")].resize(25);
    all_node_Cs[p_tree->find_descendant("B")].resize(25);
    all_node_Ls[p_tree->find_descendant("B")].resize(25);
    all_node_Cs[p_tree->find_descendant("AB")].resize(25);
    all_node_Ls[p_tree->find_descendant("AB")].resize(25);

    pupko_reconstructor::reconstruct_leaf_node(p_tree->find_descendant("A"), fam, all_node_Cs, all_node_Ls, &calc);
    pupko_reconstructor::reconstruct_leaf_node(p_tree->find_descendant("B"), fam, all_node_Cs, all_node_Ls, &calc);

    auto internal_node = p_tree->find_descendant("AB");
    pupko_reconstructor::reconstruct_internal_node(internal_node, fam, all_node_Cs, all_node_Ls, &calc);
    auto L = all_node_Ls[internal_node];

    // L holds the probability of the node moving from size 3 to size n
    CHECK_EQ(25, L.size());
    CHECK_EQ(1.0, L[0]);
    CHECK_EQ(doctest::Approx(0.4117647059), L[1]);
    CHECK_EQ(doctest::Approx(0.169550173), L[2]);
    CHECK_EQ(doctest::Approx(0.0698147771), L[3]);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript" * doctest::skip(true))
{
    gene_transcript fam;
    fam.set_expression_value("A", 3);
    fam.set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    single_lambda lambda(0.005);
    matrix_cache cache(&lambda);

    cache.precalculate_matrices(set<boundaries>(), set<double>{1, 3, 7});

    user_data ud;
    ud.max_root_family_size = 8;

    clademap<int> result;
    clademap<std::vector<double>> all_node_Cs;
    clademap<std::vector<double>> all_node_Ls;
    std::function <void(const clade*)> pupko_initializer = [this, &all_node_Cs, &all_node_Ls, &ud](const clade* c) {
        pupko_reconstructor::initialize_at_node(c, all_node_Cs, all_node_Ls, 10, ud.max_root_family_size);
    };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), pupko_initializer);
    pupko_reconstructor::reconstruct_gene_transcript(&lambda, p_tree.get(), &fam, &cache, result, all_node_Cs, all_node_Ls);
    auto AB = p_tree->find_descendant("AB");
    CHECK_EQ(4, result[AB]);
}

TEST_CASE_FIXTURE(Reconstruction, "get_weighted_averages")
{
    clade c1;
    clade c2;

    clademap<int> rc1;
    rc1[&c1] = 10;
    rc1[&c2] = 2;

    clademap<int> rc2;
    rc2[&c1] = 20;
    rc2[&c2] = 8;

    auto avg = get_weighted_averages({ rc1, rc2 }, { .25, .75 });
    CHECK_EQ(17.5, avg[&c1]);
    CHECK_EQ(6.5, avg[&c2]);
}

TEST_CASE_FIXTURE(Reconstruction, "print_node_counts")
{
    gamma_model_reconstruction gmr({ .5 });
    ostringstream ost;

    auto initializer = [&gmr](const clade* c) { gmr._reconstructions["Family5"].reconstruction[c] = 5;  };
    p_tree->apply_prefix_order(initializer);

    gmr.print_node_counts(ost, order, { fam }, p_tree.get());
    STRCMP_CONTAINS("FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>", ost.str().c_str());
    CHECK_MESSAGE(ost.str().find("Family5\t11.000000\t2.000000\t5.000000\t6.000000\t5.000000\t5.000000\t5.000000") != string::npos, ost.str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_node_change")
{
    gamma_model_reconstruction gmr({ .5 });
    ostringstream ost;

    std::normal_distribution<float> dist(0, 10);
    clademap<int> size_deltas;
    p_tree->apply_prefix_order([&gmr, &dist](const clade* c) {
        if (!c->is_leaf())
            gmr._reconstructions["Family5"].reconstruction[c] = dist(randomizer_engine);
        });

    gmr.print_node_change(ost, order, { fam }, p_tree.get());
    CHECK_MESSAGE(ost.str().find("FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>") != string::npos, ost.str());
    CHECK_MESSAGE(ost.str().find("Family5\t+1\t-8\t+5\t+6\t+17\t+7\t+0") != string::npos, ost.str());
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_index_in_angle_brackets_for_non_leaf")
{
    STRCMP_EQUAL("<0>", clade_index_or_name(p_tree.get(), { p_tree.get() }).c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_name_plus_index_in_angle_brackets_for_leaf")
{
    auto a = p_tree->find_descendant("A");
    STRCMP_EQUAL("A<1>", clade_index_or_name(a, { p_tree.get(), a }).c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_branch_probabilities__shows_NA_for_invalids")
{
    gene_transcript gf("Family5", "", "");
    std::ostringstream ost;
    branch_probabilities probs;
    for (auto c : order)
        probs.set(gf, c, 0.05);
    probs.set(gf, p_tree->find_descendant("B"), branch_probabilities::invalid());
    probs.set(gf, p_tree.get(), branch_probabilities::invalid());

    print_branch_probabilities(ost, order, { fam }, probs);
    STRCMP_CONTAINS("FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>", ost.str().c_str());
    STRCMP_CONTAINS("Family5\t0.05\tN/A\t0.05\t0.05\t0.05\t0.05\tN/A\n", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_branch_probabilities__skips_families_without_reconstructions")
{
    std::ostringstream ost;
    branch_probabilities probs;

    print_branch_probabilities(ost, order, { fam }, probs);
    CHECK(ost.str().find("Family5") == string::npos);
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities" * doctest::skip(true))
{
    single_lambda lm(0.05);
    matrix_cache cache(&lm);
    cache.precalculate_matrices(set<boundaries>(), { 1,3,7 });
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree->find_descendant("AB")] = 10;
    rec._reconstructions[fam.id()][p_tree->find_descendant("ABCD")] = 12;
    CHECK_EQ(doctest::Approx(0.537681), compute_viterbi_sum(p_tree->find_descendant("AB"), fam, &rec, 24, cache, &lm)._value);
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities_returns_invalid_if_root")
{
    single_lambda lm(0.05);
    matrix_cache cache(&lm);
    cache.precalculate_matrices(set<boundaries>(), { 1,3,7 });
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree.get()] = 11;
    CHECK_FALSE(compute_viterbi_sum(p_tree.get(), fam, &rec, 24, cache, &lm)._is_valid);
}

TEST_CASE_FIXTURE(Reconstruction, "pvalues")
{
    vector<double> cd(10);
    double n = 0;
    std::generate(cd.begin(), cd.end(), [&n]() mutable { return n += 0.01; });
    CHECK_EQ(0.5, pvalue(0.05, cd));
    CHECK_EQ(0.0, pvalue(0.0001, cd));
    CHECK_EQ(0.9, pvalue(0.099, cd));
}

TEST_CASE_FIXTURE(Reconstruction, "pvalues 2")
{
    double v = .35;
    std::vector<double> conddist = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };
    double actual = pvalue(v, conddist);
    CHECK_EQ(doctest::Approx(3.0 / 9.0), actual);
}

TEST_CASE_FIXTURE(Inference, "gamma_model_prune" * doctest::skip(true))
{
    vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    single_lambda lambda(0.005);
    matrix_cache cache(&lambda);

    _user_data.rootdist[1] = 2;
    _user_data.rootdist[2] = 2;
    _user_data.rootdist[3] = 2;
    _user_data.rootdist[4] = 2;
    _user_data.rootdist[5] = 1;
    root_equilibrium_distribution dist(_user_data.rootdist);

    gamma_model model(&lambda, &families, { 0.01, 0.05 }, { 0.1, 0.5 }, NULL);

    vector<double> cat_likelihoods;
    CHECK(model.prune(families[0], dist, cache, &lambda, p_tree.get(), cat_likelihoods));

    CHECK_EQ(2, cat_likelihoods.size());
    CHECK_EQ(doctest::Approx(-23.04433), log(cat_likelihoods[0]));
    CHECK_EQ(doctest::Approx(-16.68005), log(cat_likelihoods[1]));
}

TEST_CASE_FIXTURE(Inference, "gamma_model_prune_returns_false_if_saturated" * doctest::skip(true))
{
    vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    single_lambda lambda(0.9);
    matrix_cache cache(&lambda);

    vector<double> cat_likelihoods;

    gamma_model model(&lambda, &families, { 1.0,1.0 }, { 0.1, 0.5 }, NULL);

    CHECK(!model.prune(families[0], _user_data.prior, cache, &lambda, p_tree.get(), cat_likelihoods));
}

TEST_CASE("Inference: matrix_cache_key_handles_floating_point_imprecision")
{
    set<matrix_cache_key> keys;
    double t = 0.0;
    for (int i = 0; i < 31; i++)
    {
        t += 0.1;
        matrix_cache_key key(boundaries(0,t), 0.3);
        keys.insert(key);
    }
    CHECK_EQ(31, keys.size());

    matrix_cache_key key(boundaries(0, 3.0), 0.3);
    CHECK_EQ(1, keys.count(key));
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
    data.max_family_size = 5;
    single_lambda lambda(0.05);
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
    data.max_family_size = 5;
    single_lambda lambda(0.05);
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
    CHECK_THROWS_WITH_AS(p_tree->get_lambda_index(), "Requested lambda index from branch length tree", runtime_error);

}

TEST_CASE("Clade: get_branch_length_throws_from_lambda_tree")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7", true));
    CHECK_EQ(7, p_tree->get_lambda_index());
    CHECK_THROWS_WITH_AS(p_tree->get_branch_length(), "Requested branch length from lambda tree", runtime_error);

}

TEST_CASE("Clade: lambda_tree_root_index_is_1_if_not_specified")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:2)", true));
    CHECK_EQ(1, p_tree->get_lambda_index());
}

TEST_CASE("Clade: parse_newick_throws_exception_for_invalid_lambdas_in_tree")
{
    CHECK_THROWS_WITH_AS(parse_newick("(A:1,B:0):2", true), "Invalid lambda index set for B", runtime_error);
    CHECK_THROWS_WITH_AS(parse_newick("(A:-1,B:2)", true), "Invalid lambda index set for A", runtime_error);
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

TEST_CASE("Inference: multiple_lambda_returns_correct_values")
{
    ostringstream ost;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    map<string, int> key;
    key["A"] = 5;
    key["B"] = 3;
    multiple_lambda ml(key, { .03, .05, .07, .011, .013, .017 });
    CHECK_EQ(.017, ml.get_value_for_clade(p_tree->find_descendant("A")));
    CHECK_EQ(.011, ml.get_value_for_clade(p_tree->find_descendant("B")));
}

TEST_CASE("Simulation: specified_distribution__select_root_size__returns_exact_selection")
{
    user_data data;
    for (int i = 0; i < 20; ++i)
        data.rootdist[i] = 1;

    root_equilibrium_distribution sd(data.rootdist);
    for (size_t i = 0; i < 20; ++i)
        CHECK_EQ(i, sd.select_root_size(i));
}

TEST_CASE("Simulation: gamma_model_get_simulation_lambda_uses_multiplier_based_on_category_probability")
{
    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    single_lambda lam(0.05);
    gamma_model m(&lam, NULL, gamma_categories, multipliers, NULL);
    vector<double> results(100);
    generate(results.begin(), results.end(), [&m]() {
        unique_ptr<single_lambda> new_lam(dynamic_cast<single_lambda*>(m.get_simulation_lambda()));
        return new_lam->get_single_lambda();
        });

    CHECK_EQ(doctest::Approx(0.057), accumulate(results.begin(), results.end(), 0.0) / 100.0);

}

TEST_CASE("Inference: model_vitals")
{
    mock_model model;
    single_lambda lambda(0.05);
    model.set_lambda(&lambda);
    std::ostringstream ost;
    model.write_vital_statistics(ost, new clade("A", 5), 0.01);
    STRCMP_CONTAINS("Model mockmodel Final Likelihood (-lnL): 0.01", ost.str().c_str());
    STRCMP_CONTAINS("Lambda:            0.05", ost.str().c_str());
    STRCMP_CONTAINS("Maximum possible lambda for this topology: 0.2", ost.str().c_str());
    STRCMP_CONTAINS("No attempts made", ost.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "gene_transcript_reconstrctor__print_increases_decreases_by_family__adds_flag_for_significance")
{
    ostringstream insignificant;
    base_model_reconstruction bmr;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 5;
    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    order.clear();
    bmr.print_increases_decreases_by_family(insignificant, order, { gf }, { 0.03 }, 0.01);
    STRCMP_CONTAINS("myid\t0.03\tn", insignificant.str().c_str());

    ostringstream significant;
    bmr.print_increases_decreases_by_family(insignificant, order, { gf }, { 0.03 }, 0.05);
    STRCMP_CONTAINS("myid\t0.03\ty", insignificant.str().c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "print_increases_decreases_by_family__prints_significance_level_in_header")
{
    ostringstream ost;
    base_model_reconstruction bmr;
    gene_transcript gf;

    bmr.print_increases_decreases_by_family(ost, order, { gf }, { 0.07 }, 0.00001);
    STRCMP_CONTAINS("#FamilyID\tpvalue\tSignificant at 1e-05\n", ost.str().c_str());
}

TEST_CASE("Reconstruction: base_model_print_increases_decreases_by_family")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;
    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    base_model_reconstruction bmr;
    bmr.print_increases_decreases_by_family(empty, order, {}, {}, 0.05);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 5;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    gf.set_expression_value("B", 2);

    ostringstream ost;
    bmr.print_increases_decreases_by_family(ost, order, { gf }, { 0.07 }, 0.05);
    STRCMP_CONTAINS("#FamilyID\tpvalue\tSignificant at 0.05\n", ost.str().c_str());
    STRCMP_CONTAINS("myid\t0.07\tn\n", ost.str().c_str());
}

TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_family")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;
    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    vector<double> multipliers({ .2, .75 });
    vector<gamma_bundle*> bundles; //  ({ &bundle });
    vector<double> em;
    gamma_model_reconstruction gmr(em);
    gmr.print_increases_decreases_by_family(empty, order, {}, {}, 0.05);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    gmr._reconstructions["myid"].reconstruction[p_tree->find_descendant("AB")] = 5;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    gf.set_expression_value("B", 2);

    ostringstream ost;
    gmr.print_increases_decreases_by_family(ost, order, { gf }, { 0.07 }, 0.05);
    STRCMP_CONTAINS("#FamilyID\tpvalue\tSignificant at 0.05", ost.str().c_str());
    STRCMP_CONTAINS("myid\t0.07\tn", ost.str().c_str());
}

TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    ostringstream empty;

    vector<double> multipliers({ .2, .75 });
    vector<gamma_bundle*> bundles; //  ({ &bundle });
    vector<double> em;
    gamma_model_reconstruction gmr(em);

    gmr.print_increases_decreases_by_clade(empty, order, {});
    STRCMP_EQUAL("#Taxon_ID\tIncrease\tDecrease\n", empty.str().c_str());

    gmr._reconstructions["myid"].reconstruction[p_tree->find_descendant("AB")] = 5;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    gf.set_expression_value("B", 2);

    ostringstream ost;
    gmr.print_increases_decreases_by_clade(ost, order, { gf });
    STRCMP_CONTAINS("#Taxon_ID\tIncrease\tDecrease", ost.str().c_str());
    STRCMP_CONTAINS("A<0>\t1\t0", ost.str().c_str());
    STRCMP_CONTAINS("B<1>\t0\t1", ost.str().c_str());
}

TEST_CASE("Reconstruction: base_model_print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    clade invalid;
    cladevector order{ p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    ostringstream empty;

    base_model_reconstruction bmr;

    bmr.print_increases_decreases_by_clade(empty, order, {});
    STRCMP_EQUAL("#Taxon_ID\tIncrease\tDecrease\n", empty.str().c_str());

    //    bmr._reconstructions["myid"][p_tree->find_descendant("A")] = 4;
    //    bmr._reconstructions["myid"][p_tree->find_descendant("B")] = -3;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 5;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    gf.set_expression_value("B", 2);

    ostringstream ost;
    bmr.print_increases_decreases_by_clade(ost, order, { gf });
    STRCMP_CONTAINS("#Taxon_ID\tIncrease\tDecrease", ost.str().c_str());
    STRCMP_CONTAINS("A<0>\t1\t0", ost.str().c_str());
    STRCMP_CONTAINS("B<1>\t0\t1", ost.str().c_str());
}

TEST_CASE("Inference: lambda_epsilon_optimizer")
{
    const double initial_epsilon = 0.01;
    error_model err;
    err.set_probabilities(0, { .0, .99, initial_epsilon });
    err.set_probabilities(1, { initial_epsilon, .98, initial_epsilon });

    mock_model model;

    single_lambda lambda(0.05);
    user_data ud;
    lambda_epsilon_optimizer optimizer(&model, &err, ud, &lambda, 10, 3);
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

TEST_CASE("Inference: lambda_per_family" * doctest::skip(true))
{
    randomizer_engine.seed(10);
    user_data ud;

    ud.max_root_family_size = 10;
    ud.max_family_size = 10;
    ud.gene_families.resize(1);
    ud.prior = root_equilibrium_distribution(size_t(ud.max_root_family_size));

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

TEST_CASE_FIXTURE(Inference, "estimator_compute_pvalues" * doctest::skip(true))
{
    input_parameters params;
    single_lambda s(0.05);
    matrix_cache cache(&s);

    pvalue_parameters p = { _user_data.p_tree,  _user_data.p_lambda, _user_data.max_family_size, _user_data.max_root_family_size, cache };

    auto values = compute_pvalues(p, _user_data.gene_families, 3);
    CHECK_EQ(1, values.size());
    CHECK_EQ(doctest::Approx(0.0), values[0]);
}

TEST_CASE_FIXTURE(Inference, "gamma_lambda_optimizer updates model alpha and lambda")
{
    _user_data.max_root_family_size = 10;
    _user_data.prior = root_equilibrium_distribution(size_t(_user_data.max_root_family_size));

//    gamma_model m(_user_data.p_lambda, _user_data.p_tree, &_user_data.gene_families, 10, _user_data.max_root_family_size, 4, 0.25, NULL);
    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    gamma_model m(_user_data.p_lambda, &_user_data.gene_families, gamma_categories, multipliers, NULL);

    gamma_lambda_optimizer optimizer(_user_data.p_lambda, &m, _user_data, 7, 1);
    vector<double> values{ 0.01, 0.25 };
    optimizer.calculate_score(&values[0]);
    CHECK_EQ(doctest::Approx(0.25), m.get_alpha());
    CHECK_EQ(doctest::Approx(0.01), dynamic_cast<single_lambda *>(m.get_lambda())->get_single_lambda());
}

TEST_CASE("Inference: inference_optimizer_scorer__calculate_score__translates_nan_to_inf")
{
    single_lambda lam(0.05);
    mock_model m;
    m.set_invalid_likelihood();
    double val;
    user_data ud;
    sigma_optimizer_scorer opt(&lam, &m, ud, 0, 0);
    CHECK(std::isinf(opt.calculate_score(&val)));
}

TEST_CASE_FIXTURE(Inference, "poisson_scorer_optimizes_correct_value")
{
    poisson_scorer scorer(_user_data.gene_families);
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

TEST_CASE_FIXTURE(Inference, "poisson_scorer__lnlPoisson")
{
    poisson_scorer scorer(_user_data.gene_families);
    double lambda = 0.05;
    CHECK_EQ(doctest::Approx(3.095732), scorer.lnLPoisson(&lambda));
}

TEST_CASE_FIXTURE(Inference, "poisson_scorer__lnlPoisson_skips_incalculable_family_sizes")
{
    _user_data.gene_families.push_back(gene_transcript("TestFamily2", "", ""));
    _user_data.gene_families[1].set_expression_value("A", 3);
    _user_data.gene_families[1].set_expression_value("B", 175);

    poisson_scorer scorer(_user_data.gene_families);
    double lambda = 0.05;
    CHECK_EQ(doctest::Approx(9.830344), scorer.lnLPoisson(&lambda));
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
    STRCMP_CONTAINS("Best matches are:            0.05,0.03", ost.str().c_str());
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

TEST_CASE("Simulation, specified_distribution__with_rootdist_creates_matching_vector")
{
    std::map<int, int> m;
    m[2] = 3;
    m[4] = 1;
    m[8] = 1;
    root_equilibrium_distribution rd(m);
    CHECK_EQ(rd.select_root_size(0), 2);
    CHECK_EQ(rd.select_root_size(1), 2);
    CHECK_EQ(rd.select_root_size(2), 2);
    CHECK_EQ(rd.select_root_size(3), 4);
    CHECK_EQ(rd.select_root_size(4), 8);
    CHECK_EQ(rd.select_root_size(5), 0);
}

TEST_CASE("Simulation, specified_distribution__pare")
{
    randomizer_engine.seed(10);

    std::map<int, int> m;
    m[2] = 5;
    m[4] = 3;
    m[8] = 3;
    root_equilibrium_distribution rd(m);
    rd.resize(5);
    CHECK_EQ(rd.select_root_size(0), 2);
    CHECK_EQ(rd.select_root_size(1), 2);
    CHECK_EQ(rd.select_root_size(2), 2);
    CHECK_EQ(rd.select_root_size(3), 4);
    CHECK_EQ(rd.select_root_size(4), 8);
    CHECK_EQ(rd.select_root_size(5), 0);
}

TEST_CASE("Simulation, simulate_processes" * doctest::skip(true))
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    mock_model m;
    m.set_lambda(&lam);

    user_data ud;
    ud.p_tree = p_tree.get();
    ud.p_lambda = &lam;
    ud.prior = root_equilibrium_distribution(size_t(100));
    ud.max_family_size = 101;
    ud.max_root_family_size = 101;

    input_parameters ip;
    ip.nsims = 100;
    simulator sim(ud, ip);
    vector<simulated_family> results(1);
    sim.simulate_processes(&m, results);
    CHECK_EQ(100, results.size());
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
    std::vector<lambda*> cache(100);
    CHECK_EQ(0.0, LikelihoodRatioTest::get_likelihood_for_diff_lambdas(gf, p_tree.get(), 0, 0, cache, &opt, 12, 12));
}

TEST_CASE("LikelihoodRatioTest, compute_for_diff_lambdas" * doctest::skip(true))
{
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    user_data data;
    data.p_lambda = &lam;
    data.p_tree = p_tree.get();
    data.gene_families.resize(1);
    data.gene_families[0].set_expression_value("A", 5);
    data.gene_families[0].set_expression_value("B", 9);
    data.max_root_family_size = 12;
    data.max_family_size = 12;
    vector<int> lambda_index(data.gene_families.size(), -1);
    vector<double> pvalues(data.gene_families.size());
    vector<lambda*> lambdas(100);
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

    int res = context.run(); // run

    if (context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests

    return res;
}
