#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <numeric>
#include <cmath>
#include <memory>

#include "doctest.h"
#include "easylogging++.h"

#ifdef HAVE_VECTOR_EXP
#include "mkl.h"
#endif

#include "clade.h"
#include "probability.h"
#include "matrix_cache.h"
#include "gene_transcript.h"
#include "error_model.h"
#include "simulator.h"
#include "lambda.h"
#include "DiffMat.h"

using namespace std;
using namespace Eigen;

extern std::mt19937 randomizer_engine;

/* Useful links
1) http://www.rskey.org/gamma.htm # explanation for lgamma
2) http://www.physics.unlv.edu/~pang/cp_c.html # c code
*/

/* Necessary for old C implementation of what now is lgamma
#define M_SQRT_2PI		2.5066282746310002416123552393401042  // sqrt(2pi)
*/

/* Necessary for old C implementation of what now is lgamma 
double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
 24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
 -5.395239384953e-6 };
*/

/* Old C implementation of what now is lgamma */
/*
double gammaln(double a)
{
    int n;
    double p = __Qs[0];
    double a_add_5p5 = a + 5.5;
    for (n = 1; n <= 6; n++) p += __Qs[n] / (a + n);
    return (a + 0.5)*log(a_add_5p5) - (a_add_5p5)+log(M_SQRT_2PI*p / a);
}
*/


#define GAMMA_CACHE_SIZE 1024

vector<double> lgamma_cache;

matrix chooseln_cache(100);

inline double lgamma2(double n)
{
    if (n >= 0 && n < GAMMA_CACHE_SIZE && (n - int(n) < 0.00000000001))
        return lgamma_cache.at(int(n));

    return lgamma(n);
}

void init_lgamma_cache()
{
    lgamma_cache.resize(GAMMA_CACHE_SIZE);
    for (int i = 0; i < GAMMA_CACHE_SIZE; ++i)
    {
        lgamma_cache[i] = lgamma(i);
    }

    for (int i = 0; i < chooseln_cache.size(); ++i)
        for (int j = 0; j < chooseln_cache.size(); ++j)
            chooseln_cache.set(i, j, lgamma2(i + 1) - lgamma2(j + 1) - lgamma2(i - j + 1));
}

double chooseln(double n, double r)
{
  if (r == 0 || (n == 0 && r == 0)) return 0;
  else if (n <= 0 || r <= 0) return log(0);

  if (n >= 0 && r >= 0 && n < chooseln_cache.size() && r < chooseln_cache.size() && (n - int(n) < 0.00000000001) && (r - int(r) < 0.00000000001))
      return chooseln_cache.get(n, r);

  return lgamma2(n + 1) - lgamma2(r + 1) - lgamma2(n - r + 1);
}

/* END: Math tools ----------------------- */

//! Calculates the probabilities of a given node for a given family size.
//! The probability of a leaf node is given from the family size at that node (1 for the family
//! size, and 0 for all other sizes). Internal nodes are calculated based on the probabilities
//! of the descendants.
//! Results are stored in the probabilities vector, 0 - max possible size
void compute_node_probability(const clade* node,
    const gene_transcript& gene_transcript,
    const error_model* p_error_model,
    std::map<const clade*, VectorXd>& probabilities,
    const lambda* p_sigma,
    const DiffMat& diff_mat)
{
    if (node->is_leaf()) {
        double species_size = gene_transcript.get_expression_value(node->get_taxon_name());

        if (p_error_model != NULL)
        {
            auto error_model_probabilities = p_error_model->get_probs(species_size);
            int offset = species_size - ((p_error_model->n_deviations() - 1) / 2);
            for (size_t i = 0; i < error_model_probabilities.size(); ++i)
            {
                if (offset + int(i) < 0)
                    continue;

                probabilities[node][offset + i] = error_model_probabilities[i];
            }
        }
        else
        {
            // cout << "Leaf node " << node->get_taxon_name() << " has " << _probabilities[node].size() << " probabilities" << endl;
            probabilities[node] = VectorPos_bounds(species_size, DISCRETIZATION_RANGE, std::pair<double, double>(0, probabilities[node].size()));
        }
    }
    else  {
        auto& node_probs = probabilities[node];
        node_probs = VectorXd::Constant(DISCRETIZATION_RANGE, 1);

        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            double sigma = p_sigma->get_value_for_clade(*it);
            MatrixXd m = ConvProp_bounds((*it)->get_branch_length(), sigma * sigma / 2, diff_mat, pair<double, double>(0.0, gene_transcript.get_max_expression_value() * 1.5));

            auto result = m * probabilities[*it];
            for (VectorXd::Index i = 0; i < node_probs.size(); i++) {
                node_probs[i] *= result[i];
            }
        }
    }
}


/* END: Likelihood computation ---------------------- */

std::vector<int> uniform_dist(int n_draws, int min, int max) {
    
    std::random_device rd; // seed
    std::default_random_engine generator(rd()); // seeding generator
    std::uniform_int_distribution<int> distribution(min, max); // initializing uniform generator
    std::vector<int> uniform_vec(n_draws); // for storing results
    
    for (int i = 0; i < n_draws; ++i) {
        int number = distribution(generator); // drawing from uniform generator by plugging in random number
        //cout << "Number is: " << number << endl;
        uniform_vec[i] = number;
    }
        
    return uniform_vec;
}

vector<double> compute_family_probabilities(pvalue_parameters p, vector<simulated_family>& sizes, int root_family_size)
{
    vector<double> result(sizes.size());

    // Allocate space to calculate all of the families simultaneously
    vector<clademap<VectorXd>> pruners(sizes.size());
    for_each(pruners.begin(), pruners.end(), [&p](clademap<VectorXd>& pruner) {
        // vector of lk's at tips must go from 0 -> _max_possible_family_size, so we must add 1
        for_each(p.p_tree->reverse_level_begin(), p.p_tree->reverse_level_end(), [&p, &pruner](const clade* node) 
            { 
                pruner[node] = VectorXd::Zero(DISCRETIZATION_RANGE);
            });
        });

    // get a gene family for each clademap
    vector<gene_transcript> families(sizes.size());
    transform(sizes.begin(), sizes.end(), families.begin(), [](const simulated_family& s) {
        gene_transcript f;
        for (auto& it : s.values) {
            if (it.first->is_leaf())
            {
                f.set_expression_value(it.first->get_taxon_name(), it.second);
            }
        }
        return f;
    });

    // do math
#pragma omp parallel for
    for (size_t i = 0; i < result.size(); ++i)
    {
        for (auto it = p.p_tree->reverse_level_begin(); it != p.p_tree->reverse_level_end(); ++it)
            compute_node_probability(*it, families[i], NULL, pruners[i], p.p_lambda, p.diff_mat);
        result[i] = *std::max_element(pruners[i].at(p.p_tree).data(), pruners[i].at(p.p_tree).data() + pruners[i].at(p.p_tree).size());
    }
    return result;
}

/*! Create a sorted vector of probabilities by generating random trees 
    \param p_tree The structure of the tree to generate
    \param number_of_simulations The number of random probabilities to return
    \param root_family_size The count of the family at the root. All other family sizes will be generated randomly based on this
    \param max_family_size The maximum possible family size (Used to cut off probability calculations at a reasonable value)
    \param root_family_size The maximum possible family size at the root
    \param lambda The rate of change of the family

    \returns a sorted vector of probabilities of the requested number of randomly generated trees
*/
std::vector<double> get_random_probabilities(pvalue_parameters p, int number_of_simulations, int root_family_size)
{
    vector<simulated_family> families(number_of_simulations);

    generate(families.begin(), families.end(), [p, root_family_size]() { return create_simulated_family(p.p_tree, p.p_lambda, root_family_size, p.diff_mat); });

    auto result = compute_family_probabilities(p, families, root_family_size);

    sort(result.begin(), result.end());

    return result;
}

size_t adjust_for_error_model(size_t c, const error_model *p_error_model)
{
    if (p_error_model == nullptr)
        return c;

    if (c >= p_error_model->get_max_family_size())
    {
        throw runtime_error("Trying to simulate leaf family size that was not included in error model");
    }
    auto probs = p_error_model->get_probs(c);

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rnd = distribution(randomizer_engine);
    if (rnd < probs[0])
    {
        c--;
    }
    else if (rnd >(1 - probs[2]))
    {
        c++;
    }

    return c;
}

double pvalue(double v, const vector<double>& conddist)
{
    int idx = conddist.size() - 1;

    auto bound = std::upper_bound(conddist.begin(), conddist.end(), v);
    if (bound != conddist.end())
    {
        idx = bound - conddist.begin();
    }
    return  idx / (double)conddist.size();
}

double find_best_pvalue(const gene_transcript& fam, const VectorXd& root_probabilities, const std::vector<std::vector<double> >& conditional_distribution)
{    
    vector<double> pvalues(root_probabilities.size());
    int max_size_to_check = rint(fam.get_max_expression_value() * 1.25);
    for (int j = 0; j < max_size_to_check; ++j)
    {
        pvalues[j] = pvalue(root_probabilities[j], conditional_distribution[j]);
    }
    auto idx = std::max_element(pvalues.begin(), pvalues.end());
    LOG(TRACE) << "PValue for " << fam.id() << " : " << *idx << " found at family size " << idx - pvalues.begin();

    return *idx;
}

//! Compute pvalues for each family based on the given lambda
vector<double> compute_pvalues(pvalue_parameters p, const std::vector<gene_transcript>& families, int number_of_simulations)
{
    LOG(INFO) << "Computing pvalues...";
    LOG(WARNING) << "PValues are not calculating accurate randomizations yet";

    std::vector<std::vector<double> > conditional_distribution(p.max_root_family_size);
    for (int i = 0; i < p.max_root_family_size; ++i)
    {
        conditional_distribution[i] = get_random_probabilities(p, number_of_simulations, i+1);
    }
    VLOG(1) << "Conditional distributions calculated";

//    auto observed_max_likelihoods = compute_family_probabilities(p, families);
    vector<double> observed_max_likelihoods(families.size());

    // Allocate space to calculate all of the families simultaneously
    vector<clademap<VectorXd>> pruners(families.size());
    for (auto& pruner : pruners)
    {
        // vector of lk's at tips must go from 0 -> _max_possible_family_size, so we must add 1
        auto fn = [&](const clade* node) { pruner[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
        for_each(p.p_tree->reverse_level_begin(), p.p_tree->reverse_level_end(), fn);
    }

#pragma omp parallel for
    for (size_t i = 0; i < families.size(); ++i)
    {
        for (auto it = p.p_tree->reverse_level_begin(); it != p.p_tree->reverse_level_end(); ++it)
            compute_node_probability(*it, families[i], NULL, pruners[i], p.p_lambda, p.diff_mat);
    }

    vector<double> result(families.size());
    for (size_t i = 0; i < families.size(); ++i)
    {
        result[i] = find_best_pvalue(families[i], pruners[i].at(p.p_tree), conditional_distribution);
    }


    LOG(INFO) << "done!\n";

    return result;
}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by calling \ref compute_node_probability on all nodes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_prune(const gene_transcript& gf, const DiffMat& diff_mat, const lambda* p_lambda, const error_model* p_error_model, const clade* p_tree, double lambda_multiplier, int max_root_family_size, int max_family_size)
{
    unique_ptr<lambda> multiplier(p_lambda->multiply(lambda_multiplier));
    clademap<VectorXd> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto compute_func = [&](const clade* c) { compute_node_probability(c, gf, p_error_model, probabilities, multiplier.get(), diff_mat); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), compute_func);

    return vector<double>(probabilities.at(p_tree).data(), probabilities.at(p_tree).data() + probabilities.at(p_tree).size()); // likelihood of the whole tree = multiplication of likelihood of all nodes
}

TEST_CASE("VectorPos_bounds")
{
    auto actual = VectorPos_bounds(7.0, 20, pair<double, double>(0, 10));
    vector<double> expected{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.33, 0.57, 0, 0, 0, 0, 0 };
    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(expected[i]), actual[i]);
    }
}

TEST_CASE("VectorPos_bounds at right edge")
{
    auto actual = VectorPos_bounds(10.0, 20, pair<double, double>(0, 10));
    vector<double> expected{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.9 };
    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(expected[i]), actual[i]);
    }
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_correctly")
{
    ostringstream ost;
    gene_transcript family;
    family.set_expression_value("A", 3);
    family.set_expression_value("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);

    std::map<const clade*, VectorXd> probabilities;

    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, NULL, probabilities, &lambda, DiffMat::instance());
    auto& actual = probabilities[A];

    vector<double> expected(DISCRETIZATION_RANGE);
    expected[2] = 0.014925;
    expected[3] = 0.980075;

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], "At index " + to_string(i));
    }

    auto B = p_tree->find_descendant("B");
    compute_node_probability(B, family, NULL, probabilities, &lambda, DiffMat::instance());
    actual = probabilities[B];

    expected = vector<double>(DISCRETIZATION_RANGE);
    expected[5] = 0.02985;
    expected[6] = 0.96515;

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], "At index " + to_string(i));
    }
}

TEST_CASE("Inference: likelihood_computer_sets_root_nodes_correctly" * doctest::skip(true))
{
    ostringstream ost;
    gene_transcript family;
    family.set_expression_value("A", 3);
    family.set_expression_value("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);

    std::map<const clade*, VectorXd> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto AB = p_tree->find_descendant("AB");
    compute_node_probability(p_tree->find_descendant("A"), family, NULL, probabilities, &lambda, DiffMat::instance());
    compute_node_probability(p_tree->find_descendant("B"), family, NULL, probabilities, &lambda, DiffMat::instance());
    compute_node_probability(AB, family, NULL, probabilities, &lambda, DiffMat::instance());

    auto& actual = probabilities[AB];

    vector<double> log_expected{ -19.7743, -11.6688, -5.85672, -5.66748, -6.61256, -8.59725, -12.2301, -16.4424, -20.9882, -25.7574,
        -30.6888, -35.7439, -40.8971, -46.1299, -51.4289, -56.7837, -62.1863, -67.6304, -73.1106, -78.6228
    };

    CHECK_EQ(log_expected.size(), actual.size());
    for (size_t i = 0; i < log_expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(log_expected[i]), log(actual[i]));
    }
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_from_error_model_if_provided" * doctest::skip(true))
{
    ostringstream ost;

    gene_transcript family;
    family.set_expression_value("A", 3);
    family.set_expression_value("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    single_lambda lambda(0.03);

    string input = "maxcnt: 20\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "20 0.2 0.6 0.2\n";
    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    std::map<const clade*, VectorXd> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, &model, probabilities, &lambda, DiffMat::instance());
    auto& actual = probabilities[A];

    vector<double> expected{ 0, 0, 0.2, 0.6, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        //cout << actual[i] << endl;
        CHECK_EQ(expected[i], actual[i]);
    }
}

TEST_CASE("find_best_pvalue")
{
    gene_transcript gt;
    gt.set_id("TestFamily1");
    gt.set_expression_value("A", 1);
    gt.set_expression_value("B", 2);

    std::vector<std::vector<double> > conditional_distribution(3);
    conditional_distribution[0] = vector<double>({ .1, .2, .3, .4 });
    conditional_distribution[1] = vector<double>({ .1, .2, .3, .4 });
    conditional_distribution[2] = vector<double>({ .1, .2, .3, .4 });
    Vector3d root_probabilities;
    root_probabilities[0] = .15;
    auto value = find_best_pvalue(gt, root_probabilities, conditional_distribution);
    CHECK_EQ(doctest::Approx(0.25), value);
}

TEST_CASE("find_best_pvalue_skips_values_outside_of_range")
{
    gene_transcript gt;
    gt.set_id("TestFamily1");
    gt.set_expression_value("A", 1);
    gt.set_expression_value("B", 2);

    std::vector<std::vector<double> > conditional_distribution(3);
    conditional_distribution[0] = vector<double>({ .1, .2, .3, .4 });
    conditional_distribution[1] = vector<double>({ .1, .2, .3, .4 });
    conditional_distribution[2] = vector<double>({ .1, .2, .3, .4 });
    Vector3d root_probabilities;
    root_probabilities[0] = .15;
    root_probabilities[2] = .35;
    auto value = find_best_pvalue(gt, root_probabilities, conditional_distribution);
    CHECK_EQ(doctest::Approx(0.25), value);
}


TEST_CASE("find_best_pvalue_selects_largest_value_in_range")
{
    gene_transcript gt;
    gt.set_id("TestFamily1");
    gt.set_expression_value("A", 1);
    gt.set_expression_value("B", 2);

    std::vector<std::vector<double> > conditional_distribution(3);
    conditional_distribution[0] = vector<double>({ .1, .2, .3, .4 });
    conditional_distribution[1] = vector<double>({ .1, .2, .3, .4 });
    conditional_distribution[2] = vector<double>({ .1, .2, .3, .4 });
    Vector3d root_probabilities;
    root_probabilities[0] = .15;
    root_probabilities[1] = .25;
    auto value = find_best_pvalue(gt, root_probabilities, conditional_distribution);
    CHECK_EQ(doctest::Approx(0.5), value);
}

TEST_CASE("compute_family_probabilities")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

    single_lambda lambda(0.03);
    pvalue_parameters p = { p_tree.get(),  &lambda, 20, 15, DiffMat::instance() };

    vector<simulated_family> v(1);
    v[0].values[p_tree.get()] = 5;

    gene_transcript fam;
    fam.set_id("Family5");
    fam.set_expression_value("A", 11);
    fam.set_expression_value("B", 2);
    fam.set_expression_value("C", 5);
    fam.set_expression_value("D", 6);

    // note we do not use an error model for creating family sizes. See architecture decision #6
    p.p_tree->apply_prefix_order([&v, &fam](const clade* c)
        {
            if (c->is_leaf())
                v[0].values[c] = fam.get_expression_value(c->get_taxon_name());
            else
                v[0].values[c] = 5;
        });
    auto result = compute_family_probabilities(p, v, 5);

    CHECK_EQ(1, result.size());
    CHECK_EQ(doctest::Approx(0.0000000001), result[0]);
}

