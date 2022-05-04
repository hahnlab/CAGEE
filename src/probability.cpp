#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <numeric>
#include <cmath>
#include <memory>
#include <iterator>

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
#include "sigma.h"
#include "DiffMat.h"
#include "core.h"

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

Eigen::MatrixXd chooseln_cache(100,100);

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

    for (int i = 0; i < chooseln_cache.rows(); ++i)
        for (int j = 0; j < chooseln_cache.cols(); ++j)
            chooseln_cache(i, j) = lgamma2(i + 1) - lgamma2(j + 1) - lgamma2(i - j + 1);
}

double chooseln(double n, double r)
{
  if (r == 0 || (n == 0 && r == 0)) return 0;
  else if (n <= 0 || r <= 0) return log(0);

  if (n >= 0 && r >= 0 && n < chooseln_cache.size() && r < chooseln_cache.size() && (n - int(n) < 0.00000000001) && (r - int(r) < 0.00000000001))
      return chooseln_cache(int(n), int(r));

  return lgamma2(n + 1) - lgamma2(r + 1) - lgamma2(n - r + 1);
}


/* END: Math tools ----------------------- */


int upper_bound_calculator::get_max_bound(const vector<gene_transcript>& transcripts) const
{
    vector<int> bounds(transcripts.size());
    transform(transcripts.begin(), transcripts.end(), bounds.begin(), [this](const gene_transcript& gf) {
        return get(gf);
        });

    return *max_element(bounds.begin(), bounds.end());
}

class upper_bound_calculator_linear_space : public upper_bound_calculator
{
public:
    virtual int get(const gene_transcript& gt) const override
    {
        int val = gt.get_max_expression_value() * MATRIX_SIZE_MULTIPLIER;

        int remainder = val % BOUNDING_STEP_SIZE;
        if (remainder == 0 && val > 0) return val;

        return val + BOUNDING_STEP_SIZE - remainder;
    }
};

class upper_bound_calculator_log_space : public upper_bound_calculator
{
public:
    virtual int get(const gene_transcript& gt) const override
    {
        return max(1.0, ceil(gt.get_max_expression_value() * MATRIX_SIZE_MULTIPLIER));
    }
};

upper_bound_calculator* upper_bound_calculator::create()
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    return new upper_bound_calculator_log_space();
#else
    return new upper_bound_calculator_linear_space();
#endif
}


void print_probabilities(const std::map<const clade*, VectorXd>& probabilities, const clade *node)
{
    auto& probs = probabilities.at(node);
    cout << "node " << node->get_taxon_name() << ":";
    for (int i = 0; i < probs.size(); ++i)
        if (probs[i] != 0)
            cout << i << "=" << probs[i] << ' ';
    cout << endl;
}

void compute_node_probability(const clade* node,
    const gene_transcript& gene_transcript,
    const error_model* p_error_model,
    std::map<const clade*, VectorXd>& probabilities,
    const sigma_squared* p_sigma,
    const matrix_cache& cache,
    int upper_bound)
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
            probabilities[node] = VectorPos_bounds(species_size, DISCRETIZATION_RANGE, boundaries(0, upper_bound));
            //print_probabilities(probabilities, node);
        }
    }
    else  {
        auto& node_probs = probabilities[node];
        node_probs = VectorXd::Constant(DISCRETIZATION_RANGE, 1);

        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            const MatrixXd& m = cache.get_matrix((*it)->get_branch_length(), p_sigma->get_named_value(*it, gene_transcript), upper_bound);

            VectorXd result = m * probabilities[*it];
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


double get_value(const gene_transcript& t, const VectorXd& likelihood, int upper_bound)
{
    Index max_likelihood_index;
    likelihood.maxCoeff(&max_likelihood_index);
    return (float(max_likelihood_index) / likelihood.size()) * upper_bound;
}

inference_pruner::inference_pruner(const matrix_cache& cache,
    const sigma_squared* sigma,
    const error_model* p_error_model,
    const clade* p_tree,
    double sigma_multiplier) :
        _cache(cache),  _p_sigsqd(sigma), _p_error_model(p_error_model), _p_tree(p_tree), _sigma_multiplier(sigma_multiplier)
    {
    }


clademap<VectorXd> inference_pruner::compute_all_probabilities(const gene_transcript& gf, int upper_bound)
{
    unique_ptr<sigma_squared> multiplier(_p_sigsqd->multiply(_sigma_multiplier));
    clademap<VectorXd> probabilities;

    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(_p_tree->reverse_level_begin(), _p_tree->reverse_level_end(), init_func);

    auto compute_func = [&](const clade* c) { compute_node_probability(c, gf, _p_error_model, probabilities, multiplier.get(), _cache, upper_bound); };
    for_each(_p_tree->reverse_level_begin(), _p_tree->reverse_level_end(), compute_func);

    return probabilities;

}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by calling \ref compute_node_probability on all nodes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_pruner::prune(const gene_transcript& gf, int upper_bound)
{
    auto probabilities = compute_all_probabilities(gf, upper_bound);

    return vector<double>(probabilities.at(_p_tree).data(), probabilities.at(_p_tree).data() + probabilities.at(_p_tree).size());
}

clademap<double> inference_pruner::reconstruct(const gene_transcript& gf, int upper_bound)
{
    auto probabilities = compute_all_probabilities(gf, upper_bound);

    clademap<double> reconstruction;
    for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
    {
        reconstruction[*it] = (*it)->is_leaf() ? gf.get_expression_value((*it)->get_taxon_name()) : get_value(gf, probabilities[*it], upper_bound);
    }
    return reconstruction;
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
    family.set_expression_value("A", 18.7);
    family.set_expression_value("B", 17.3);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    sigma_squared lambda(0.03);
    matrix_cache cache;

    std::map<const clade*, VectorXd> probabilities;

    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, NULL, probabilities, &lambda, cache, 60);
    auto& actual = probabilities[A];

    vector<double> expected(DISCRETIZATION_RANGE);

    if (MATRIX_SIZE_MULTIPLIER == 1.5)
    {
        expected[93] = 4.8133125;
        expected[94] = 0.1616875;
    }
    else
    {
        expected[62] = 3.2448055556;
        expected[63] = 0.0718611111;
    }

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], "At index " + to_string(i));
    }

    auto B = p_tree->find_descendant("B");
    compute_node_probability(B, family, NULL, probabilities, &lambda, cache, 60);
    actual = probabilities[B];

    expected = vector<double>(DISCRETIZATION_RANGE);

    if (MATRIX_SIZE_MULTIPLIER == 1.5)
    {
        expected[86] = 4.6391875;
        expected[87] = 0.3358125;
    }
    else
    {
        expected[57] = 2.0618611111;
        expected[58] = 1.2548055556;
    }

    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], "At index " + to_string(i));
    }
}

TEST_CASE("Inference: likelihood_computer_sets_root_nodes_correctly")
{
    ostringstream ost;
    gene_transcript family;
    family.set_expression_value("A", 14.7);
    family.set_expression_value("B", 22.3);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    double prob = 0.005;
    VectorXd equal_probs = VectorXd::Constant(DISCRETIZATION_RANGE, prob);
    MatrixXd doubler = MatrixXd::Identity(DISCRETIZATION_RANGE, DISCRETIZATION_RANGE) * 2;
    sigma_squared lambda(0.03);
    matrix_cache cache;
    cache.set_matrix(1, 0.03, 40, doubler);
    cache.set_matrix(3, 0.03, 40, doubler);
    std::map<const clade*, VectorXd> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(DISCRETIZATION_RANGE); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto AB = p_tree->find_descendant("AB");
    probabilities[p_tree->find_descendant("A")] = equal_probs;
    probabilities[p_tree->find_descendant("B")] = equal_probs;

    compute_node_probability(AB, family, NULL, probabilities, &lambda, cache, 40);

    auto& actual = probabilities[AB];
    double expected = (2 * prob) * (2 * prob);

    CHECK_EQ(DISCRETIZATION_RANGE, actual.size());
    for (Eigen::Index i = 0; i < actual.size(); ++i)
    {
        CHECK_MESSAGE(doctest::Approx(expected) == actual[i], "At index " + to_string(i));
    }
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_from_error_model_if_provided")
{
    ostringstream ost;

    gene_transcript family;
    family.set_expression_value("A", 17.4);
    family.set_expression_value("B", 22.9);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    sigma_squared lambda(0.03);
    matrix_cache cache;

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
    compute_node_probability(A, family, &model, probabilities, &lambda, cache, 40);
    VectorXd& actual = probabilities[A];

    vector<double> expected(DISCRETIZATION_RANGE);
    expected[16] = 0.2;
    expected[17] = 0.6;
    expected[18] = 0.2;

    REQUIRE(expected.size() == actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], "At index " + to_string(i));
    }
}

TEST_CASE("find_best_pvalue")
{
    gene_transcript gt("TestFamily1", "", "");
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
    gene_transcript gt("TestFamily1", "", "");
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
    gene_transcript gt("TestFamily1", "", "");
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

TEST_CASE("inference_pruner: check stats of returned probabilities")
{
    ostringstream ost;

    gene_transcript fam;
    fam.set_expression_value("A", 3);
    fam.set_expression_value("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    sigma_squared ss(10.0);
    matrix_cache cache;
    cache.precalculate_matrices(ss.get_values(), set<double>{1, 3, 7}, 20);
    inference_pruner pruner(cache, &ss, nullptr, p_tree.get(), 1.0);

    auto actual = pruner.prune(fam, 20);

    size_t sz = actual.size();
    CHECK_EQ(sz, DISCRETIZATION_RANGE);
    auto mean = accumulate(actual.begin(), actual.end(), 0.0) / sz;
    auto variance = accumulate(actual.begin(), actual.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));
        });
    double max = *max_element(actual.begin(), actual.end());

    CHECK_EQ(doctest::Approx(-5.542), log(mean));
    CHECK_EQ(doctest::Approx(-10.62318), log(variance));
    CHECK_EQ(doctest::Approx(-4.37886), log(max));

}

/// Bounds are next largest integer if in log space.
TEST_CASE("Bounds returns next integer of the largest value times MATRIX_SIZE_MULTIPLIER")
{
    upper_bound_calculator_log_space calc;
    gene_transcript gt;
    gt.set_expression_value("A", .3);
    gt.set_expression_value("B", .8);
    CHECK_EQ(3, calc.get(gt));

    gt.set_expression_value("A", 1.1);
    CHECK_EQ(4, calc.get(gt));

    gt.set_expression_value("B", 1.8);
    CHECK_EQ(6, calc.get(gt));
}

TEST_CASE("Bounds never returns less than 1")
{
    upper_bound_calculator_log_space calc;
    gene_transcript gt;
    gt.set_expression_value("A", 0.00005);
    gt.set_expression_value("B", 0.00004);
    CHECK_EQ(1, calc.get(gt));
}

TEST_CASE("Bounds never returns less than 1 even if all values are very small")
{
    upper_bound_calculator_log_space calc;
    gene_transcript gt;
    gt.set_expression_value("A", 0.0000000002);
    gt.set_expression_value("B", 0.0000000005);
    CHECK_EQ(1, calc.get(gt));
}

TEST_CASE("Bounds never returns less than 1 even if all values are zero")
{
    gene_transcript gt;
    gt.set_expression_value("A", 0);
    gt.set_expression_value("B", 0);
    upper_bound_calculator_log_space calc;
    CHECK_EQ(1, calc.get(gt));
}

TEST_CASE("Bounds never returns less than 20")
{
    upper_bound_calculator_linear_space calc;
    gene_transcript gt;
    gt.set_expression_value("A", 5);
    gt.set_expression_value("B", 4);
    CHECK_EQ(20, calc.get(gt));
}

TEST_CASE("Bounds never returns less than 20 even if all values are less than .3")
{
    upper_bound_calculator_linear_space calc;
    gene_transcript gt;
    gt.set_expression_value("A", 0.254007);
    gt.set_expression_value("B", 0.1);
    CHECK_EQ(20, calc.get(gt));
}

TEST_CASE("get_value")
{
    gene_transcript gt;
    gt.set_expression_value("A", 0.254007);
    gt.set_expression_value("B", 0.1);

    VectorXd v = VectorXd::Zero(200);
    v(4) = 8;
    v(6) = 12;
    v(112) = 22;
    v(158) = 9;

    CHECK_EQ(56, get_value(gt, v, 100));
}