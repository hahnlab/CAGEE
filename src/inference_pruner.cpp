#include <algorithm>

#include "doctest.h"
#include "easylogging++.h"

#ifdef HAVE_VECTOR_EXP
#include "mkl.h"
#endif

#include "inference_pruner.h"

#include "clade.h"
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
            VectorPos_bounds(species_size, boundaries(0, upper_bound), probabilities[node]);
            //print_probabilities(probabilities, node);
        }
    }
    else  {
        auto& node_probs = probabilities[node];
        node_probs.fill(1);

        for (auto it = node->descendant_begin(); it != node->descendant_end(); ++it) {
            const MatrixXd& m = cache.get_matrix((*it)->get_branch_length(), p_sigma->get_named_value(*it, gene_transcript), upper_bound);

            VectorXd result = m * probabilities[*it];
            for (VectorXd::Index i = 0; i < node_probs.size(); i++) {
                node_probs[i] *= result[i];
            }
        }
    }
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

node_reconstruction get_value(const VectorXd& likelihood, int upper_bound)
{
    auto midpoint = [&likelihood, upper_bound](int idx) {
        double bin_start = (double(idx) / likelihood.size()) * upper_bound;
        double bin_width = double(upper_bound) / double(likelihood.size());
        return bin_start + bin_width / 2.0;
    };
    Index max_likelihood_index;
    likelihood.maxCoeff(&max_likelihood_index);
    VectorXd cumsum(likelihood.size());
    partial_sum(likelihood.begin(), likelihood.end(), cumsum.begin());

    VectorXd nml = cumsum / cumsum.tail(1)(0);
    auto low_bound = std::upper_bound(nml.begin(), nml.end(), 0.025);
    auto high_bound = std::upper_bound(nml.begin(), nml.end(), 0.975);

    node_reconstruction nr;
    nr.most_likely_value = midpoint(max_likelihood_index);
    nr.credible_interval.first = midpoint(low_bound - nml.begin());
    nr.credible_interval.second = midpoint(high_bound - nml.begin());
    return nr;
}

inference_pruner::inference_pruner(const matrix_cache& cache,
    const sigma_squared* sigma,
    const error_model* p_error_model,
    const clade* p_tree,
    double sigma_multiplier) :
        _cache(cache),  _p_sigsqd(sigma), _p_error_model(p_error_model), _p_tree(p_tree), _sigma_multiplier(sigma_multiplier)
{
    auto init_func = [&](const clade* node) { _probabilities[node] = _cache.create_vector(); };
    for_each(_p_tree->reverse_level_begin(), _p_tree->reverse_level_end(), init_func);
}


void inference_pruner::compute_all_probabilities(const gene_transcript& gf, int upper_bound)
{
    unique_ptr<sigma_squared> multiplier(_p_sigsqd->multiply(_sigma_multiplier));

    auto compute_func = [gf, this, upper_bound, &multiplier](const clade* c) { compute_node_probability(c, gf, _p_error_model, _probabilities, multiplier.get(), _cache, upper_bound); };
    for_each(_p_tree->reverse_level_begin(), _p_tree->reverse_level_end(), compute_func);

}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by calling \ref compute_node_probability on all nodes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_pruner::prune(const gene_transcript& gf, int upper_bound)
{
    compute_all_probabilities(gf, upper_bound);

    return vector<double>(_probabilities[_p_tree].begin(), _probabilities[_p_tree].end());
}

clademap<node_reconstruction> inference_pruner::reconstruct(const gene_transcript& gf, int upper_bound)
{
    compute_all_probabilities(gf, upper_bound);

    clademap<node_reconstruction> reconstruction;
    for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
    {
        if ((*it)->is_leaf())
        {
            node_reconstruction nr;
            nr.most_likely_value = gf.get_expression_value((*it)->get_taxon_name());
            reconstruction[*it] = nr;
        }
        else
        {
            reconstruction[*it] = get_value(_probabilities[*it], upper_bound);
        }
    }
    return reconstruction;
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_correctly")
{
    int Npts = 200;
    ostringstream ost;
    gene_transcript family;
    family.set_expression_value("A", 18.7);
    family.set_expression_value("B", 17.3);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    sigma_squared lambda(0.03);
    matrix_cache cache;

    std::map<const clade*, VectorXd> probabilities;

    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(Npts); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, NULL, probabilities, &lambda, cache, 60);
    auto& actual = probabilities[A];

    vector<double> expected(Npts);

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
        string x = string("At index ") + to_string(i);
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], x);
    }

    auto B = p_tree->find_descendant("B");
    compute_node_probability(B, family, NULL, probabilities, &lambda, cache, 60);
    actual = probabilities[B];

    expected = vector<double>(Npts);

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
        string x = string("At index ") + to_string(i);
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], x);
    }
}

TEST_CASE("Inference: likelihood_computer_sets_root_nodes_correctly")
{
    int Npts = 200;
    ostringstream ost;
    gene_transcript family;
    family.set_expression_value("A", 14.7);
    family.set_expression_value("B", 22.3);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    double prob = 0.005;
    VectorXd equal_probs = VectorXd::Constant(Npts, prob);
    MatrixXd doubler = MatrixXd::Identity(Npts, Npts) * 2;
    sigma_squared lambda(0.03);
    matrix_cache cache;
    cache.set_matrix(1, 0.03, 40, doubler);
    cache.set_matrix(3, 0.03, 40, doubler);
    std::map<const clade*, VectorXd> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(Npts); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto AB = p_tree->find_descendant("AB");
    probabilities[p_tree->find_descendant("A")] = equal_probs;
    probabilities[p_tree->find_descendant("B")] = equal_probs;

    compute_node_probability(AB, family, NULL, probabilities, &lambda, cache, 40);

    auto& actual = probabilities[AB];
    double expected = (2 * prob) * (2 * prob);

    CHECK_EQ(Npts, actual.size());
    for (Eigen::Index i = 0; i < actual.size(); ++i)
    {
        string x = string("At index ") + to_string(i);
        CHECK_MESSAGE(doctest::Approx(expected) == actual[i], x);
    }
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_from_error_model_if_provided")
{
    int Npts = 200;
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
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(Npts); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, family, &model, probabilities, &lambda, cache, 40);
    VectorXd& actual = probabilities[A];

    vector<double> expected(Npts);
    expected[16] = 0.2;
    expected[17] = 0.6;
    expected[18] = 0.2;

    REQUIRE(expected.size() == actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        string x = string("At index ") + to_string(i);
        CHECK_MESSAGE(doctest::Approx(expected[i]) == actual[i], x);
    }
}

TEST_CASE("inference_pruner: check stats of returned probabilities")
{
    int Npts = 200;
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
    CHECK_EQ(sz, Npts);
    auto mean = accumulate(actual.begin(), actual.end(), 0.0) / sz;
    auto variance = accumulate(actual.begin(), actual.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));
        });
    double max = *max_element(actual.begin(), actual.end());

    CHECK_EQ(doctest::Approx(-5.542), log(mean));
    CHECK_EQ(doctest::Approx(-10.62318), log(variance));
    CHECK_EQ(doctest::Approx(-4.37886), log(max));

}

TEST_CASE("get_value most likely value")
{
    VectorXd v = VectorXd::Zero(200);
    v(4) = 8;
    v(6) = 12;
    v(112) = 22;
    v(158) = 9;

    CHECK_EQ(doctest::Approx(56.25), get_value(v, 100).most_likely_value);
}

TEST_CASE("get_value credible interval")
{
    std::normal_distribution<float> dist(10, 3);

    VectorXd v(20);
    generate(v.begin(), v.end(), [&dist]() { return dist(randomizer_engine);  });
    auto i = get_value(v, 20);
    CHECK_EQ(doctest::Approx(0.5), i.credible_interval.first);
    CHECK_EQ(doctest::Approx(19.5), i.credible_interval.second);

}