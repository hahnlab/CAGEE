#include <algorithm>

#include "doctest.h"
#include "easylogging++.h"

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
    int upper_bound)    // upper_bound is the maximum value we expect ever to hit in our expression values. It is NOT the size of our diffusion matrix!
{
    if (node->is_leaf()) {
        double expression_value = gene_transcript.get_expression_value(node->get_taxon_name());

        if (p_error_model != NULL)
        {          
                    int Npts = probabilities[node].size();
                    double nx = (Npts - 1) * (expression_value -0) / double(upper_bound );
                    int ix = floor(nx);
                    size_t l =ix;
                    double ux = nx - ix;
                    double scale =(Npts - 1) / upper_bound; 
                    std::vector<double> v= p_error_model->get_probs(l,ux,scale);
                    Eigen::VectorXd v3 = Eigen::VectorXd::Map(v.data(), v.size());
                    probabilities[node]= v3;
                    
        }
        else
        {
            // cout << "Leaf node " << node->get_taxon_name() << " has " << _probabilities[node].size() << " probabilities" << endl;
            VectorPos_bounds(expression_value, boundaries(0, upper_bound), probabilities[node]);
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

inline void CHECK_VECTORS_EQ(const std::vector<double>& expected, const VectorXd& actual)
{
    CHECK_EQ(expected.size(), actual.size());
    CHECK(mismatch(expected.begin(), expected.end(), actual.begin(), [](double a, double b) { return doctest::Approx(a) == b;  }).first == expected.end());
}

TEST_CASE("compute_node_probability sets leaf probabilities correctly")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    gene_transcript gt;
    gt.set_expression_value("A", 7.0);

    clademap<VectorXd> probabilities;
    matrix_cache cache;
    sigma_squared ss(0.03);
    int upper_bound = 10;

    auto A = p_tree->find_descendant("A");
    auto& actual = probabilities[A];
    actual = VectorXd::Zero(20);

    compute_node_probability(A, gt, nullptr, probabilities, &ss, cache, upper_bound);
    vector<double> expected{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.33, 0.57, 0, 0, 0, 0, 0 };
    CHECK_VECTORS_EQ(expected, actual);
}

TEST_CASE("compute_node_probability sets leaf probabilities correctly with error model")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    gene_transcript gt;
    gt.set_expression_value("A", 7.0);

    clademap<VectorXd> probabilities;
    matrix_cache cache;
    sigma_squared ss(0.03);
    int upper_bound = 10;

    auto A = p_tree->find_descendant("A");
    auto& actual = probabilities[A];

    int Npts = 20;
    actual = VectorXd::Zero(Npts);

    error_model err;
    double nx = double(upper_bound) / Npts;
    for (auto i = nx; i < upper_bound; i = i + nx)
    {
        double l = i;
        err.set_probabilities(0.5, 0.4, l, upper_bound);
    }
    compute_node_probability(A, gt, &err, probabilities, &ss, cache, upper_bound);
    vector<double> expected{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.33, 0.57, 0, 0, 0, 0, 0 };
    CHECK_VECTORS_EQ(expected, actual);
}

TEST_CASE("Inference: likelihood_computer_sets_leaf_nodes_correctly")
{
    int Npts = 200;
    ostringstream ost;
    gene_transcript gt;
    gt.set_expression_value("A", 18.7);
    gt.set_expression_value("B", 17.3);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    sigma_squared ss(0.03);
    matrix_cache cache;

    std::map<const clade*, VectorXd> probabilities;

    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(Npts); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto A = p_tree->find_descendant("A");
    compute_node_probability(A, gt, NULL, probabilities, &ss, cache, 60);
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
    compute_node_probability(B, gt, NULL, probabilities, &ss, cache, 60);
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
    gene_transcript gt;
    gt.set_expression_value("A", 14.7);
    gt.set_expression_value("B", 22.3);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    double prob = 0.005;
    VectorXd equal_probs = VectorXd::Constant(Npts, prob);
    MatrixXd doubler = MatrixXd::Identity(Npts, Npts) * 2;
    sigma_squared ss(0.03);
    matrix_cache cache;
    cache.set_matrix(1, 0.03, 40, doubler);
    cache.set_matrix(3, 0.03, 40, doubler);
    std::map<const clade*, VectorXd> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node] = VectorXd::Zero(Npts); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto AB = p_tree->find_descendant("AB");
    probabilities[p_tree->find_descendant("A")] = equal_probs;
    probabilities[p_tree->find_descendant("B")] = equal_probs;

    compute_node_probability(AB, gt, NULL, probabilities, &ss, cache, 40);

    auto& actual = probabilities[AB];
    double expected = (2 * prob) * (2 * prob);

    CHECK_EQ(Npts, actual.size());
    for (Eigen::Index i = 0; i < actual.size(); ++i)
    {
        string x = string("At index ") + to_string(i);
        CHECK_MESSAGE(doctest::Approx(expected) == actual[i], x);
    }
}

TEST_CASE("inference_pruner: check stats of returned probabilities")
{
    int Npts = 200;
    ostringstream ost;

    gene_transcript rt;
    rt.set_expression_value("A", 3);
    rt.set_expression_value("B", 6);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    sigma_squared ss(10.0);
    matrix_cache cache;
    cache.precalculate_matrices(ss.get_values(), set<double>{1, 3, 7}, 20);
    inference_pruner pruner(cache, &ss, nullptr, p_tree.get(), 1.0);

    auto actual = pruner.prune(rt, 20);

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

TEST_CASE("verify a complex credible interval")
{
    VectorXd v(100);
    v << 3.354626279025118532e-04,
        4.619598156723461941e-04,        6.320163806627527312e-04,        8.590462575532045947e-04,        1.160028999668890866e-03,
        1.556270993742927145e-03,        2.074271881760087182e-03,        2.746693614949631285e-03,        3.613423262850326730e-03,
        4.722712790075706685e-03,        6.132369489282128196e-03,        7.910959733930246832e-03,        1.013897644996755200e-02,
        1.290990763670983067e-02,        1.633113002139567962e-02,        2.052453933573554490e-02,        2.562681777342115663e-02,
        3.178923110480684489e-02,        3.917684398136177248e-02,        4.796704348995193407e-02,        5.834726929084104591e-02,
        7.051186479651316841e-02,        8.465798862252993384e-02,        1.009805593325092954e-01,        1.196662491095593917e-01,
        1.408865925476422143e-01,        1.647903336550449738e-01,        1.914951950146631665e-01,        2.210793147324477681e-01,
        2.535726555739642452e-01,        2.889487423304250568e-01,        3.271171235486371454e-01,        3.679169779747537561e-01,
        4.111122905071876166e-01,        4.563890040359159239e-01,        5.033545103172082369e-01,        5.515397744971644034e-01,
        6.004042952285054691e-01,        6.493439884813196894e-01,        6.977019528530103987e-01,        7.447819337574211884e-01,
        7.898641609263440388e-01,        8.322230966367535343e-01,        8.711465097193803464e-01,        9.059551911095097276e-01,
        9.360225578954148862e-01,        9.607933603497511577e-01,        9.798007140447178021e-01,        9.926807281297391761e-01,
        9.991840897954091805e-01,        9.991840897954091805e-01,        9.926807281297390650e-01,        9.798007140447175800e-01,
        9.607933603497510466e-01,        9.360225578954145531e-01,        9.059551911095095056e-01,        8.711465097193799023e-01,
        8.322230966367533123e-01,        7.898641609263437058e-01,        7.447819337574206333e-01,        6.977019528530101766e-01,
        6.493439884813191343e-01,        6.004042952285052470e-01,        5.515397744971638483e-01,        5.033545103172076818e-01,
        4.563890040359156464e-01,        4.111122905071871170e-01,        3.679169779747534785e-01,        3.271171235486367013e-01,
        2.889487423304248903e-01,        2.535726555739638566e-01,        2.210793147324474073e-01,        1.914951950146630000e-01,
        1.647903336550446685e-01,        1.408865925476421033e-01,        1.196662491095591696e-01,        1.009805593325090733e-01,
        8.465798862252993384e-02,        7.051186479651304351e-02,        5.834726929084096958e-02,        4.796704348995182998e-02,
        3.917684398136172391e-02,        3.178923110480677550e-02,        2.562681777342109765e-02,        2.052453933573554490e-02,
        1.633113002139563452e-02,        1.290990763670982026e-02,        1.013897644996752598e-02,        7.910959733930231219e-03,
        6.132369489282106512e-03,        4.722712790075693674e-03,        3.613423262850320224e-03,        2.746693614949621311e-03,
        2.074271881760083713e-03,        1.556270993742921724e-03,        1.160028999668888915e-03,        8.590462575532045947e-04,
        6.320163806627504544e-04,        4.619598156723453268e-04,        3.354626279025118532e-04;

    auto i = get_value(v, 1);    
    
    CHECK_EQ(doctest::Approx(0.255), i.credible_interval.first);
    CHECK_EQ(doctest::Approx(0.745), i.credible_interval.second);
}