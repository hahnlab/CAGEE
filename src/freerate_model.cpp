#include <algorithm>
#include <numeric>
#include <iomanip>

#include <boost/math/tools/minima.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "freerate_model.h"
#include "user_data.h"
#include "inference_pruner.h"
#include "matrix_cache.h"
#include "prior.h"
#include "sigma.h"

using namespace std;

double compute_distribution_mean(const user_data& user_data);

void initialize_leaves(const user_data &ud, std::vector<clademap<optional_probabilities>> &thing, const matrix_cache &cache)
{
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade *c)
             {
        if (c->is_leaf())
        {
            for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
            {
                thing[i][c].reserve(cache.create_vector());
                thing[i][c].initialize(ud.gene_transcripts[i].get_expression_value(c->get_taxon_name()), ud.bounds);
            }
        } });
}


double get_log_likelihood(matrix_cache& cache, 
                            const vector<gene_transcript>& gene_transcripts,
                            const vector<double>& priors, 
                            vector<clademap<optional_probabilities>>& probs,
                            const clade* c, 
                            const boundaries& bounds,
                            double sigsqd)
{
    std::vector<double> all_transcripts_likelihood(gene_transcripts.size());
    vector<vector<double>> partial_likelihoods(gene_transcripts.size());
    sigma_squared ss(sigsqd);
    set<double> branch_lengths;
    c->apply_to_descendants([&](const clade* d) { branch_lengths.insert(d->get_branch_length()); });
    cache.precalculate_matrices(ss.get_values(),  branch_lengths, bounds);

    for(auto& m : probs) {
        m[c].reserve(cache.create_vector());
    }

#pragma omp parallel for
    for (int i = 0; i < (int)gene_transcripts.size(); ++i) {    
        compute_node_probability(c, gene_transcripts[i], nullptr, nullptr, probs[i], &ss, cache, bounds);
        auto& p = probs[i][c].probabilities();
        partial_likelihoods[i] = vector<double>(p.begin(), p.end());
    }
    
    // At this point partial_likelihoods has all the data we need to compute the parent of c
    std::transform(partial_likelihoods.begin(), partial_likelihoods.end(), all_transcripts_likelihood.begin(), [&priors](const vector<double>& a) {
        return log(compute_prior_likelihood(a, priors));
    });

    double result = std::accumulate(all_transcripts_likelihood.begin(), all_transcripts_likelihood.end(), 0.0);
    if (std::isinf(result))
        result = 10000;
    return result;
}

freerate_model::freerate_model() :
    model(nullptr, nullptr, nullptr)
{
}

double freerate_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma)
{
    double distmean = compute_distribution_mean(ud);
    matrix_cache cache;
    vector<double> priors;
    try
    {
        priors = get_priors(cache, ud.bounds, ud.p_prior);
    }
    catch (std::domain_error& e)
    {
        LOG(DEBUG) << e.what();
        LOG(WARNING) << "Prior not valid for this sigma and data set";
        return -log(0);
    }
    
    vector<clademap<optional_probabilities>> probs(ud.gene_transcripts.size());
    initialize_leaves(ud, probs, cache);

    clademap<sigma_squared*> sigmas;
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* c) {
        if (!c->is_leaf())
        {
            using boost::math::tools::brent_find_minima;

            std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
                return -get_log_likelihood(cache, ud.gene_transcripts, priors, probs, c, ud.bounds, sigsqd);
            }, 0.0, distmean*100, 5);

            LOG(INFO) << "Node " << c->get_ape_index() << " Sigma^2:" << r.first;
            LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << r.second;

            sigmas[c] = new sigma_squared(r.first);
        }
    });

    _p_sigma = sigmas[ud.p_tree];
    return 0.0;
}

sigma_optimizer_scorer* freerate_model::get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups)
{
    return nullptr;
}

class freerate_reconstruction : public reconstruction
{  
public:
    freerate_reconstruction(const user_data& ud, matrix_cache* p_calc) : reconstruction(ud.gene_transcripts) {};
    virtual ~freerate_reconstruction() {};

    virtual node_reconstruction get_internal_node_value(const gene_transcript& gf, const clade* c) const override
    {
        return node_reconstruction();
    }
};

reconstruction* freerate_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    return new freerate_reconstruction(ud, p_calc);
}


TEST_CASE("freerate_model")
{
    new freerate_model();
}   

TEST_CASE("freerate_model optimizes a branch length")
{
    freerate_model m;

    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("uniform", 0.0, 1.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    m.infer_transcript_likelihoods(ud, nullptr);
}   

TEST_CASE("get_log_likelihood")
{
    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("uniform", 0.0, 100.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    matrix_cache cache;
    auto priors = get_priors(cache, ud.bounds, ud.p_prior);
    vector<clademap<optional_probabilities>> thing(1);
    auto init_func = [&](const clade* node) { thing[0][node].reserve(cache.create_vector()); };
    ud.p_tree->apply_to_descendants(init_func);
    init_func(ud.p_tree);

    initialize_leaves(ud, thing, cache);
    auto lnl = get_log_likelihood(cache, ud.gene_transcripts, priors, thing, ud.p_tree, ud.bounds, 1.0);
    CHECK_EQ(doctest::Approx(-7.98967), lnl);
}

TEST_CASE("get_log_likelihood can calculate likelihoods for subtrees")
{
    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("uniform", 0.0, 100.0);
    ud.p_tree = parse_newick("((A:1,B:3):7,(C:11,D:17):23);");
    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);
    ud.gene_transcripts[0].set_expression_value("C", 3);
    ud.gene_transcripts[0].set_expression_value("D", 4);

    matrix_cache cache;
    auto priors = get_priors(cache, ud.bounds, ud.p_prior);
    vector<clademap<optional_probabilities>> thing(1);
    auto init_func = [&](const clade* node) { thing[0][node].reserve(cache.create_vector()); };
    ud.p_tree->apply_to_descendants(init_func);
    init_func(ud.p_tree);

    initialize_leaves(ud, thing, cache);

    auto lnl = get_log_likelihood(cache, ud.gene_transcripts, priors, thing, ud.p_tree->find_descendant("AB"), ud.bounds, 1.0);
    CHECK_EQ(doctest::Approx(-8.21231), lnl);

    lnl = get_log_likelihood(cache, ud.gene_transcripts, priors, thing, ud.p_tree, ud.bounds, 1.0);
    CHECK_EQ(doctest::Approx(-6.96722), lnl);
}

TEST_CASE("get_log_likelihood can calculate likelihoods for subtrees with different priors")
{
    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("uniform", 0.0, 100.0);
    ud.p_tree = parse_newick("((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;");

    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 0.5);
    ud.gene_transcripts[0].set_expression_value("B", 1);
    ud.gene_transcripts[0].set_expression_value("C", 1.5);
    ud.gene_transcripts[0].set_expression_value("D", 2);
    ud.gene_transcripts[0].set_expression_value("E", 2.5);

    freerate_model m;
    m.infer_transcript_likelihoods(ud, nullptr);
}
