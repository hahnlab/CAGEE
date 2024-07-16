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
#include "optimizer_scorer.h"   // for definition of compute_distribution_mean

using namespace std;

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

clademap<prior> node_priors(const clade* p_tree, const vector<gene_transcript>& gene_transcripts, bool ratios)
{
    clademap<vector<double>> averages;
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&](const clade* c) {
        if (!c->is_leaf())
        {
            averages[c] = vector<double>(gene_transcripts.size());
            int child_count = 0;
            c->apply_to_descendants([&](const clade* d) 
            { 
                if (averages.find(d) != averages.end())
                {
                    // child has already been processed, so add its (already averaged) values to the current node
                    transform(averages[d].begin(), averages[d].end(), averages[c].begin(), averages[c].begin(), std::plus<double>());
                }
                else
                {
                    // child is a leaf, so add its transcript values to the current node
                    transform(gene_transcripts.begin(), gene_transcripts.end(), averages[c].begin(), averages[c].begin(), [d, ratios](const gene_transcript &gf, double val) {
                        double newval = gf.get_expression_value(d->get_taxon_name());
                        return val + (ratios ? newval : exp(newval));
                    });
                }
                child_count++;  
            });
            transform(averages[c].begin(), averages[c].end(), averages[c].begin(), [child_count](double val) { return val / child_count; });  
        }
    });
    clademap<prior> result;
    for(auto& a : averages)
    {
        string dist = ratios ? "normal" : "gamma";
        result[a.first] = estimate_distribution(dist, a.second);
    }
    return result;
}

double get_log_likelihood(matrix_cache& cache, 
                            const vector<gene_transcript>& gene_transcripts,
                            prior p,
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
    auto priors = get_priors(cache, bounds, &p);
    
    std::transform(partial_likelihoods.begin(), partial_likelihoods.end(), all_transcripts_likelihood.begin(), [&priors](const vector<double>& a) {
        return log(compute_prior_likelihood(a, priors));
    });

    double result = std::accumulate(all_transcripts_likelihood.begin(), all_transcripts_likelihood.end(), 0.0);
    if (std::isinf(result))
        result = 10000;
    return result;
}

freerate_model::freerate_model(bool values_are_ratios) :
    model(nullptr, nullptr, nullptr), _values_are_ratios(values_are_ratios)
{
}

double freerate_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma)
{
    matrix_cache cache;
    
    vector<clademap<optional_probabilities>> probs(ud.gene_transcripts.size());
    initialize_leaves(ud, probs, cache);
    auto priors = node_priors(ud.p_tree, ud.gene_transcripts, _values_are_ratios);
    for (auto& a : priors)
    {
        LOG(DEBUG) << "Node " << a.first->get_ape_index() << " Prior: " << a.second;
    }   
    double distmean = compute_distribution_mean(ud);
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* c) {
        if (!c->is_leaf())
        {
            using boost::math::tools::brent_find_minima;

            std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
                _monitor.Event_InferenceAttempt_Started();
                return -get_log_likelihood(cache, ud.gene_transcripts, priors[c], probs, c, ud.bounds, sigsqd);
            }, 0.0, distmean*100, 5);

            LOG(INFO) << "Node " << c->get_ape_index() << " Sigma^2:" << r.first;
            LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << r.second;

            _sigmas[c] = r;
        }
    });

    _p_sigma = new sigma_squared(_sigmas[ud.p_tree].first);
    return _sigmas[ud.p_tree].second;
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

void freerate_model::write_extra_vital_statistics(std::ostream& ost) 
{
    ost << "Computed sigma2 by node:\n";
    vector<pair<int,double>> sorted;
    transform(_sigmas.begin(), _sigmas.end(), back_inserter(sorted), [](const pair<const clade*, pair<double, double>>& a) { return make_pair(a.first->get_ape_index(), a.second.first); });
    sort(sorted.begin(), sorted.end(), [](const pair<int, double>& a, const pair<int, double>& b) { return a.first < b.first; });
    for (auto& a : sorted)
    {
        ost << a.first << ":" << a.second << "\n";
    }
}


TEST_CASE("freerate_model")
{
    new freerate_model(false);
}   

TEST_CASE("freerate_model optimizes a branch length")
{
    freerate_model m(false);

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
    auto lnl = get_log_likelihood(cache, ud.gene_transcripts, *ud.p_prior, thing, ud.p_tree, ud.bounds, 1.0);
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

    auto lnl = get_log_likelihood(cache, ud.gene_transcripts, *ud.p_prior, thing, ud.p_tree->find_descendant("AB"), ud.bounds, 1.0);
    CHECK_EQ(doctest::Approx(-8.21231), lnl);

    lnl = get_log_likelihood(cache, ud.gene_transcripts, *ud.p_prior, thing, ud.p_tree, ud.bounds, 1.0);
    CHECK_EQ(doctest::Approx(-6.96722), lnl);
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

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

    freerate_model m(false);
    m.infer_transcript_likelihoods(ud, nullptr);
}

TEST_CASE("write_extra_vital_statistics writes a sigma for each internal node")
{
    user_data ud;
    ud.p_sigma = NULL;
    ud.p_prior = new prior("gamma", 0.375, 1600.0);
    ud.p_tree = parse_newick("((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;");

    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 0.5);
    ud.gene_transcripts[0].set_expression_value("B", 1);
    ud.gene_transcripts[0].set_expression_value("C", 1.5);
    ud.gene_transcripts[0].set_expression_value("D", 2);
    ud.gene_transcripts[0].set_expression_value("E", 2.5);

    freerate_model m(false);
    m.infer_transcript_likelihoods(ud, nullptr);

    ostringstream ost;
    m.write_extra_vital_statistics(ost);
    CHECK_STREAM_CONTAINS(ost, "Computed sigma2 by node:");
    CHECK_STREAM_CONTAINS(ost, "6:42.8082");
    CHECK_STREAM_CONTAINS(ost, "7:42.8082");
    CHECK_STREAM_CONTAINS(ost, "8:42.8082");
    CHECK_STREAM_CONTAINS(ost, "9:42.8082");
}

TEST_CASE("node_priors")
{
    unique_ptr<clade> t(parse_newick("((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;"));

    vector<gene_transcript> gene_transcripts;
    gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    gene_transcripts.push_back(gene_transcript("TestFamily2", "", ""));
    gene_transcripts[0].set_expression_value("A", 0.5);
    gene_transcripts[0].set_expression_value("B", 1);
    gene_transcripts[0].set_expression_value("C", 1.5);
    gene_transcripts[0].set_expression_value("D", 2);
    gene_transcripts[0].set_expression_value("E", 3);
    gene_transcripts[1].set_expression_value("A", 8);
    gene_transcripts[1].set_expression_value("B", 4);
    gene_transcripts[1].set_expression_value("C", 4.5);
    gene_transcripts[1].set_expression_value("D", 3.5);
    gene_transcripts[1].set_expression_value("E", 2.5);

    auto result = node_priors(t.get(), gene_transcripts, false);
    CHECK_EQ(4, result.size());
    ostringstream ost;
    for(auto& a : result)
    {
        ost << a.first->get_ape_index() << ":" << a.second << "|";
    }   
    CHECK_STREAM_CONTAINS(ost, "6:gamma:1.08613:194.18");
    CHECK_STREAM_CONTAINS(ost, "7:gamma:16.6708:1.09132");
    CHECK_STREAM_CONTAINS(ost, "8:gamma:1.01672:396.977");
    CHECK_STREAM_CONTAINS(ost, "9:gamma:1.00577:755.62");
}
