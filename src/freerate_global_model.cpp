#include <algorithm>
#include <numeric>
#include <iomanip>

#include <boost/math/tools/minima.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "freerate_global_model.h"
#include "user_data.h"
#include "inference_pruner.h"
#include "matrix_cache.h"
#include "prior.h"
#include "sigma.h"
#include "optimizer_scorer.h"   // for definition of compute_distribution_mean

using namespace std;
using namespace Eigen;

typedef vector<optional_probabilities> ASV;

ASV get_taxon_asv(std::string taxon, const user_data &ud, const matrix_cache &cache)
{
    ASV result;
    result.resize(ud.gene_transcripts.size());
    for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
    {
        result[i].reserve(cache.create_vector());
        result[i].initialize(ud.gene_transcripts[i].get_expression_value(taxon), ud.bounds);
    }
    return result;
}

const clade * get_sibling(const clade *c)
{
    if (c->is_root())
        return nullptr;
    auto p = c->get_parent();
    auto sib = find_if(p->descendant_begin(), p->descendant_end(), [&](clade *d) {
        return d != c;
    });
    return sib == p->descendant_end() ? nullptr : *sib;
}

double get_log_likelihood(const ASV& partial_likelihoods, const std::vector<double>& priors)
{
    std::vector<double> all_transcripts_likelihood(partial_likelihoods.size());

    transform(partial_likelihoods.begin(), partial_likelihoods.end(), all_transcripts_likelihood.begin(), [&](const optional_probabilities& a) {
        vector<double> t(a.probabilities().size());
        copy(a.probabilities().begin(), a.probabilities().end(), t.begin());
        return log(compute_prior_likelihood(t, priors));
    });

    double result = std::accumulate(all_transcripts_likelihood.begin(), all_transcripts_likelihood.end(), 0.0);
    if (std::isinf(result))
        result = 10000;
    return result;
}


freerate_global_model::freerate_global_model(bool values_are_ratios, std::string initial_values, bool initial_values_are_weights) :
    model(nullptr, nullptr, nullptr), _values_are_ratios(values_are_ratios), _initial_values(initial_values), _initial_values_are_weights(initial_values_are_weights)
{
}

// This has a lot of elements in common with compute_node_probability. Can they be refactored?
double compute_node_likelihood(const clade* d,
                    clademap<ASV>& probs, clademap<std::pair<double, double>>& sigmas,
                    const vector<double>& priors_by_bin, double sigsqd, matrix_cache& cache, const user_data& ud)
{
    auto p = d->get_parent();
    auto sib = get_sibling(d);
    cache.precalculate_matrices({ sigsqd, sigmas[sib].first }, { d->get_branch_length(), sib->get_branch_length() }, ud.bounds);

    probs[p] = ASV(ud.gene_transcripts.size());

    for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i)
    {
        probs[p][i].reserve(cache.create_vector());
    }

#pragma omp parallel for
    for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i)
    {
        if (probs[d][i].hasValue() || probs[sib][i].hasValue())
        {
            probs[p][i].setOne();
        }

        if (probs[d][i].hasValue())
        {
            const MatrixXd& m = cache.get_matrix(d->get_branch_length(), sigsqd);
            VectorXd result = m * probs[d][i].probabilities();
            probs[p][i].multiply_elements(result);
        }
        if (probs[sib][i].hasValue())
        {
            const MatrixXd& m = cache.get_matrix(sib->get_branch_length(), sigmas[sib].first);
            VectorXd result = m * probs[sib][i].probabilities();
            probs[p][i].multiply_elements(result);
        }
    }

    auto lnl = get_log_likelihood(probs[p], priors_by_bin);
    LOG(DEBUG) << "Node " << d->get_ape_index() << " Sigma^2:" << sigsqd << " Score (-lnL): " << lnl;
    return -lnl;
}

/* Init all branches with a sigma
To calculate AF:
        ASV[B] = get_ASV(childB branch length, sigma[B], childB init)
        optimize ss:
                ASV[A] = get_ASV(childA branch length, ss, childA init)
                ASV[F] = ASV[AF] * ASV[BF]
                likelihood = ASV[F] * prior


To calculate BF:
        ASV[BF] = optimized value of get_ASV(childA bl, childB bl, sigma, ASV[AF], childB init)
To calculate GF:
        ASV[GC] = get_ASV(childC branch length, childF branch length, sigma_init, childC init, ASV[AF] * ASV[BF])
        ASV[GF] = optimized value of get_ASV(childC bl, childF bl, sigma, ASV[GC], childC init)
 */
double freerate_global_model::optimize_sigmas(const user_data& ud, const clademap<prior>& priors)
{
    matrix_cache cache;

    clademap<ASV> probs;

    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* p) {
        if (p->is_leaf())
        {
            probs[p] = get_taxon_asv(p->get_taxon_name(), ud, cache);
        }
    });
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* p) {
        if (!p->is_leaf())
        {
            auto priors_by_bin = get_priors(cache, ud.bounds, &priors.at(p));
            for_each(p->descendant_begin(), p->descendant_end(), [&](const clade* d) {
                using boost::math::tools::brent_find_minima;

                LOG(DEBUG) << "Optimizing node " << d->get_ape_index() << " with sibling sigma^2 " << _sigmas[get_sibling(d)].first;
                std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
                    return compute_node_likelihood(d, probs, _sigmas, priors_by_bin, sigsqd, cache, ud);
                }, 0.0, 20.0, 10);
                _sigmas[d] = r;
                LOG(INFO) << "Sigma^2:" << r.first << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << r.second;
            });
        }
    });

    auto root_priors = get_priors(cache, ud.bounds, &priors.at(ud.p_tree));
    auto result = -get_log_likelihood(probs.at(ud.p_tree), root_priors);
    LOG(INFO) << "Final tree score (-lnL): " << std::setw(15) << std::setprecision(14) << result;
    return result;
}

void freerate_global_model::initialize_sigmas(const clade* p_tree, double distmean)
{
    unique_ptr<clade> p_weight_tree;
    if (!_initial_values.empty())
        p_weight_tree.reset(parse_newick(_initial_values));


    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&](const clade* p) {
        if (p->is_root() || !p_weight_tree)
            _sigmas[p] = pair<double,double>(distmean,-1);
        else 
        {
            auto weight = 0.0;
            auto matched_clade = p_weight_tree->find_descendant(p->get_taxon_name());
            if (matched_clade)
            {   
                weight = matched_clade->get_branch_length();
                LOG(DEBUG) << "Found a weight for " << p->get_taxon_name() << ". It is " << weight;
            }
            else
            {
                vector<string> descendants(distance(p->descendant_begin(), p->descendant_end()));
                transform(p->descendant_begin(), p->descendant_end(), descendants.begin(), [](const clade* c) { return c->get_taxon_name(); });
                auto d1 = p_weight_tree->find_descendant(descendants[0]);
                auto d2 = p_weight_tree->find_descendant(descendants[1]);
                weight = distance(d1, d2);

                LOG(DEBUG) << "There is no weight for " << p->get_taxon_name() << ". Computed distance between " << descendants[0] << " and " << descendants[1] << " is " << weight;
            }
            if (_initial_values_are_weights)
                weight = distmean * weight;
            _sigmas[p] =  pair<double,double>(weight,-1);  
        }
        string name = p->is_leaf() ? " (" + p->get_taxon_name() + ")" : "";
        LOG(DEBUG) << "Initial sigma for Node " << p->get_ape_index() << name << " set to " << _sigmas[p].first;

    });
    LOG(DEBUG) << "Initial sigmas set to " << distmean;
}

double freerate_global_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma)
{
    _root_ape_index = ud.p_tree->get_ape_index();

    auto priors = compute_tree_priors(ud.p_tree, ud.gene_transcripts, _values_are_ratios);
    for (auto& a : priors)
    {
        LOG(DEBUG) << "Node " << a.first->get_ape_index() << " Prior: " << a.second;
    }
    double distmean = compute_distribution_mean(ud);
    LOG(DEBUG) << "Distribution mean: " << distmean;
    
    initialize_sigmas(ud.p_tree, distmean);

    double score = 100;
    for (int i = 0; i<20; ++i)
    {
        LOG(INFO) << "Iteration " << i+1;
        _monitor.Event_InferenceAttempt_Started();
        double lnl = optimize_sigmas(ud, priors);
        if (abs(score - lnl) < 10)
            break;
        score = lnl;
    }
    return score;
}

sigma_optimizer_scorer* freerate_global_model::get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups)
{
    return nullptr;
}

class freerate_global_reconstruction : public reconstruction
{
public:
    freerate_global_reconstruction(const user_data& ud, matrix_cache* p_calc) : reconstruction(ud.gene_transcripts) {};
    virtual ~freerate_global_reconstruction() {};

    std::map<const gene_transcript*, clademap<node_reconstruction>> _reconstructions;

    virtual node_reconstruction get_internal_node_value(const gene_transcript& gf, const clade* c) const override
    {
        return node_reconstruction();
    }
};

void compute_transcript_likelihoods(const clade* c, int i, const matrix_cache& cache, clademap<ASV>& probs, const clademap<std::pair<double, double>>& sigmas, size_t transcript_count)
{
    if (c->is_root())
        return;

    auto p = c->get_parent();
    probs[p] = ASV(transcript_count);
    for (int i = 0; i < (int)transcript_count; ++i)
    {
        probs[p][i].reserve(cache.create_vector());
    }

    auto sib = get_sibling(c);

    if (probs[c][i].hasValue() || probs[sib][i].hasValue())
    {
        probs[p][i].setOne();
    }

    if (probs[c][i].hasValue())
    {
        const MatrixXd& m = cache.get_matrix(c->get_branch_length(), sigmas.at(c).first);
        VectorXd result = m * probs[c][i].probabilities();
        probs[p][i].multiply_elements(result);
    }
    if (probs[sib][i].hasValue())
    {
        const MatrixXd& m = cache.get_matrix(sib->get_branch_length(), sigmas.at(sib).first);
        VectorXd result = m * probs[sib][i].probabilities();
        probs[p][i].multiply_elements(result);
    }
}

extern node_reconstruction get_value(const VectorXd& likelihood, boundaries bounds);

reconstruction* freerate_global_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    auto result = new freerate_global_reconstruction(ud, p_calc);

    clademap<ASV> probs;
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* p) {
        if (p->is_leaf())
        {
            probs[p] = get_taxon_asv(p->get_taxon_name(), ud, *p_calc);
        }
    });

    for (size_t i = 0; i< ud.gene_transcripts.size(); ++i)
    {
        auto compute_func = [&](const clade* c) { compute_transcript_likelihoods(c, i, *p_calc, probs, _sigmas, ud.gene_transcripts.size()); };
        for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), compute_func);
    }

    for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
    {
        auto &rc = result->_reconstructions[&ud.gene_transcripts[i]];
        ud.p_tree->apply_prefix_order([&rc](const clade* c) {
            rc[c].most_likely_value = 0;
            });
    }

    clademap<node_reconstruction> reconstruction;
    for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
    {
        cladevector internal_nodes_with_values;
        copy_if(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), back_inserter(internal_nodes_with_values), [this, &probs, i](const clade* c)
        {
            return !c->is_leaf() && probs.at(c)[i].hasValue();
        });

        for (auto& c : internal_nodes_with_values)
        {
            reconstruction[c] = get_value(probs[c][i].probabilities(), ud.bounds);
        }
        result->_reconstructions[&ud.gene_transcripts[i]] = reconstruction;
    }
    return result;
}

void freerate_global_model::write_extra_vital_statistics(std::ostream& ost)
{
    ost << "Computed sigma2 by node:\n";
    vector<pair<int,double>> sorted;
    transform(_sigmas.begin(), _sigmas.end(), back_inserter(sorted), [](const pair<const clade*, pair<double, double>>& a) { return make_pair(a.first->get_ape_index(), a.second.first); });
    sort(sorted.begin(), sorted.end(), [](const pair<int, double>& a, const pair<int, double>& b) { return a.first < b.first; });
    for (auto& a : sorted)
    {
        if (a.first == _root_ape_index)
            continue;

        ost << a.first << ":" << a.second << "\n";
    }
}


TEST_CASE("freerate_model")
{
    new freerate_global_model(false, "", false);
}

TEST_CASE("freerate_model optimizes a branch length")
{
    freerate_global_model m(false, "", false);

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
    ASV asv(1);
    Eigen::VectorXd v(5);
    v.setOnes();
    asv[0].set(v);
    auto lnl = get_log_likelihood(asv, vector<double>({ 1.0, 1.0, 1.0, 1.0, 1.0 }));
    CHECK_EQ(doctest::Approx(1.60944), lnl);
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("write_extra_vital_statistics writes a sigma for each internal node")
{
//    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
//    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToStandardOutput, "true");
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

    freerate_global_model m(false, "", false);
    m.infer_transcript_likelihoods(ud, nullptr);

    ostringstream ost;
    m.write_extra_vital_statistics(ost);
    CHECK_STREAM_CONTAINS(ost, "Computed sigma2 by node:");
    CHECK_STREAM_CONTAINS(ost, "7:0.000817127");
    CHECK_STREAM_CONTAINS(ost, "8:2.23993");
    CHECK_STREAM_CONTAINS(ost, "9:0.595146");

    // Check that the root is not included
    CHECK_MESSAGE(ost.str().find("6:") == std::string::npos, ost.str());
    //el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToStandardOutput, "false");
}

TEST_CASE("get_parent_likelihood")
{
    user_data ud;
    ud.p_sigma = nullptr;
    ud.p_prior = new prior("uniform", 0.0, 100.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    matrix_cache cache;
    clademap<ASV> probs;
    const clade* p = ud.p_tree;
    const clade* sib = get_sibling(*p->descendant_begin());
    const clade* d = *p->descendant_begin();
    probs[d] = get_taxon_asv(d->get_taxon_name(), ud, cache);
    probs[sib] = get_taxon_asv(sib->get_taxon_name(), ud, cache);
    vector<double> priors_by_bin = get_priors(cache, ud.bounds, ud.p_prior);
    double sigsqd = 1.0;
    clademap<std::pair<double, double>> sigmas;
    sigmas[sib]= pair<double,double>(1.5, -1.0);

    double likelihood = compute_node_likelihood(d, probs, sigmas, priors_by_bin, sigsqd, cache, ud);

    CHECK_EQ(doctest::Approx(8.07594), likelihood);
}

TEST_CASE("compute_transcript_likelihoods")
{
    // Setup user_data and tree
    user_data ud;
    ud.p_sigma = nullptr;
    ud.p_prior = new prior("uniform", 0.0, 10.0);
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.bounds = boundaries(0, 5);

    // Add a gene transcript with values for A and B
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1.0);
    ud.gene_transcripts[0].set_expression_value("B", 2.0);

    // Prepare matrix_cache and ASV probabilities
    matrix_cache cache;
    cache.precalculate_matrices({1,1}, ud.p_tree->get_branch_lengths(), ud.bounds);
    clademap<ASV> probs;
    const clade* p = ud.p_tree;
    const clade* d = *p->descendant_begin();
    const clade* sib = get_sibling(d);

    probs[d] = get_taxon_asv(d->get_taxon_name(), ud, cache);
    probs[sib] = get_taxon_asv(sib->get_taxon_name(), ud, cache);
    probs[p] = ASV(ud.gene_transcripts.size());
    for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i)
    {
        probs[p][i].reserve(cache.create_vector());
    }
    
    // Setup sigmas for both children
    clademap<std::pair<double, double>> sigmas;
    sigmas[d] = std::make_pair(1.0, -1.0);
    sigmas[sib] = std::make_pair(1.0, -1.0);

    // Compute transcript likelihoods for gene 0
    compute_transcript_likelihoods(d, 0, cache, probs, sigmas, ud.gene_transcripts.size());

    // Check that parent probabilities are set (should not be empty)
    CHECK(probs[p][0].hasValue());
    // Optionally, check that the probabilities vector is not all zeros
    auto v = probs[p][0].probabilities();
    bool nonzero = false;
    for (int i = 0; i < v.size(); ++i) {
        if (v[i] != 0.0) {
            nonzero = true;
            break;
        }
    }
    CHECK(nonzero);

}

TEST_CASE("reconstruct_ancestral_states")
{
    user_data ud;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);
    ud.bounds = boundaries(0, 20);
    sigma_squared sl(0.05);
    ud.p_sigma = &sl;

    std::vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 4);
    
    freerate_global_model model(false,"", false);
    model.initialize_sigmas(ud.p_tree, 0.5);

    matrix_cache calc;
    calc.precalculate_matrices({1,.5}, ud.p_tree->get_branch_lengths(), ud.bounds);

    auto rec = dynamic_cast<freerate_global_reconstruction*>(model.reconstruct_ancestral_states(ud, &calc));

    CHECK_EQ(1, rec->_reconstructions.size());
}