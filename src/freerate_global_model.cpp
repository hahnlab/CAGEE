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


freerate_global_model::freerate_global_model(bool values_are_ratios) :
    model(nullptr, nullptr, nullptr), _values_are_ratios(values_are_ratios)
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
    return -lnl;
}

void optimize_root_sigma(const user_data& ud, const clademap<ASV>& probs, const clademap<prior>& priors, 
    clademap<std::pair<double, double>>& _sigmas, matrix_cache& cache)
{
    using boost::math::tools::brent_find_minima;

    auto priors_by_bin = get_priors(cache, ud.bounds, &priors.at(ud.p_tree));
    std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
        return -get_log_likelihood(probs.at(ud.p_tree), priors_by_bin);

    }, 0.0, _sigmas[ud.p_tree].first*100.0, 10);
    _sigmas[ud.p_tree] = r;
    LOG(INFO) << "Node " << ud.p_tree->get_ape_index() << " Sigma^2:" << r.first;
    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << r.second;
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
void freerate_global_model::optimize_sigmas(const user_data& ud, const clademap<prior>& priors)
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

                std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
                    return compute_node_likelihood(d, probs, _sigmas, priors_by_bin, sigsqd, cache, ud);
                }, 0.0, _sigmas[p].first*100.0, 10);
                _sigmas[d] = r;
                LOG(INFO) << "Node " << d->get_ape_index() << " Sigma^2:" << r.first;
                LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << r.second;
            });
        }
    });

    optimize_root_sigma(ud, probs, priors, _sigmas, cache);
}

double freerate_global_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma)
{
    auto priors = compute_tree_priors(ud.p_tree, ud.gene_transcripts, _values_are_ratios);
    for (auto& a : priors)
    {
        LOG(DEBUG) << "Node " << a.first->get_ape_index() << " Prior: " << a.second;
    }
    double distmean = compute_distribution_mean(ud);

    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* p) {
        _sigmas[p] = pair<double,double>(distmean,-1);
    });
    LOG(DEBUG) << "Initial sigmas set to " << distmean;

    double score = 100;
    for (int i = 0; i<20; ++i)
    {
        LOG(INFO) << "Iteration " << i+1;
        optimize_sigmas(ud, priors);
        if (abs(score - _sigmas[ud.p_tree].second) < 10)
            break;
        score = _sigmas[ud.p_tree].second;
    }
    _p_sigma = new sigma_squared(_sigmas[ud.p_tree].first);
    return _sigmas[ud.p_tree].second;
}

sigma_optimizer_scorer* freerate_global_model::get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups)
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

reconstruction* freerate_global_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    auto result = new freerate_reconstruction(ud, p_calc);

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
        ost << a.first << ":" << a.second << "\n";
    }
}


TEST_CASE("freerate_model")
{
    new freerate_global_model(false);
}

TEST_CASE("freerate_model optimizes a branch length")
{
    freerate_global_model m(false);

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

    freerate_global_model m(false);
    m.infer_transcript_likelihoods(ud, nullptr);

    ostringstream ost;
    m.write_extra_vital_statistics(ost);
    CHECK_STREAM_CONTAINS(ost, "Computed sigma2 by node:");
    CHECK_STREAM_CONTAINS(ost, "6:0.000766305");
    CHECK_STREAM_CONTAINS(ost, "7:0.0590542");
    CHECK_STREAM_CONTAINS(ost, "8:0.0590542");
    CHECK_STREAM_CONTAINS(ost, "9:0.595233");
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
