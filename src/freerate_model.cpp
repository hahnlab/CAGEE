#include <algorithm>
#include <numeric>

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

double get_log_likelihood(matrix_cache& cache, 
                            const vector<gene_transcript>& gene_transcripts,
                            const vector<double>& priors, 
                            const clade* c, 
                            const boundaries& bounds,
                            double sigsqd)
{
    cout << "Candidate: " << sigsqd << "\n";
    cout << "Log likelihood:" << -log(sigsqd) << "\n";
    return log(sigsqd);
    std::vector<double> all_transcripts_likelihood(gene_transcripts.size());
    vector<vector<double>> partial_likelihoods(gene_transcripts.size());
    sigma_squared ss(sigsqd);
    cache.precalculate_matrices(ss.get_values(),  c->get_branch_lengths(), bounds);
    inference_pruner pruner(cache, &ss, nullptr, nullptr, c, bounds);
    transform(gene_transcripts.begin(), gene_transcripts.end(), partial_likelihoods.begin(), [&](const gene_transcript& gt) {
        return pruner.prune(gt);
    }); 
    
    transform(partial_likelihoods.begin(), partial_likelihoods.end(), all_transcripts_likelihood.begin(), [&priors](const vector<double>& a) {
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
    clademap<optional_probabilities> probabilities;

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
    
    for_each(ud.p_tree->reverse_level_begin(), ud.p_tree->reverse_level_end(), [&](const clade* c) {
        if (!c->is_leaf())
        {
            using boost::math::tools::brent_find_minima;

            std::pair<double, double> r = brent_find_minima([&](double sigsqd) {
                return -get_log_likelihood(cache, ud.gene_transcripts, priors, c, ud.bounds, sigsqd);
            }, 0.0, distmean*100, 5);
            cout << "Optimal sigma^2 for clade " << c->get_taxon_name() << " is " << r.first << " with likelihood " << r.second;

        }
    });

    return 0.0;
}

sigma_optimizer_scorer* freerate_model::get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups)
{
    return nullptr;
}

reconstruction* freerate_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    return nullptr;
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

TEST_CASE("sfdsf")
{

}