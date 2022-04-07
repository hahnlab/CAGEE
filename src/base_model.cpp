#include <cmath>
#include <numeric>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>
#include <iomanip>

#include "doctest.h"

#include "base_model.h"
#include "rootdist_estimator.h"
#include "matrix_cache.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "optimizer_scorer.h"
#include "simulator.h"
#include "sigma.h"
#include "io.h"
#include "DiffMat.h"
#include "transcript_reconstructor.h"
#include "proportional_variance.h"

using namespace std;
namespace pv = proportional_variance;

extern mt19937 randomizer_engine;

base_model::base_model(sigma_squared* p_lambda, const vector<gene_transcript>* p_gene_families,
    error_model *p_error_model) :
    model(p_lambda, p_gene_families, p_error_model)
{

}

vector<size_t> build_reference_list(const vector<gene_transcript>& families)
{
    const size_t invalid_index = std::numeric_limits<size_t>::max();

    vector<size_t> reff;
    const int num_families = families.size();
    reff.resize(num_families, invalid_index);
    for (int i = 0; i < num_families; ++i) {
        if (reff[i] != invalid_index) continue;

        reff[i] = i;

        for (int j = i + 1; j < num_families; ++j) {
            if (reff[j] == invalid_index)
            {
                if (families[i].species_size_match(families[j]))
                {
                    reff[j] = i;
                }
            }
        }
    }

    return reff;
}

set<int> get_all_bounds(const vector<gene_transcript>& transcripts)
{
    vector<int> boundses(transcripts.size());
    transform(transcripts.begin(), transcripts.end(), boundses.begin(), [](const gene_transcript& gf) {
        return get_upper_bound(gf);
        });

    return set<int>(boundses.begin(), boundses.end());

}

inline double computational_space_prior(double val, const gamma_distribution<double>& prior)
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    return exp(val) * gammapdf(exp(val), prior);
#else
    return gammapdf(val, prior);
#endif

}
double compute_prior_likelihood(const vector<double>& partial_likelihood, const gene_transcript& t, const gamma_distribution<double>& prior)
{
    std::vector<double> full(partial_likelihood.size());
    double bound = get_upper_bound(t);
    for (size_t j = 0; j < partial_likelihood.size(); ++j) {
        full[j] = std::log(partial_likelihood[j]);

        double x = (double(j) + 0.5) * bound / (DISCRETIZATION_RANGE - 1);

        full[j] += log(computational_space_prior(x, prior));

        if (isnan(full[j]))
               full[j] = -numeric_limits<double>::infinity();
    }

#ifdef USE_MAX_PROBABILITY
    return *max_element(full.begin(), full.end()); // get max (CAFE's approach)
#else
    return accumulate(full.begin(), full.end(), 0.0, [](double a, double b) { return isinf(b) ? a : a+b; }); // sum over all sizes (Felsenstein's approach)
#endif
}

double base_model::infer_family_likelihoods(const user_data& ud, const sigma_squared *p_sigma, const gamma_distribution<double>& prior) {
    //TIMED_FUNC(timerObj);
    _monitor.Event_InferenceAttempt_Started();

    if (!p_sigma->is_valid())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    results.resize(ud.gene_families.size());
    std::vector<double> all_families_likelihood(ud.gene_families.size());

    matrix_cache calc;
    calc.precalculate_matrices(p_sigma->get_values(), get_all_bounds(ud.gene_families), ud.p_tree->get_branch_lengths());

    vector<vector<double>> partial_likelihoods(ud.gene_families.size());
#pragma omp parallel for
    for (int i = 0; i < ud.gene_families.size(); ++i) {
        if (references[i] == i)
            partial_likelihoods[i] = inference_prune(ud.gene_families.at(i), calc, p_sigma, _p_error_model, ud.p_tree, 1.0);
            // probabilities of various family sizes
    }

    // prune all the families with the same lambda
#pragma omp parallel for
    for (int i = 0; i < ud.gene_families.size(); ++i) {

        all_families_likelihood[i] = compute_prior_likelihood(partial_likelihoods[references[i]], ud.gene_families[i], prior);
        
        results[i] = family_info_stash(ud.gene_families.at(i).id(), 0.0, 0.0, 0.0, all_families_likelihood[i], false);
    }
    double final_likelihood = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;

    return final_likelihood;
}

void base_model::write_family_likelihoods(std::ostream& ost)
{
    ost << "#FamilyID\tLikelihood of Family" << endl;
    for (const auto& r : results)
    {
        ost << r.family_id << "\t" << r.posterior_probability << endl;
    }
}

sigma_optimizer_scorer* base_model::get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior)
{
    if (data.p_lambda != NULL)  // already have a lambda, nothing we want to optimize
        return nullptr;

    _p_sigma = initialize_search_sigma(data.p_lambda_tree, sample_groups);

    if (_p_error_model && !data.p_error_model)
    {
        return new sigma_optimizer_scorer(this, data, prior, _p_sigma, _p_error_model);
    }
    else
    {
        return new sigma_optimizer_scorer(this, data, prior, _p_sigma);
    }
}

#define EPSILON_RANGES

reconstruction* base_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    LOG(INFO) << "Starting reconstruction processes for Base model";

    auto result = new base_model_reconstruction();

    p_calc->precalculate_matrices(_p_sigma->get_values(), get_all_bounds(ud.gene_families), ud.p_tree->get_branch_lengths());


    for (size_t i = 0; i < ud.gene_families.size(); ++i)
    {
        auto &rc = result->_reconstructions[ud.gene_families[i].id()];
        ud.p_tree->apply_prefix_order([&rc](const clade* c) {
            rc[c] = 0;
            });
    }

//    pupko_reconstructor::pupko_data data(ud.gene_families.size(), ud.p_tree, ud.max_family_size, ud.max_root_family_size);
    transcript_reconstructor tr(_p_sigma, ud.p_tree, p_calc);

    for (int i = 0; i< ud.gene_families.size(); ++i)
    {
        result->_reconstructions[ud.gene_families[i].id()] = tr.reconstruct_gene_transcript(ud.gene_families[i]);
    }

//    size_t success = count_if(tr.v_all_node_Ls.begin(), tr.v_all_node_Ls.end(), [this, &ud](const clademap<Eigen::VectorXd>& L)
//        { 
//            return *max_element( L.at(ud.p_tree).begin(), L.at(ud.p_tree).end()) > 0;
//        });

//    if (success != ud.gene_families.size())
//    {
//        LOG(WARNING) << "Failed to reconstruct " << ud.gene_families.size() - success << " families" << endl;
//    }

    LOG(INFO) << "Done!\n";

    return result;
}

sigma_squared* base_model::get_simulation_lambda()
{
    return _p_sigma->multiply(simulation_lambda_multiplier);
}

double base_model_reconstruction::get_node_count(const gene_transcript& family, const clade *c) const
{
    if (c->is_leaf())
        return family.get_expression_value(c->get_taxon_name());

    if (_reconstructions.find(family.id()) == _reconstructions.end())
        throw runtime_error("Family " + family.id() + " not found in reconstruction");

    return _reconstructions.at(family.id()).at(c);
}

TEST_CASE("compute_prior_likelihood combines prior and inference correctly")
{
    gene_transcript gt;
    gt.set_expression_value("A", 12);
    gt.set_expression_value("B", 24);
    gamma_distribution<double> prior(0.75, 1/30.0);
    vector<double> inf{ 0.1, 0.2, 0.3};

    double actual = compute_prior_likelihood(inf, gt, prior);
#ifdef MODEL_GENE_EXPRESSION_LOGS
#ifdef USE_MAX_PROBABILITY
    CHECK_EQ(doctest::Approx(-35.7683), actual);
#else
    CHECK_EQ(doctest::Approx(-158.5437), actual);
#endif
#else
    CHECK_EQ(doctest::Approx(-5.584), actual);
#endif
}

TEST_CASE("base optimizer guesses sigma only")
{
    user_data ud;
    ud.p_lambda = NULL;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    base_model model(ud.p_lambda, NULL, NULL);

    unique_ptr<sigma_optimizer_scorer> opt(model.get_sigma_optimizer(ud, vector<string>(), std::gamma_distribution<double>(1, 2)));
    auto guesses = opt->initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(doctest::Approx(0.7071).epsilon(0.00001), guesses[0]);
    delete model.get_sigma();
}
