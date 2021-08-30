#include <cmath>
#include <numeric>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>
#include <iomanip>

#include "base_model.h"
#include "gene_family_reconstructor.h"
#include "matrix_cache.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "optimizer_scorer.h"
#include "simulator.h"
#include "sigma.h"
#include "io.h"
#include "DiffMat.h"

using namespace std;

extern mt19937 randomizer_engine;

base_model::base_model(sigma* p_lambda, const vector<gene_transcript>* p_gene_families,
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

set<pair<double, double>> get_all_bounds(const vector<gene_transcript>& transcripts)
{
    vector<pair<double, double>> boundses(transcripts.size());
    transform(transcripts.begin(), transcripts.end(), boundses.begin(), [](const gene_transcript& gf) {
        return bounds(gf);
        });

    return set<pair<double, double>>(boundses.begin(), boundses.end());

}

double base_model::infer_family_likelihoods(const user_data& ud, const sigma *p_sigma) {
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
    calc.precalculate_matrices(p_sigma->get_lambdas(), get_all_bounds(ud.gene_families), ud.p_tree->get_branch_lengths());

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

        auto& partial_likelihood = partial_likelihoods[references[i]];
        std::vector<double> full(partial_likelihood.size());

        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = ud.p_prior->compute(ud.gene_families.at(i), j);

            full[j] = std::log(partial_likelihood[j]) + std::log(eq_freq);
        }

        //        all_families_likelihood[i] = accumulate(full.begin(), full.end(), 0.0); // sum over all sizes (Felsenstein's approach)
        all_families_likelihood[i] = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
                                                                             // cout << i << " contribution " << scientific << all_families_likelihood[i] << endl;
        
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

inference_optimizer_scorer *base_model::get_lambda_optimizer(const user_data& data)
{
    if (data.p_lambda != NULL)  // already have a lambda, nothing we want to optimize
        return nullptr;

    initialize_lambda(data.p_lambda_tree);

    if (_p_error_model && !data.p_error_model)
    {
        return new lambda_epsilon_optimizer(this, _p_error_model, data, _p_lambda);
    }
    else
    {
        return new sigma_optimizer_scorer(_p_lambda, this, data);
    }
}

#define EPSILON_RANGES

reconstruction* base_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    LOG(INFO) << "Starting reconstruction processes for Base model";

    auto result = new base_model_reconstruction();

    p_calc->precalculate_matrices(_p_lambda->get_lambdas(), get_all_bounds(ud.gene_families), ud.p_tree->get_branch_lengths());

    pupko_reconstructor::pupko_data data(ud.gene_families.size(), ud.p_tree, ud.max_family_size, ud.max_root_family_size);

    for (size_t i = 0; i < ud.gene_families.size(); ++i)
    {
        clademap<int> &rc = result->_reconstructions[ud.gene_families[i].id()];
        ud.p_tree->apply_prefix_order([&rc](const clade* c) {
            rc[c] = 0;
            });
    }

#pragma omp parallel for
    for (int i = 0; i< ud.gene_families.size(); ++i)
    {
        pupko_reconstructor::reconstruct_gene_transcript(_p_lambda, ud.p_tree, &ud.gene_families[i], p_calc, result->_reconstructions[ud.gene_families[i].id()], data.C(i), data.L(i));
    }

    size_t success = count_if(data.v_all_node_Ls.begin(), data.v_all_node_Ls.end(), [this, &ud](const clademap<std::vector<double>>& L)
        { 
            return *max_element( L.at(ud.p_tree).begin(), L.at(ud.p_tree).end()) > 0;
        });

    if (success != ud.gene_families.size())
    {
        LOG(WARNING) << "Failed to reconstruct " << ud.gene_families.size() - success << " families" << endl;
    }

    LOG(INFO) << "Done!\n";

    return result;
}

sigma* base_model::get_simulation_lambda()
{
    return _p_lambda->multiply(simulation_lambda_multiplier);
}

double base_model_reconstruction::get_node_count(const gene_transcript& family, const clade *c) const
{
    if (c->is_leaf())
        return family.get_expression_value(c->get_taxon_name());

    if (_reconstructions.find(family.id()) == _reconstructions.end())
        throw runtime_error("Family " + family.id() + " not found in reconstruction");

    return _reconstructions.at(family.id()).at(c);
}
