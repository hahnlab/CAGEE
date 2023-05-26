#include <cmath>
#include <set>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "doctest.h"
#include "easylogging++.h"

#include "execute.h"
#include "core.h"
#include "user_data.h"
#include "chisquare.h"
#include "optimizer_scorer.h"
#include "matrix_cache.h"
#include "optimizer.h"
#include "reconstruction.h"
#include "error_model.h"
#include "arguments.h"
#include "root_equilibrium_distribution.h"
#include "sigma.h"
#include "DiffMat.h"
#include "replicate_model.h"

using namespace std;

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };


estimator::estimator(user_data& d, const input_parameters& ui) : action(d, ui)
{
#ifdef USE_MAX_PROBABILITY
    LOG(INFO) << "Maximum probability calculation";
#else
    LOG(INFO) << "Summation probability calculation";
#endif

}

// void estimator::write_error_model_if_specified(const input_parameters& my_input_parameters, const model * p_model)
// {
//     if (my_input_parameters.use_error_model)
//     {
//         ofstream errmodel(filename("error_model", _user_input.output_prefix));
//         if (data.p_error_model)
//         {
//             /// user specified an error model, write that out to the results directory
//             write_error_model_file(errmodel, *data.p_error_model);
//         }
//         else
//         {
//             /// user did not specify an error model, write out the estimated one or a default
//             p_model->write_error_model(200, errmodel);
//         }
//     }
// }

void check_tree(const user_data& data)
{
    if (data.p_replicate_model)
    {
        data.p_replicate_model->verify_replicates(data.p_tree, data.gene_transcripts[0]);
    }
    else
    {
        data.p_tree->apply_prefix_order([data](const clade* c) {
            if (c->is_leaf())
            {
                data.gene_transcripts[0].get_expression_value(c->get_taxon_name());
            }
        });
    }

}

void estimator::compute(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    std::vector<double> model_likelihoods(models.size());
    for (size_t i = 0; i < models.size(); ++i) {
        LOG(INFO) << "Inferring processes";

        double result = models[i]->infer_transcript_likelihoods(data, models[i]->get_sigma());
        std::ofstream results_file(filename("results", my_input_parameters.output_prefix));
        models[i]->write_vital_statistics(results_file, data.p_tree, result, my_input_parameters);

        write_error_model_if_specified(my_input_parameters, models[i]);

        model_likelihoods[i] = result;
    }

    if (model_likelihoods.size() == 2)
    {
        LOG(INFO) << "PValue = " << (1.0 - chi2cdf(2 * (model_likelihoods[1] - model_likelihoods[0]), 1.0));
    }
}

bool compare_result(const optimizer::result& a, const optimizer::result& b)
{
    return a.score < b.score;
}

void estimator::estimate_missing_variables(std::vector<model *>& models, user_data& data)
{
    if (data.p_tree == NULL)
    {
        throw runtime_error("No tree specified for lambda estimation");
    }
    for (model* p_model : models) {
        unique_ptr<sigma_optimizer_scorer> scorer(p_model->get_sigma_optimizer(data, _user_input.sample_groups));
        if (scorer.get() == nullptr)
            continue;   // nothing to be optimized

        optimizer opt(scorer.get());

        auto result = opt.optimize(_user_input.optimizer_params);
        scorer->finalize(&result.values[0]);

        LOG(INFO) << p_model->get_monitor();
    }
    if (data.p_sigma == nullptr)
        data.p_sigma = models[0]->get_sigma();

}

void estimator::execute(std::vector<model *>& models)
{
    check_tree(data);

    string dir = _user_input.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);
    
    try
    {
        estimate_missing_variables(models, data);

        compute(models, _user_input);

        for (model* p_model : models) {

            /// For Gamma models, we tried using the most rapidly changing lambda multiplier here, but that
            /// caused issues in the pvalue calculation. It should be best to use the original lambda
            /// instead
            matrix_cache cache;

            std::unique_ptr<reconstruction> rec(p_model->reconstruct_ancestral_states(data, &cache));

            rec->write_results(_user_input.output_prefix, data.p_tree, data.gene_transcripts, _user_input.count_all_changes);
        }
    }
    catch (const OptimizerInitializationFailure& e )
    {
        initialization_failure_advice(cerr, data.gene_transcripts);
        throw;
    }
}

void initialization_failure_advice(std::ostream& ost, const std::vector<gene_transcript>& gt)
{
    std::vector<std::pair<std::string, double>> m;
    transform(gt.begin(), gt.end(), std::inserter(m, m.end()),
        [](const gene_transcript& gf) { return std::make_pair(gf.id(), gf.species_size_differential()); });
    auto compare = [](const std::pair<string, double>& a, const std::pair<string, double>& b) { return a.second > b.second; };
    sort(m.begin(), m.end(), compare);
    if (m.size() > 20)
        m.resize(20);

    ost << "\nTranscripts with largest size differentials:\n";
    for (auto& t : m)
        ost << t.first << ": " << t.second << "\n";
    ost << "\nYou may want to try removing the top few transcripts with the largest difference\nbetween the max and min counts and then re-run the analysis.\n\n";
}

TEST_CASE("check_tree checks tree against a gene family to make sure all leaf node species are available")
{
    user_data ud;
    ud.gene_transcripts.resize(1);
    gene_transcript gt;
    ud.gene_transcripts[0].set_expression_value("Dog", 0.254007);
    ud.gene_transcripts[0].set_expression_value("Rat", 0.1);

    unique_ptr<clade> p_tree(parse_newick("(Dog:1,Rat:3):7"));
    ud.p_tree = p_tree.get();
    check_tree(ud);
}

TEST_CASE("check_tree throws error on unknown species")
{
    user_data ud;
    ud.gene_transcripts.resize(1);
    gene_transcript gt;
    ud.gene_transcripts[0].set_expression_value("Dog", 0.254007);
    ud.gene_transcripts[0].set_expression_value("Cow", 0.1);

    unique_ptr<clade> p_tree(parse_newick("(Dog:1,Rat:3):7"));
    ud.p_tree = p_tree.get();

    CHECK_THROWS_WITH(check_tree(ud), "Rat was not found in transcript ");
}