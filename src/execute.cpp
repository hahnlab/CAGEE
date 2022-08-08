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
#include "transcript_reconstructor.h"
#include "error_model.h"
#include "likelihood_ratio.h"
#include "arguments.h"
#include "root_equilibrium_distribution.h"
#include "sigma.h"
#include "DiffMat.h"

using namespace std;

double __Qs[] = { 1.000000000190015, 76.18009172947146, -86.50532032941677,
24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
-5.395239384953e-6 };

estimator::estimator(user_data& d, const input_parameters& ui) : action(d, ui)
{
    auto tokens = tokenize_str(ui.prior_params_or_default(), ':');
    if (tokens[0] != "gamma")
        throw std::runtime_error("Prior must be given in the form gamma:k:theta");

    auto k = stof(tokens[1]), theta = stof(tokens[2]);
    LOG(INFO) << "Using gamma prior with k=" << k << ", theta=" << theta << ")";
    _prior = gamma_distribution<double>(k, theta);

#ifdef USE_MAX_PROBABILITY
    LOG(INFO) << "Maximum probability calculation";
#else
    LOG(INFO) << "Summation probability calculation";
#endif

}

void estimator::write_error_model_if_specified(const input_parameters& my_input_parameters, const model * p_model)
{
    if (my_input_parameters.use_error_model)
    {
        ofstream errmodel(filename(p_model->name() + "_error_model", _user_input.output_prefix));
        if (data.p_error_model)
        {
            /// user specified an error model, write that out to the results directory
            write_error_model_file(errmodel, *data.p_error_model);
        }
        else
        {
            /// user did not specify an error model, write out the estimated one or a default
            p_model->write_error_model(200, errmodel);
        }
    }
}

void check_tree(const clade* p_tree, gene_transcript& t)
{
    p_tree->apply_prefix_order([p_tree, t](const clade* c) {
        if (c->is_leaf())
        {
            t.get_expression_value(c->get_taxon_name());
        }
        });
}

void estimator::compute(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    std::vector<double> model_likelihoods(models.size());
    for (size_t i = 0; i < models.size(); ++i) {
        LOG(INFO) << "Inferring processes for " << models[i]->name() << " model";

        double result = models[i]->infer_family_likelihoods(data, models[i]->get_sigma(), _prior);
        std::ofstream results_file(filename(models[i]->name() + "_results", my_input_parameters.output_prefix));
        models[i]->write_vital_statistics(results_file, data.p_tree, result);

        std::ofstream likelihoods_file(filename(models[i]->name() + "_family_likelihoods", my_input_parameters.output_prefix));
        models[i]->write_family_likelihoods(likelihoods_file);

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
        unique_ptr<sigma_optimizer_scorer> scorer(p_model->get_sigma_optimizer(data, _user_input.sample_groups, _prior));
        if (scorer.get() == nullptr)
            continue;   // nothing to be optimized

        optimizer opt(scorer.get());

        auto result = opt.optimize(_user_input.optimizer_params);
        scorer->finalize(&result.values[0]);

        LOG(INFO) << p_model->get_monitor();
    }
    if (data.p_lambda == nullptr)
        data.p_lambda = models[0]->get_sigma();

}

void estimator::estimate_lambda_per_family(model *p_model, ostream& ost)
{
    auto families = data.gene_families;
    vector<sigma_squared*> result(data.gene_families.size());
    std::transform(data.gene_families.begin(), data.gene_families.end(), result.begin(),
        [this, p_model](gene_transcript& fam)
    {
#ifndef SILENT
            LOG(INFO) << "Estimating for " << fam.id() << endl;
#endif
        vector<gene_transcript> v({ fam });
        vector<model *> models{ p_model };
        throw std::runtime_error("Not implemented yet");
        //p_model->set_families(&v);
        data.p_lambda = nullptr;
        estimate_missing_variables(models, data);
        return p_model->get_sigma();
    });
    std::transform(data.gene_families.begin(), data.gene_families.end(), result.begin(),
        ostream_iterator<string>(ost, "\n"),
        [](const gene_transcript& fam, sigma_squared* lambda)
    {
        return fam.id() + '\t' + lambda->to_string();
    });
    for (auto r : result) delete r; // TODO: use unique_ptrs in result for exception safety

}

/*! Calls estimate_lambda_per_family if the user has set that parameter, otherwise
    calls \ref estimate_missing_variables; \ref compute; \ref compute_pvalues, and \ref model::reconstruct_ancestral_states */
void estimator::execute(std::vector<model *>& models)
{
    check_tree(data.p_tree, data.gene_families[0]);

    string dir = _user_input.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);
    
    if (_user_input.lambda_per_family)
    {
        auto p_model = models[0];   // no support for multiple models
        std::ofstream results_file(filename(p_model->name() + "_lambda_per_family", _user_input.output_prefix));
        estimate_lambda_per_family(p_model, results_file);
    }
    else
    {
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

#ifdef RUN_LHRTEST
                LikelihoodRatioTest::lhr_for_diff_lambdas(data, p_model);
#endif
                rec->write_results(p_model->name(), _user_input.output_prefix, data.p_tree, data.gene_families, _user_input.pvalue);
            }
        }
        catch (const OptimizerInitializationFailure& e )
        {
            initialization_failure_advice(cerr, data.gene_families);
            throw;
        }
    }
}

//Calculate the difference between the Max and Min count for each family, report the 20 families with the largest difference.
void initialization_failure_advice(std::ostream& ost, const std::vector<gene_transcript>& families)
{
    std::vector<std::pair<std::string, double>> m;
    transform(families.begin(), families.end(), std::inserter(m, m.end()),
        [](const gene_transcript& gf) { return std::make_pair(gf.id(), gf.species_size_differential()); });
    auto compare = [](const std::pair<string, double>& a, const std::pair<string, double>& b) { return a.second > b.second; };
    sort(m.begin(), m.end(), compare);
    if (m.size() > 20)
        m.resize(20);

    ost << "\nFamilies with largest size differentials:\n";
    for (auto& t : m)
        ost << t.first << ": " << t.second << "\n";
    ost << "\nYou may want to try removing the top few families with the largest difference\nbetween the max and min counts and then re-run the analysis.\n\n";
}


TEST_CASE("create_rootdist throws error if nothing set")
{
    user_data ud;
    input_parameters ip;
    estimator e(ud, ip);
    CHECK_EQ(gamma_distribution<double>(0.375, 1600.0), e.prior());
}

TEST_CASE("check_tree checks tree against a gene family to make sure all leaf node species are available")
{
    gene_transcript gt;
    gt.set_expression_value("Dog", 0.254007);
    gt.set_expression_value("Rat", 0.1);

    unique_ptr<clade> p_tree(parse_newick("(Dog:1,Rat:3):7"));

    check_tree(p_tree.get(), gt);
}

TEST_CASE("check_tree throws error on unknown species")
{
    gene_transcript gt;
    gt.set_expression_value("Dog", 0.254007);
    gt.set_expression_value("Cow", 0.1);

    unique_ptr<clade> p_tree(parse_newick("(Dog:1,Rat:3):7"));

    CHECK_THROWS_WITH(check_tree(p_tree.get(), gt), "Rat was not found in gene family ");
}