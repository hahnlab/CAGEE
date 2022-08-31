#include <memory>
#include <cmath>
#include <algorithm>

#include "core.h"
#include "matrix_cache.h"
#include "user_data.h"
#include "chisquare.h"
#include "optimizer_scorer.h"
#include "root_equilibrium_distribution.h"
#include "optimizer.h"
#include "sigma.h"
#include "DiffMat.h"
#include "inference_pruner.h"

namespace LikelihoodRatioTest
{
    using namespace std;

    clade * update_branchlength(const clade * p_tree, double bl_augment, int t)
    {
        return new clade(*p_tree, nullptr, [&](const clade& c) {
            return c.get_branch_length() + (c.get_branch_length() + bl_augment * t);  });
    }

    double get_likelihood_for_diff_lambdas(const gene_transcript & gf, const clade * p_tree, const clade * p_lambda_tree, 
        int lambda_index, 
        std::vector<sigma_squared*> & lambda_cache,
        optimizer *opt,
        int upper_bound)
    {
        const double bl_augment = 0.5;
        unique_ptr<clade> adjusted_tree(update_branchlength(p_tree, bl_augment, lambda_index));
        if (lambda_cache[lambda_index] == nullptr)
        {
            auto result = opt->optimize(optimizer_parameters());
            if (p_lambda_tree)
                lambda_cache[lambda_index] = new sigma_squared(map<string, int>(), result.values, sigma_type::lineage_specific);
            else
                lambda_cache[lambda_index] = new sigma_squared(result.values[0]);
        }

        matrix_cache cache;
        inference_pruner pruner(cache, lambda_cache[lambda_index], nullptr, adjusted_tree.get(), 1.0);
        auto probs = pruner.prune(gf, upper_bound);
        return *max_element(probs.begin(), probs.end());
    }

    void compute_for_diff_lambdas_i(const user_data & data,
        std::vector<int> & lambda_index,
        std::vector<double> & pvalues,
        std::vector<sigma_squared*> & lambda_cache,
        optimizer* p_opt
    )
    {
        auto references = build_reference_list(data.gene_families);
        int upper_bound = upper_bound_from_transcript_values(data.gene_families);
        matrix_cache cache;
        inference_pruner pruner(cache, data.p_lambda, data.p_error_model, data.p_tree, 1.0);
        for (size_t i = 0; i < data.gene_families.size(); i += 1)
        {
            auto& pitem = data.gene_families[i];
            if (references[i] != i) continue;

            //cache.precalculate_matrices(get_lambda_values(data.p_lambda), data.p_tree->get_branch_lengths());
            auto values = pruner.prune(pitem, upper_bound);
            double maxlh1 = *max_element(values.begin(), values.end());
            double prev = -1;
            double next = get_likelihood_for_diff_lambdas(pitem, data.p_tree, data.p_lambda_tree, 0, lambda_cache, p_opt, upper_bound);
            int j = 1;
            for (; prev < next; j++)
            {
                prev = next;
                next = get_likelihood_for_diff_lambdas(pitem, data.p_tree, data.p_lambda_tree, j, lambda_cache, p_opt, upper_bound);
            }
            pvalues[i] = (prev == maxlh1) ? 1 : 2 * (log(prev) - log(maxlh1));
            lambda_index[i] = j - 2;
        }
    }

    void likelihood_ratio_report(std::ostream & ost, const std::vector<gene_transcript> & families,
        const clade * pcafe,
        const std::vector<double> & pvalues,
        const std::vector<int> & plambda,
        const std::vector<sigma_squared*> & lambda_cache)
    {
        for (size_t i = 0; i < families.size(); ++i)
        {
            ost << families[i].id() << "\t";
            pcafe->write_newick(ost, [](const clade* c) { return c->get_taxon_name(); });
            auto l = lambda_cache[plambda[i]];
            cout << "(" << plambda[i] << ", " << *l << ")\t" << pvalues[i] << "\t" << (pvalues[1] == 1 ? 1 : 1 - chi2cdf(pvalues[i], 1)) << endl;
        }
    }

    void lhr_for_diff_lambdas(const user_data & data, model *p_model, const std::gamma_distribution<double>& prior)
    {
        std::vector<sigma_squared*> lambda_cache(100);

        cout << "Running Likelihood Ratio Test 2....\n";

        std::vector<double> pvalues(data.gene_families.size());
        std::vector<int> lambdas(data.gene_families.size());

        auto scorer = new sigma_optimizer_scorer(p_model, data, prior, data.p_lambda);

        optimizer opt(scorer);
        opt.quiet = true;
        compute_for_diff_lambdas_i(data, lambdas, pvalues, lambda_cache, &opt);

        likelihood_ratio_report(cout, data.gene_families, data.p_tree, pvalues, lambdas, lambda_cache);
    }
}