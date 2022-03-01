class clade;
class gene_transcript;
class user_data;
class sigma_squared;
class optimizer_scorer;
class model;

#include <vector>

namespace LikelihoodRatioTest {
    clade* update_branchlength(const clade* lambda_tree, double bl_augment, int t);

    double get_likelihood_for_diff_lambdas(const gene_transcript& gf, const clade* p_tree,
        const clade* p_lambda_tree, int t, std::vector<sigma_squared*>& lambda_cache, optimizer* p_opt);

    void compute_for_diff_lambdas_i(const user_data& data,
        std::vector<int>& lambda_index,
        std::vector<double>& pvalues,
        std::vector<sigma_squared*>& lambda_cache,
        optimizer* p_opt
    );

    void lhr_for_diff_lambdas(const user_data& data, model *p_model);
}

