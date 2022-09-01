#include "core.h"

//! @brief Represents a model of species change in which lambda values are expected to belong to a gamma distribution
//! \ingroup gamma
class gamma_model : public model {
private:
    //! Gamma
    std::vector<double> _lambda_multipliers;

    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    std::vector<std::vector<double>> _category_likelihoods;

    double _alpha;

    std::vector<double> get_posterior_probabilities(std::vector<double> cat_likelihoods);
public:

    //! Calculate gamma categories and lambda multipliers based on category count and a fixed alpha
    gamma_model(sigma_squared* p_lambda, std::vector<gene_transcript>* p_gene_families, int n_gamma_cats, double fixed_alpha, error_model *p_error_model);

    //! Specify gamma categories and lambda multipliers explicitly
    gamma_model(sigma_squared* p_lambda, std::vector<gene_transcript>* p_gene_families, std::vector<double> gamma_categories, std::vector<double> multipliers, error_model *p_error_model);

    void set_alpha(double alpha);
    double get_alpha() const { return _alpha; }

    void write_probabilities(std::ostream& ost);

    //! Randomly select one of the multipliers to apply to the simulation
    virtual sigma_squared* get_simulation_lambda() override;

    double infer_family_likelihoods(const user_data& ud, const sigma_squared*p_lambda, const std::gamma_distribution<double>& prior) override;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups, const std::gamma_distribution<double>& prior) override;

    virtual std::string name() const override {
        return "Gamma";
    }

    virtual void write_vital_statistics(std::ostream& ost, const clade *p_tree, double final_likelihood) override;

    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *) override;

    std::size_t get_gamma_cat_probs_count() const {
        return _gamma_cat_probs.size();
    }

    std::size_t get_lambda_multiplier_count() const {
        return _lambda_multipliers.size();
    }

    std::vector<double> get_lambda_multipliers() const {
        return _lambda_multipliers;
    }

    bool can_infer() const;

    bool prune(const gene_transcript& family, 
        const std::gamma_distribution<double>& prior, 
        const matrix_cache& calc, 
        const sigma_squared*p_sigma,
        const clade *p_tree, 
        std::vector<double>& category_likelihoods, 
        int upper_bound);
};

