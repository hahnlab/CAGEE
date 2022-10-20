#include "core.h"

//! @brief Represents a model of species change in which lambda values are expected to belong to a gamma distribution
//! \ingroup gamma
class gamma_model : public model {
private:
    //! Gamma
    std::vector<double> _sigma_multipliers;

    std::vector<double> _gamma_cat_probs; // each item is the probability of belonging to a given gamma category

    std::vector<std::vector<double>> _category_likelihoods;

    double _alpha;

    std::vector<double> get_posterior_probabilities(std::vector<double> cat_likelihoods);
public:

    //! Calculate gamma categories and lambda multipliers based on category count and a fixed alpha
    gamma_model(sigma_squared* p_sigma, std::vector<gene_transcript>* p_gene_transcripts, int n_gamma_cats, double fixed_alpha, error_model *p_error_model);

    //! Specify gamma categories and lambda multipliers explicitly
    gamma_model(sigma_squared* p_sigma, std::vector<gene_transcript>* p_gene_transcripts, std::vector<double> gamma_categories, std::vector<double> multipliers, error_model *p_error_model);

    void set_alpha(double alpha);
    double get_alpha() const { return _alpha; }

    void write_probabilities(std::ostream& ost);

    //! Randomly select one of the multipliers to apply to the simulation
    virtual sigma_squared* get_simulation_sigma() override;

    double infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma, const std::gamma_distribution<double>& prior) override;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups, const std::gamma_distribution<double>& prior) override;

    virtual void write_extra_vital_statistics(std::ostream& ost);

    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *) override;

    std::size_t get_gamma_cat_probs_count() const {
        return _gamma_cat_probs.size();
    }

    std::size_t get_sigma_multiplier_count() const {
        return _sigma_multipliers.size();
    }

    std::vector<double> get_sigma_multipliers() const {
        return _sigma_multipliers;
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

