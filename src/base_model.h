#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

/*! @brief A Base model can simulate families or estimate sigma and error models.

    The estimation method creates a best guess for the missing sigma and, optionally,
    epsilon values, by using an @optimizer to provide guesses for the missing values
    and calculating a score for each one.

    \defgroup base Base Model

 */
class base_model : public model {
    double simulation_sigma_multiplier = 1.0;

public:
    //! Computation or estimation constructor
    base_model(sigma_squared* p_sigma, const std::vector<gene_transcript>* p_gene_transcripts,
        error_model *p_error_model);

    virtual double infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma, const std::gamma_distribution<double>& prior) override;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups, const std::gamma_distribution<double>& prior) override;

    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc);

    virtual sigma_squared* get_simulation_sigma();

};



#endif
