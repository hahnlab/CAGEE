#ifndef GAMMA_CORE_H
#define GAMMA_CORE_H

#include <memory>

#include "core.h"
#include "discretized_gamma.h"

class prior;
class likelihood_combiner;

//! @brief Represents a model of species change in which lambda values are expected to belong to a gamma distribution
//! \ingroup gamma
class gamma_model : public model {
private:
    discretized_gamma _gamma;
    const int _n_gamma_cats;
    std::map<const gene_transcript*, std::vector<double>> _transcript_category_likelihoods;
public:

    //! Calculate gamma categories and lambda multipliers based on category count and a fixed alpha
    gamma_model(sigma_squared* p_sigma, std::vector<gene_transcript>* p_gene_transcripts, int n_gamma_cats, double fixed_alpha, error_model *p_error_model);
    ~gamma_model();

    void set_alpha(double alpha);
    double get_alpha() const { return _gamma.get_alpha(); }

    void write_probabilities(std::ostream& ost);

    //! Randomly select one of the multipliers to apply to the simulation
    virtual sigma_squared* get_simulation_sigma() override;

    double infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma) override;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups) override;

    virtual void write_extra_vital_statistics(std::ostream& ost);

    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *) override;

    bool can_infer() const;

    virtual std::string get_name() const override { return "Gamma"; }
};

#endif // GAMMA_CORE_H
