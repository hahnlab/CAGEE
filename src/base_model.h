#ifndef BASE_MODEL_H
#define BASE_MODEL_H

#include "core.h"

class gene_family_reconstructor;
class matrix_cache;
class gene_transcript;
class lambda;
class error_model;
class root_equilibrium_distribution;
class reconstruction;

/*! @brief A Base model can simulate families or estimate lambdas and error models.

    The estimation method creates a best guess for the missing lambda and, optionally,
    epsilon values, by using an @optimizer to provide guesses for the missing values
    and calculating a score for each one.

    \defgroup base Base Model

 */
class base_model : public model {
    double simulation_lambda_multiplier = 1.0;

public:
    //! Computation or estimation constructor
    base_model(lambda* p_lambda, const clade *p_tree, const std::vector<gene_transcript>* p_gene_families,
        int max_family_size, int max_root_family_size, error_model *p_error_model);

    virtual double infer_family_likelihoods(const root_equilibrium_distribution& prior, const lambda *p_lambda) override;

    virtual std::string name() const {
        return "Base";
    }

    virtual void write_family_likelihoods(std::ostream& ost);

    virtual inference_optimizer_scorer *get_lambda_optimizer(const user_data& data);

    virtual reconstruction* reconstruct_ancestral_states(const std::vector<gene_transcript>& families, matrix_cache *p_calc, root_equilibrium_distribution* p_prior);

    virtual void prepare_matrices_for_simulation(matrix_cache& cache);

    virtual lambda* get_simulation_lambda();

};

//! \ingroup base Base Model
class base_model_reconstruction : public reconstruction
{
public:

    base_model_reconstruction()
    {

    }

    std::map<std::string, clademap<int>> _reconstructions;

    double get_node_count(const gene_transcript& gf, const clade* c) const override;

};

#endif
