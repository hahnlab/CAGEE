#ifndef BOUNDED_BROWNIAN_MOTION_MODEL_H
#define BOUNDED_BROWNIAN_MOTION_MODEL_H

#include "core.h"
#include "optimizer_scorer.h"

//! @brief Scorer that optimizes for sigma
//! \ingroup optimizer
class sigma_optimizer_scorer : public optimizer_scorer
{
    const clade* _p_tree;
    const vector<gene_family>& _families;

public:
    sigma_optimizer_scorer(const clade* p_tree, const vector<gene_family>& families) : _p_tree(p_tree), _families(families)
    {
    }

    virtual std::vector<double> initial_guesses() override;

    virtual double calculate_score(const double* values) override;
};

class bounded_brownian_motion_model : public model
{
public:
    bounded_brownian_motion_model(lambda* p_lambda,
        const clade* p_tree,
        const std::vector<gene_family>* p_gene_families,
        int max_family_size,
        int max_root_family_size,
        error_model* p_error_model) : model(p_lambda, p_tree, p_gene_families, max_family_size, max_root_family_size, p_error_model) {}
private:
    virtual void prepare_matrices_for_simulation(matrix_cache& cache) {}

    virtual double infer_family_likelihoods(const root_equilibrium_distribution& prior, const lambda* p_lambda) { return 0.0; }

    virtual std::string name() const { return "Bounded Brownian Motion"; }
    virtual void write_family_likelihoods(std::ostream& ost) {}

    //! Based on the model parameters, attempts to reconstruct the most likely counts of each family at each node
    virtual reconstruction* reconstruct_ancestral_states(const vector<gene_family>& families, matrix_cache* p_calc, root_equilibrium_distribution* p_prior) { return nullptr; }

    virtual optimizer_scorer* get_sigma_optimizer(const user_data& data);
};
#endif


