#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>
#include <map>
#include <random>

class user_data;
class model;
class gamma_model;
class sigma;
class error_model;
class root_distribution_gamma;

/// @brief Base class for use by the optimizer
//! \ingroup optimizer
class optimizer_scorer {
public:
    virtual ~optimizer_scorer() {}

    virtual std::vector<double> initial_guesses() = 0;

    virtual double calculate_score(const double *values) = 0;
};

//! @brief  Scorer that holds a model and calls its inference method
//! for scoring
//! \ingroup optimizer
class inference_optimizer_scorer : public optimizer_scorer {
protected:
    virtual void prepare_calculation(const double *values) = 0;
    virtual void report_precalculation() = 0;

    sigma*_p_sigma;
    model *_p_model;
    const user_data& _user_data;
    const std::gamma_distribution<double>&  _prior;

public:
    inference_optimizer_scorer(sigma*p_sigma, model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior) :
        _p_sigma(p_sigma),
        _p_model(p_model),
        _user_data(user_data),
        _prior(prior),
        quiet(false)
    {
#ifdef SILENT
        quiet = true;
#endif
    }

    virtual ~inference_optimizer_scorer() {}

    double calculate_score(const double *values) ;

    virtual void finalize(double *result) = 0;

    bool quiet;
};

//! @brief Scorer that optimizes for sigma
//! \ingroup optimizer
class sigma_optimizer_scorer : public inference_optimizer_scorer
{
    double _tree_length;
    double _species_variance;
    error_model* _p_error_model = nullptr;
    std::vector<double> current_epsilon_guesses;
    bool optimize_sigma, optimize_epsilon, optimize_gamma;

public:
    // sigma only
    sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma* p_lambda);

    // sigma and epsilon
    sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma* p_lambda, error_model* p_error_model);

    // alpha and sigma
    sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma* p_lambda);

    // alpha only
    sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior);

    void force_tree_length_and_variance(double tl, double v) {
        _tree_length = tl;
        _species_variance = v;
    }

    std::string description() const;
    std::vector<double> initial_guesses() override;

    virtual void finalize(double *results) override;

    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;
};

#endif
