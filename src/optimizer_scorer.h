#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>
#include <map>
#include <random>

class user_data;
class model;
class gamma_model;
class sigma_squared;
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

//! @brief Scorer that optimizes for sigma
//! \ingroup optimizer
class sigma_optimizer_scorer : public optimizer_scorer
{
    error_model* _p_error_model = nullptr;
    std::vector<double> current_epsilon_guesses;
    bool optimize_sigma, optimize_epsilon, optimize_gamma;
    sigma_squared* _p_sigma;
    model* _p_model;
    const user_data& _user_data;
    const std::gamma_distribution<double>& _prior;

    double* _distribution_mean = nullptr;
public:
    // sigma only
    sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma_squared* p_lambda);

    // sigma and epsilon
    sigma_optimizer_scorer(model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma_squared* p_lambda, error_model* p_error_model);

    // alpha and sigma
    sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior, sigma_squared* p_lambda);

    // alpha only
    sigma_optimizer_scorer(gamma_model* p_model, const user_data& user_data, const std::gamma_distribution<double>& prior);

    void force_distribution_mean(double tl, double v);

    std::string description() const;
    std::vector<double> initial_guesses() override;
    double calculate_score(const double* values) override;

    virtual void finalize(double *results);

    virtual void prepare_calculation(const double *values);
    virtual void report_precalculation();

#ifdef SILENT
    bool quiet = true;
#else
    bool quiet = false;
#endif
};

#endif
