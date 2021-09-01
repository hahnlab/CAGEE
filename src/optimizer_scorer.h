#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>
#include <map>

class user_data;
class model;
class sigma;
class error_model;

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

public:
    inference_optimizer_scorer(sigma*p_sigma, model* p_model, const user_data& user_data) :
        _p_sigma(p_sigma),
        _p_model(p_model),
        _user_data(user_data),
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

public:
    sigma_optimizer_scorer(sigma*p_lambda, model* p_model, const user_data& user_data, double tree_length, double species_variance) :
        inference_optimizer_scorer(p_lambda, p_model, user_data),
        _tree_length(tree_length),
        _species_variance(species_variance)
    {
    }

    sigma_optimizer_scorer(sigma* p_lambda, model* p_model, const user_data& user_data);
    
    std::vector<double> initial_guesses() override;

    virtual void finalize(double *results) override;

    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;
};


/// @brief Scorer that optimizes lambdas and epsilons together
//! \ingroup optimizer
class lambda_epsilon_optimizer : public inference_optimizer_scorer
{
    sigma_optimizer_scorer _lambda_optimizer;
    error_model* _p_error_model;
    std::vector<double> current_guesses;
public:
    lambda_epsilon_optimizer(
        model* p_model,
        error_model *p_error_model,
        const user_data& user_data,
        sigma*p_lambda,
        double tree_length, 
        double species_variance) :
        inference_optimizer_scorer(p_lambda, p_model, user_data),
        _lambda_optimizer(p_lambda, p_model, user_data, tree_length, species_variance),
        _p_error_model(p_error_model)
    {
    }

    lambda_epsilon_optimizer(
        model* p_model,
        error_model* p_error_model,
        const user_data& user_data,
        sigma* p_lambda) :
        inference_optimizer_scorer(p_lambda, p_model, user_data),
        _lambda_optimizer(p_lambda, p_model, user_data),
        _p_error_model(p_error_model)
    {
    }

    std::vector<double> initial_guesses() override;

    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;

    virtual void finalize(double *results) override;
};

class gamma_model;

//! @brief Scorer that optimizes for alpha
//! \ingroup optimizer
//! \ingroup gamma
class gamma_optimizer : public inference_optimizer_scorer {
    gamma_model *_p_gamma_model;
public:
    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual void finalize(double * result) override;
    gamma_optimizer(gamma_model* p_model, const user_data& user_data);

    double get_alpha() const;
};

//! @brief Scorer that optimizes for both lambda and alpha
//! \ingroup optimizer
//! \ingroup gamma
class gamma_lambda_optimizer : public inference_optimizer_scorer
{
    virtual void prepare_calculation(const double *values) override;
    virtual void report_precalculation() override;
    sigma_optimizer_scorer _lambda_optimizer;
    gamma_optimizer _gamma_optimizer;
public:
    gamma_lambda_optimizer(sigma*p_lambda, gamma_model * p_model, const user_data& user_data, double tree_length, double species_variance);
    gamma_lambda_optimizer(sigma* p_lambda, gamma_model* p_model, const user_data& user_data);

    std::vector<double> initial_guesses() override;

    /// results consists of the desired number of lambdas and one alpha value
    void finalize(double *results) override;
};




#endif
