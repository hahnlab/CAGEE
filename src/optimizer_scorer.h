#ifndef OPTIMIZER_SCORER_H
#define OPTIMIZER_SCORER_H

#include <vector>
#include <map>

class error_model;
class lambda;
class model;
class root_equilibrium_distribution;
class clade;
class base_model;

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

    lambda *_p_lambda;
    model *_p_model;
    const root_equilibrium_distribution *_p_distribution;

public:
    inference_optimizer_scorer(lambda *p_lambda, model* p_model, const root_equilibrium_distribution *p_distribution) :
        _p_lambda(p_lambda),
        _p_model(p_model),
        _p_distribution(p_distribution),
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






#endif
