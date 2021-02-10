#include <iomanip>
#include <iostream>
#include <random>

#include "easylogging++.h"

#include "core.h"
#include "optimizer_scorer.h"
#include "clade.h"
#include "lambda.h"
#include "gamma.h"
#include "error_model.h"
#include "root_equilibrium_distribution.h"

#define GAMMA_INITIAL_GUESS_EXPONENTIAL_DISTRIBUTION_LAMBDA 1.75

extern std::mt19937 randomizer_engine;

using namespace std;

double inference_optimizer_scorer::calculate_score(const double *values)
{
    prepare_calculation(values);

    if (!quiet)
    {
        report_precalculation();
    }

    double score = _p_model->infer_family_likelihoods(*_p_distribution, _p_lambda);

    if (std::isnan(score)) score = -log(0);

    return score;
}

