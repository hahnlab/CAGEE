#ifndef ROOTDIST_ESTIMATOR_H
#define ROOTDIST_ESTIMATOR_H

#include <vector>
#include "optimizer_scorer.h"

class gene_transcript;

class poisson_scorer : public optimizer_scorer
{
    std::vector<double> leaf_family_sizes;
public:
    poisson_scorer(const std::vector<gene_transcript>& gene_families);

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual double calculate_score(const double * values) override;

    double lnLPoisson(const double* plambda);
};

class gamma_scorer : public optimizer_scorer
{
    std::vector<double> leaf_family_sizes;
public:
    gamma_scorer(const std::vector<gene_transcript>& gene_families);

    // Inherited via optimizer_scorer
    virtual std::vector<double> initial_guesses() override;
    virtual double calculate_score(const double* values) override;
};


double poisspdf(double x, double lambda);
double gammapdf(double value, const std::gamma_distribution<double>& dist);

std::vector<double> get_prior_rfsize_poisson_lambda(int min_family_size, int max_family_size, double poisson_lambda);

#endif
