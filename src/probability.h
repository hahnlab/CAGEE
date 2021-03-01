#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <tuple>
#include <set>

#include <iostream>
#include <set>
#include "io.h"
#include "lambda.h"
#include "clade.h"

class matrix;
class matrix_cache;
class gene_transcript;

struct pvalue_parameters
{
    const clade* p_tree;
    const lambda* p_lambda;
    const int max_family_size;
    const int max_root_family_size;
    const matrix_cache& cache;
};


double chooseln(double n, double k);

/* START: Uniform distribution */
std::vector<int> uniform_dist(int n_draws, int min, int max);
/* END: Uniform distribution - */

std::vector<double> get_random_probabilities(pvalue_parameters p, int number_of_simulations, int root_family_size);
size_t adjust_for_error_model(size_t c, const error_model *p_error_model);

double pvalue(double v, const vector<double>& conddist);

//! computes a pvalue for each family. Returns a vector of pvalues matching the list of families
std::vector<double> compute_pvalues(pvalue_parameters p, const std::vector<gene_transcript>& families, int number_of_simulations);

#endif
