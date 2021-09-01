#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <Eigen/Dense>

class DiffMat;
class gene_transcript;
class sigma;
class clade;
class error_model;
class matrix_cache;

struct pvalue_parameters
{
    const clade* p_tree;
    const sigma* p_lambda;
    const int max_family_size;
    const int max_root_family_size;
    const matrix_cache& cache;
};

double chooseln(double n, double k);

std::vector<int> uniform_dist(int n_draws, int min, int max);

std::vector<double> get_random_probabilities(pvalue_parameters p, int number_of_simulations, int root_family_size);
size_t adjust_for_error_model(size_t c, const error_model *p_error_model);

double pvalue(double v, const std::vector<double>& conddist);

//! computes a pvalue for each family. Returns a vector of pvalues matching the list of families
std::vector<double> compute_pvalues(pvalue_parameters p, const std::vector<gene_transcript>& families, int number_of_simulations);

std::vector<double> inference_prune(const gene_transcript& gf, const matrix_cache& cache, const sigma* _lambda, const error_model* p_error_model, const clade* _p_tree, double _lambda_multiplier);

std::pair<double, double> bounds(const gene_transcript& gt);

#endif
