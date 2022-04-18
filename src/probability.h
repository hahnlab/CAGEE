#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <Eigen/Dense>

class DiffMat;
class gene_transcript;
class sigma_squared;
class clade;
class error_model;
class matrix_cache;

struct pvalue_parameters
{
    const clade* p_tree;
    const sigma_squared* p_lambda;
    const int max_family_size;
    const int max_root_family_size;
    const matrix_cache& cache;
};

double chooseln(double n, double k);

std::vector<int> uniform_dist(int n_draws, int min, int max);

size_t adjust_for_error_model(size_t c, const error_model *p_error_model);

double pvalue(double v, const std::vector<double>& conddist);

std::vector<double> inference_prune(const gene_transcript& gf, 
    const matrix_cache& cache, 
    const sigma_squared* sigma, 
    const error_model* p_error_model, 
    const clade* _p_tree, 
    double _sigma_multiplier,
    int upper_bound);

class upper_bound_calculator
{
public:
    virtual int get(const gene_transcript& gt) const = 0;

    static upper_bound_calculator* create();
};

 
#endif
