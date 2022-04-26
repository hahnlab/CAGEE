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

template<typename T>
using clademap = std::map<const clade*, T>;

class inference_pruner
{
    const matrix_cache& _cache;
    const sigma_squared* _p_sigsqd;
    const error_model* _p_error_model;
    const clade* _p_tree;
    const double _sigma_multiplier;

    clademap<Eigen::VectorXd> compute_all_probabilities(const gene_transcript& gf, int upper_bound);
public:
    inference_pruner(const matrix_cache& cache,
        const sigma_squared* sigma,
        const error_model* p_error_model,
        const clade* _p_tree,
        double _sigma_multiplier);

    std::vector<double> prune(const gene_transcript& gf, int upper_bound);
    clademap<double> reconstruct(const gene_transcript& gf, int upper_bound);
};

class upper_bound_calculator
{
public:
    virtual int get(const gene_transcript& gt) const = 0;

    static upper_bound_calculator* create();
};

 
#endif
