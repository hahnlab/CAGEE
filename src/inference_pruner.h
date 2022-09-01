#ifndef PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939
#define PROBABILITY_H_A2E01F6E_6A7D_44FB_A9C0_6512F15FF939

#include <vector>
#include <Eigen/Dense>

#include "reconstruction.h"

class DiffMat;
class gene_transcript;
class sigma_squared;
class clade;
class error_model;
class matrix_cache;

size_t adjust_for_error_model(size_t c, const error_model *p_error_model);

template<typename T>
using clademap = std::map<const clade*, T>;

class inference_pruner
{
    const matrix_cache& _cache;
    const sigma_squared* _p_sigsqd;
    const error_model* _p_error_model;
    const clade* _p_tree;
    const double _sigma_multiplier;

    clademap<Eigen::VectorXd> _probabilities;

    void compute_all_probabilities(const gene_transcript& gf, int upper_bound);
public:
    inference_pruner(const matrix_cache& cache,
        const sigma_squared* sigma,
        const error_model* p_error_model,
        const clade* _p_tree,
        double _sigma_multiplier);
    inference_pruner(const sigma_squared* sigma, const clade* p_tree, const matrix_cache* p_cache) :
        inference_pruner(*p_cache, sigma, nullptr, p_tree, 1.0)
    {

    }

    std::vector<double> prune(const gene_transcript& gf, int upper_bound);

    clademap<node_reconstruction> reconstruct(const gene_transcript& gf, int upper_bound);
};

#endif
