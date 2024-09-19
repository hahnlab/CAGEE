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
class replicate_model;
class matrix_cache;

size_t adjust_for_error_model(size_t c, const error_model *p_error_model);

template<typename T>
using clademap = std::map<const clade*, T>;

using boundaries = std::pair<double, double>;

class optional_probabilities
{
    bool has_value = false;
    Eigen::VectorXd _probabilities;
public:
    void setOne();
    void clear() { has_value = false; }
    const Eigen::VectorXd& probabilities() const;
    bool hasValue() const { return has_value; }
    void initialize(double transcript_value, boundaries bounds);
    void reserve(const Eigen::VectorXd& v) { _probabilities = v;  }
    size_t capacity() const { return _probabilities.size(); }
    void multiply_elements(const Eigen::VectorXd& multipliers);

    void set(const Eigen::VectorXd& v)
    {
        _probabilities = v;
        has_value = true;
    }
};

class inference_pruner
{
    const matrix_cache& _cache;
    const sigma_squared* _p_sigsqd;
    const error_model* _p_error_model;
    const replicate_model* _p_replicate_model;
    const clade* _p_tree;
    const boundaries _bounds;

    clademap<optional_probabilities> _probabilities;

    void compute_all_probabilities(const gene_transcript& gf);
public:
    inference_pruner(const matrix_cache& cache,
        const sigma_squared* sigma,
        const error_model* p_error_model,
        const replicate_model* p_replicate_model,
        const clade* _p_tree,
        boundaries bounds);

    inference_pruner(const sigma_squared* sigma, const clade* p_tree, const replicate_model* p_replicate_model, const matrix_cache* p_cache, boundaries bounds) :
        inference_pruner(*p_cache, sigma, nullptr, p_replicate_model, p_tree, bounds)
    {

    }

    std::vector<double> prune(const gene_transcript& gf);

    clademap<node_reconstruction> reconstruct(const gene_transcript& gf);
};

void compute_node_probability(const clade* node,
    const gene_transcript& gene_transcript,
    const error_model* p_error_model,
    const replicate_model* p_replicate_model,
    std::map<const clade*, optional_probabilities>& probabilities,
    const sigma_squared* p_sigma,
    const matrix_cache& cache,
    boundaries bounds);

#endif
