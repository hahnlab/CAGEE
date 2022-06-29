#ifndef TRANSCRIPT_RECONSTRUCTOR_H
#define TRANSCRIPT_RECONSTRUCTOR_H

#include "core.h"
#include <map>

#include <Eigen/Core>

class matrix_cache;
class sigma_squared;

/// Given a gene gamily and a tree, reconstructs the most likely values at each node on tree. Used in the base model to calculate values for each
/// gene family. Also used in a gamma bundle, one for each gamma category. Differences are represented by the lambda multiplier.
class transcript_reconstructor
{
    clademap<Eigen::VectorXd> all_node_Ls;
    const sigma_squared* _p_sigma;
    const clade* _p_tree;
    const matrix_cache* _p_cache;
public:
    transcript_reconstructor(const sigma_squared* lambda, const clade* p_tree, const matrix_cache* p_cache);

    clademap<double> reconstruct_gene_transcript(const gene_transcript& gf, int upper_bound);
};

branch_probabilities::branch_probability compute_viterbi_sum(const clade* c, const gene_transcript& family, const reconstruction* rec, const matrix_cache& cache, const sigma_squared* p_lambda, int upper_bound);

void print_branch_probabilities(std::ostream& ost, const cladevector& order, const std::vector<gene_transcript>& gene_families, const branch_probabilities& branch_probabilities);

#endif
