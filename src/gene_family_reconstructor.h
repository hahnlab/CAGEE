#ifndef RECONSTRUCTION_PROCESS_H
#define RECONSTRUCTION_PROCESS_H

#include "core.h"
#include <map>

class matrix_cache;
class sigma;

namespace pupko_reconstructor {

void reconstruct_leaf_node(const clade * c, const sigma *p_sigma, const gene_transcript& t, clademap<std::vector<double>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, const matrix_cache *_p_calc);
void reconstruct_at_node(const clade *c, const gene_transcript& t, const sigma*_lambda, clademap<std::vector<double>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, const matrix_cache* p_calc);
void reconstruct_internal_node(const clade * c, const sigma* p_sigma, const gene_transcript& t, clademap<std::vector<double>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, const matrix_cache *_p_calc);
void initialize_at_node(const clade* c, clademap<std::vector<double>>& all_node_Cs, clademap<std::vector<double>>& all_node_Ls, int max_family_size, int max_root_family_size);

/// Structure holding required memory to calculate Pupko reconstructions
struct pupko_data {
    std::vector<clademap<std::vector<double>>> v_all_node_Cs;
    std::vector<clademap<std::vector<double>>> v_all_node_Ls;

    pupko_data(size_t num_families, const clade* p_tree, int max_family_size, int max_root_family_size);

    clademap<std::vector<double>>& C(size_t i) { return v_all_node_Cs.at(i); }
    clademap<std::vector<double>>& L(size_t i) { return v_all_node_Ls.at(i); }
};


/// Given a gene gamily and a tree, reconstructs the most likely values at each node on tree. Used in the base model to calculate values for each
/// gene family. Also used in a gamma bundle, one for each gamma category. Differences are represented by the lambda multiplier.
void reconstruct_gene_transcript(const sigma* lambda, const clade *p_tree,
    const gene_transcript *gf,
    matrix_cache *p_calc,
    clademap<int>& reconstructed_states,
    clademap<std::vector<double>>& all_node_Cs,
    clademap<std::vector<double>>& all_node_Ls);
}


branch_probabilities::branch_probability compute_viterbi_sum(const clade* c, const gene_transcript& family, const reconstruction* rec, int max_family_size, const matrix_cache& cache, const sigma* p_lambda);

void print_branch_probabilities(std::ostream& ost, const cladevector& order, const std::vector<gene_transcript>& gene_families, const branch_probabilities& branch_probabilities);

#endif
