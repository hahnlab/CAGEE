#ifndef TRANSCRIPT_RECONSTRUCTOR_H
#define TRANSCRIPT_RECONSTRUCTOR_H

#include "gene_transcript.h"
#include "Eigen/Dense"

class replicate_model;

struct node_reconstruction {
    double most_likely_value;
    std::pair<double, double> credible_interval;
};

class reconstruction {
private:
    Eigen::IOFormat _vector_fmt{Eigen::IOFormat(3, 0, ", ", "", "[ ", " ]", "", "", (char)32)};

// #ifdef VECTOR_DEBUG
//     bool _vprint = true;
// #else
//     bool _vprint = false;
// #endif

public:
    void print_reconstructed_states(std::ostream& ost, transcript_vector& transcripts, const clade* p_tree);

    void print_increases_decreases_by_clade(std::ostream& ost, const clade* p_tree, transcript_vector& transcripts, bool count_all_changes);

    void write_results(std::string output_prefix, const clade* p_tree, transcript_vector& families, bool count_all_changes);

    reconstruction() : _p_replicates(nullptr) {};
    reconstruction(replicate_model* p_replicate) : _p_replicates(p_replicate) {};

    virtual ~reconstruction()
    {
    }

    node_reconstruction get_reconstructed_value(const gene_transcript& gf, const clade* c) const;

    double get_difference_from_parent(const gene_transcript& gf, const clade* c) const;

    double get_leaf_node_value(const gene_transcript& gf, const clade* c) const;

private:
    virtual void print_additional_data(transcript_vector& transcripts, std::string output_prefix) {};

    virtual void write_nexus_extensions(std::ostream& ost) {};

    virtual node_reconstruction get_internal_node_value(const gene_transcript& gf, const clade* c) const = 0;

    replicate_model* _p_replicates;
};



#endif
