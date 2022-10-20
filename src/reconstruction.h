#ifndef TRANSCRIPT_RECONSTRUCTOR_H
#define TRANSCRIPT_RECONSTRUCTOR_H

#include "gene_transcript.h"

struct node_reconstruction {
    double most_likely_value;
    std::pair<double, double> credible_interval;
};

class reconstruction {
public:
    void print_reconstructed_states(std::ostream& ost, transcript_vector& transcripts, const clade* p_tree);

    void print_increases_decreases_by_clade(std::ostream& ost, const clade* p_tree, transcript_vector& transcripts, bool count_all_changes);

    void print_family_clade_table(std::ostream& ost, const cladevector& order, transcript_vector& transcripts, const clade* p_tree,
        std::function<std::string(const gene_transcript& transcript, const clade* c)> get_family_clade_value);

    void write_results(std::string output_prefix, const clade* p_tree, transcript_vector& families, bool count_all_changes);

    virtual ~reconstruction()
    {
    }

    virtual node_reconstruction get_internal_node_value(const gene_transcript& gf, const clade* c) const = 0;

    double get_difference_from_parent(const gene_transcript& gf, const clade* c) const;


private:
    virtual void print_additional_data(transcript_vector& transcripts, std::string output_prefix) {};

    virtual void write_nexus_extensions(std::ostream& ost) {};
};



#endif
