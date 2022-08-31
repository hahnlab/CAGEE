#ifndef TRANSCRIPT_RECONSTRUCTOR_H
#define TRANSCRIPT_RECONSTRUCTOR_H

#include "gene_transcript.h"

//! The result of a model reconstruction. Should be able to (a) print reconstructed states with all available information;
/// (b) print increases and decreases by family; and (c) print increases and decreases by clade.
class reconstruction {
public:
    typedef const std::vector<gene_transcript> transcript_vector;

    void print_node_change(std::ostream& ost, const cladevector& order, transcript_vector& transcripts, const clade* p_tree);

    void print_node_values(std::ostream& ost, const cladevector& order, transcript_vector& transcripts, const clade* p_tree);

    void print_reconstructed_states(std::ostream& ost, transcript_vector& transcripts, const clade* p_tree, double test_pvalue);

    void print_increases_decreases_by_clade(std::ostream& ost, const clade* p_tree, transcript_vector& transcripts);

    void print_family_clade_table(std::ostream& ost, const cladevector& order, transcript_vector& transcripts, const clade* p_tree,
        std::function<std::string(int family_index, const clade* c)> get_family_clade_value);

    void write_results(std::string model_identifier, std::string output_prefix, const clade* p_tree, transcript_vector& families, double test_pvalue);

    virtual ~reconstruction()
    {
    }

    virtual double get_node_value(const gene_transcript& gf, const clade* c) const = 0;

    double get_difference_from_parent(const gene_transcript& gf, const clade* c);
private:
    virtual void print_additional_data(transcript_vector& transcripts, std::string output_prefix) {};

    virtual void write_nexus_extensions(std::ostream& ost) {};

};


#endif
