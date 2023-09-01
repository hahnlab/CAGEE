#ifndef USER_DATA_H
#define USER_DATA_H

#include <stddef.h>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "gene_transcript.h"

class root_equilibrium_distribution;
class clade;
class sigma_squared;
class error_model;
struct input_parameters;
class replicate_model;
class prior;

using boundaries = std::pair<double, double>;

int upper_bound_from_transcript_values(const std::vector<gene_transcript>& transcripts);

typedef struct model_params {
    std::string model_type;
    double a;
    double r;
} model_params;

/// Class holding data defined by the user, or derived from data defined by the user
class user_data {
public:
    user_data()
    {

    }
    clade *p_tree = NULL; // instead of new clade(), o.w. mem leak
    sigma_squared *p_sigma = NULL;
    clade *p_sigma_tree = NULL;
    error_model *p_error_model = NULL;
    replicate_model* p_replicate_model = NULL;
    prior* p_prior = NULL;

    std::vector<gene_transcript> gene_transcripts;
    std::vector<std::pair<float, int>> rootdist;    // first is value, second is count

    boundaries bounds;

    void read_datafiles(const input_parameters& my_input_parameters);

    //! Read in gene transcript data
    void read_gene_transcript_data(const input_parameters &my_input_parameters, clade *p_tree, std::vector<gene_transcript> *p_gene_transcripts);

    //! Read in error model file
    // void read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model);

    void read_replicate_model(const input_parameters& my_input_parameters);

    //! Read in phylogenetic tree data
    clade * read_input_tree(const input_parameters &my_input_parameters);

    //! Read in sigma tree
    clade * read_sigma_tree(const input_parameters &my_input_parameters);

    //! Read in single or multiple sigma
    sigma_squared* read_sigma(const input_parameters &my_input_parameters, clade *p_sigma_tree);

    void read_rootdist(std::string rootdist_file_path);

    void update_boundaries(bool transcripts_are_ratios);
};

#endif
