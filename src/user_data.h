#ifndef USER_DATA_H
#define USER_DATA_H

#include <stddef.h>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "gene_transcript.h"
#include "root_equilibrium_distribution.h"

class root_equilibrium_distribution;
class clade;
class sigma_squared;
class error_model;
struct input_parameters;

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

    std::vector<gene_transcript> gene_transcripts;
    std::vector<std::pair<float, int>> rootdist;    // first is value, second is count

    void read_datafiles(const input_parameters& my_input_parameters);

    //! Read in gene transcript data
    void read_gene_transcript_data(const input_parameters &my_input_parameters, clade *p_tree, std::vector<gene_transcript> *p_gene_transcripts);

    //! Read in error model file
    void read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model);

    //! Read in phylogenetic tree data
    clade * read_input_tree(const input_parameters &my_input_parameters);

    //! Read in sigma tree
    clade * read_sigma_tree(const input_parameters &my_input_parameters);

    //! Read in single or multiple sigma
    sigma_squared* read_sigma(const input_parameters &my_input_parameters, clade *p_sigma_tree);

    void read_rootdist(std::string rootdist_file_path);
};

#endif
