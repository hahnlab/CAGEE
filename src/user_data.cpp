#include "user_data.h"

#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <random>

#include "doctest.h"
#include "easylogging++.h"

#include "arguments.h"
#include "sigma.h"
#include "clade.h"
#include "error_model.h"
#include "io.h"
#include "root_equilibrium_distribution.h"

using namespace std;

extern std::mt19937 randomizer_engine;

void user_data::read_gene_transcript_data(const input_parameters &my_input_parameters, clade *p_tree, std::vector<gene_transcript> *p_transcripts) {

    try
    {
        ifstream input_file(my_input_parameters.input_file_path);
        if (!input_file.is_open())
            throw std::runtime_error("Failed to open");

        read_gene_transcripts(input_file, p_tree, *p_transcripts); // in io.cpp/io.h
    }
    catch (runtime_error& err)
    {
        throw std::runtime_error(my_input_parameters.input_file_path + ": " + err.what() + ". Exiting...");
    }
    
    LOG(DEBUG) << "Read input file " << my_input_parameters.input_file_path << ".";
}

//! Read user provided error model file (whose path is stored in input_parameters instance)
void user_data::read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model) {

    ifstream error_model_file(my_input_parameters.error_model_file_path);
    if (!error_model_file.is_open()) {
        throw std::runtime_error("Failed to open " + my_input_parameters.error_model_file_path + ". Exiting...");
    }

    read_error_model_file(error_model_file, p_error_model);

} 

void user_data::read_replicate_model(const input_parameters& my_input_parameters) {

    ifstream replicate_model_file(my_input_parameters.replicate_model_file_path);
    if (!replicate_model_file.is_open()) {
        throw std::runtime_error("Failed to open " + my_input_parameters.replicate_model_file_path + ". Exiting...");
    }

    this->p_replicate_model = read_replicate_model_file(replicate_model_file);
}

  //! Read user provided phylogenetic tree (whose path is stored in input_parameters instance)
clade * user_data::read_input_tree(const input_parameters &my_input_parameters) {
    return read_tree(my_input_parameters.tree_file_path, false);
}

clade * user_data::read_sigma_tree(const input_parameters &my_input_parameters) {
    return read_tree(my_input_parameters.sigma_tree_file_path, true);
}

sigma_squared* user_data::read_sigma(const input_parameters &my_input_parameters, clade *p_sigma_tree) {

    sigma_squared*p_sigma = NULL; // sigma is an abstract class, and so we can only instantiate it as single_sigma or multiple sigma -- therefore initializing it to NULL

                             // -l
    if (my_input_parameters.fixed_sigma > 0.0) {
        p_sigma = new sigma_squared(my_input_parameters.fixed_sigma);
    }

    // -m
    if (!my_input_parameters.fixed_multiple_sigmas.empty()) {
        map<std::string, int> node_name_to_sigma_index = p_sigma_tree->get_sigma_index_map(); // allows matching different sigma values to nodes in sigma tree
        vector<string> sigmastrings = tokenize_str(my_input_parameters.fixed_multiple_sigmas, ',');
        vector<double> sigmas(sigmastrings.size());

        transform(sigmastrings.begin(), sigmastrings.end(), sigmas.begin(), [](string const& val) { return stod(val); });

        p_sigma = new sigma_squared(node_name_to_sigma_index, sigmas, sigma_type::lineage_specific);
    }

    return p_sigma;
}

void user_data::read_rootdist(string rootdist_file_path) {

    ifstream rootdist_file(rootdist_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    if (!rootdist_file.is_open())
        throw std::runtime_error("Failed to open file '" + rootdist_file_path + "'");
    string line;
    while (getline(rootdist_file, line)) {
        istringstream ist(line);
        pair<float, int> values;
        ist >> values.first >> values.second;
        rootdist.push_back(values);
    }
}


void user_data::read_datafiles(const input_parameters& my_input_parameters)
{
    if (!my_input_parameters.log_config_file.empty())
    {
        el::Configurations conf(my_input_parameters.log_config_file.c_str());
        el::Loggers::reconfigureAllLoggers(conf);
    }
    if (my_input_parameters.verbose_logging_level > 0 )
    {
        el::Loggers::setVerboseLevel(my_input_parameters.verbose_logging_level);
    }
    /* -t */
    if (!my_input_parameters.tree_file_path.empty()) {
        p_tree = read_input_tree(my_input_parameters); // populates p_tree (pointer to phylogenetic tree)
    }

    /* -i */
    if (!my_input_parameters.input_file_path.empty()) {
        read_gene_transcript_data(my_input_parameters, p_tree, &gene_transcripts);
    }

    /* -e */
    if (!my_input_parameters.error_model_file_path.empty()) {
        p_error_model = new error_model;
        read_error_model(my_input_parameters, p_error_model);
    }

    if (!my_input_parameters.replicate_model_file_path.empty())
    {
        read_replicate_model(my_input_parameters);
    }
    /* -y */
    if (!my_input_parameters.sigma_tree_file_path.empty()) {
        p_sigma_tree = read_sigma_tree(my_input_parameters);
		p_tree->validate_sigma_tree(p_sigma_tree);
    }

    /* -l/-m (in the absence of -l, estimate) */
    p_sigma = read_sigma(my_input_parameters, p_sigma_tree);

    if (!my_input_parameters.rootdist_params.empty())
    {
        auto rootdist = tokenize_str(my_input_parameters.rootdist_params, ':');
        if (rootdist[0] == "file")
            read_rootdist(rootdist[1]);
    }
}

