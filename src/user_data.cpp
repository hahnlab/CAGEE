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

//! Read user provided gene family data (whose path is stored in input_parameters instance)
/// @param[in] my_input_parameters Parsed parameters passed to the application
/// @param[in] p_tree The tree to be used in calculations. Necessary for syncing tree data to gene family data
/// @param[out] p_gene_families Parsed data in the gene family file specified by my_input_parameters
void user_data::read_gene_family_data(const input_parameters &my_input_parameters, clade *p_tree, std::vector<gene_transcript> *p_gene_families) {

    try
    {
        ifstream input_file(my_input_parameters.input_file_path);
        if (!input_file.is_open())
            throw std::runtime_error("Failed to open");

        read_gene_families(input_file, p_tree, *p_gene_families); // in io.cpp/io.h
    }
    catch (runtime_error& err)
    {
        throw std::runtime_error(my_input_parameters.input_file_path + ": " + err.what() + ". Exiting...");
    }
    
    LOG(DEBUG) << "Read input file " << my_input_parameters.input_file_path << "." << endl;
}

//! Read user provided error model file (whose path is stored in input_parameters instance)
void user_data::read_error_model(const input_parameters &my_input_parameters, error_model *p_error_model) {

    ifstream error_model_file(my_input_parameters.error_model_file_path);
    if (!error_model_file.is_open()) {
        throw std::runtime_error("Failed to open " + my_input_parameters.error_model_file_path + ". Exiting...");
    }

    read_error_model_file(error_model_file, p_error_model);

} // GOTTA WRITE THIS!

  //! Read user provided phylogenetic tree (whose path is stored in input_parameters instance)
clade * user_data::read_input_tree(const input_parameters &my_input_parameters) {
    return read_tree(my_input_parameters.tree_file_path, false);
}

//! Read user provided lambda tree (lambda structure)
clade * user_data::read_lambda_tree(const input_parameters &my_input_parameters) {
    return read_tree(my_input_parameters.lambda_tree_file_path, true);
}

//! Read user provided single or multiple lambdas
sigma_squared* user_data::read_lambda(const input_parameters &my_input_parameters, clade *p_lambda_tree) {

    sigma_squared*p_sigma = NULL; // sigma is an abstract class, and so we can only instantiate it as single_sigma or multiple sigma -- therefore initializing it to NULL

                             // -l
    if (my_input_parameters.fixed_lambda > 0.0) {
        p_sigma = new sigma_squared(my_input_parameters.fixed_lambda);
        // call_viterbi(max_family_size, max_root_family_size, 15, p_lambda, *p_gene_families, p_tree);
    }

    // -m
    if (!my_input_parameters.fixed_multiple_lambdas.empty()) {
        map<std::string, int> node_name_to_sigma_index = p_lambda_tree->get_sigma_index_map(); // allows matching different sigma values to nodes in sigma tree
        vector<string> sigmastrings = tokenize_str(my_input_parameters.fixed_multiple_lambdas, ',');
        vector<double> sigmas(sigmastrings.size());

        transform(sigmastrings.begin(), sigmastrings.end(), sigmas.begin(), [](string const& val) { return stod(val); });

        p_sigma = new sigma_squared(node_name_to_sigma_index, sigmas, sigma_type::lineage_specific);
    }

    return p_sigma;
}

//! Populate famdist_map with root family distribution read from famdist_file_path
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
        // Populates (pointer to) vector of gene family instances, max_family_size and max_root_family_size (last two passed by reference)
        read_gene_family_data(my_input_parameters, p_tree, &gene_families);
    }

    /* -e */
    if (!my_input_parameters.error_model_file_path.empty()) {
        p_error_model = new error_model;
        read_error_model(my_input_parameters, p_error_model);
    }

    /* -y */
    if (!my_input_parameters.lambda_tree_file_path.empty()) {
        p_lambda_tree = read_lambda_tree(my_input_parameters);
		p_tree->validate_sigma_tree(p_lambda_tree);
    }

    /* -l/-m (in the absence of -l, estimate) */
    p_lambda = read_lambda(my_input_parameters, p_lambda_tree);

    if (!my_input_parameters.rootdist_params.empty())
    {
        auto rootdist = tokenize_str(my_input_parameters.rootdist_params, ':');
        if (rootdist[0] == "file")
            read_rootdist(rootdist[1]);
    }
}

