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

using namespace std;

extern std::mt19937 randomizer_engine;

//! Read user provided gene family data (whose path is stored in input_parameters instance)
/// @param[in] my_input_parameters Parsed parameters passed to the application
/// @param[in] p_tree The tree to be used in calculations. Necessary for syncing tree data to gene family data
/// @param[out] p_gene_families Parsed data in the gene family file specified by my_input_parameters
/// @param[out] max_family_size Equal to the largest family size given in the file plus 20%, or plus 50 if the largest family size is more than 250
/// @param[out] max_root_family_size Equal to 5/4 the size of the largest family size given in the file (with a minimum of 30)
void user_data::read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size, clade *p_tree, std::vector<gene_transcript> *p_gene_families) {

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
    
    // Iterating over gene families to get max gene family size
    for (std::vector<gene_transcript>::iterator it = p_gene_families->begin(); it != p_gene_families->end(); ++it) {
        double this_family_max_size = it->get_max_expression_value();

        if (max_family_size < this_family_max_size)
            max_family_size = this_family_max_size;
    }

    max_root_family_size = std::max(30.0, max_family_size*1.25);
    max_family_size = max_family_size + std::max(50.0, max_family_size / 5.0);
    // cout << "Read input file " << my_input_parameters.input_file_path << "." << endl;
    // cout << "Max (parsed) family size is: " << max_family_size << endl;
    // cout << "Max root family size is: " << max_root_family_size << endl;
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
sigma* user_data::read_lambda(const input_parameters &my_input_parameters, clade *p_lambda_tree) {

    sigma*p_lambda = NULL; // lambda is an abstract class, and so we can only instantiate it as single_lambda or multiple lambda -- therefore initializing it to NULL

                             // -l
    if (my_input_parameters.fixed_lambda > 0.0) {
        p_lambda = new sigma(my_input_parameters.fixed_lambda);
        // call_viterbi(max_family_size, max_root_family_size, 15, p_lambda, *p_gene_families, p_tree);
    }

    // -m
    if (!my_input_parameters.fixed_multiple_lambdas.empty()) {
        map<std::string, int> node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map(); // allows matching different lambda values to nodes in lambda tree
        vector<string> lambdastrings = tokenize_str(my_input_parameters.fixed_multiple_lambdas, ',');
        vector<double> lambdas(lambdastrings.size());

        // transform is like R's apply (vector lambdas takes the outputs, here we are making doubles from strings
        transform(lambdastrings.begin(), lambdastrings.end(), lambdas.begin(),
            [](string const& val) { return stod(val); } // this is the equivalent of a Python's lambda function
        );

        p_lambda = new sigma(node_name_to_lambda_index, lambdas);
    }

    return p_lambda;
}

//! Populate famdist_map with root family distribution read from famdist_file_path
void user_data::read_rootdist(string rootdist_file_path) {

    ifstream rootdist_file(rootdist_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    if (!rootdist_file.is_open())
        throw std::runtime_error("Failed to open file '" + rootdist_file_path + "'");
    string line;
    while (getline(rootdist_file, line)) {
        istringstream ist(line);
        int fam_size, fam_count;
        ist >> fam_size >> fam_count;
        rootdist[fam_size] = fam_count;
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
        read_gene_family_data(my_input_parameters, max_family_size, max_root_family_size, p_tree, &gene_families);
    }

    /* -e */
    if (!my_input_parameters.error_model_file_path.empty()) {
        p_error_model = new error_model;
        read_error_model(my_input_parameters, p_error_model);
    }

    /* -y */
    if (!my_input_parameters.lambda_tree_file_path.empty()) {
        p_lambda_tree = read_lambda_tree(my_input_parameters);
		p_tree->validate_lambda_tree(p_lambda_tree);
    }

    /* -l/-m (in the absence of -l, estimate) */
    p_lambda = read_lambda(my_input_parameters, p_lambda_tree);

    if (!my_input_parameters.rootdist.empty())
        read_rootdist(my_input_parameters.rootdist);
}

/// Root distributions are affected by three parameters: -p, -i, -f
/// If a rootdist file is specified (-f), those values will be used for the root distribution and the other flags
/// are ignored. Otherwise, if a poisson distribution is specified (-p) with a value, the root
/// distribution will be based on that poisson distribution, with the provided mean. If no poisson value
/// is specified, a family file must be given (-i) and those families will be used to calculate a
/// poisson distribution. If a Poisson distribution is used, values above a max family size
/// will be considered to be 0. The max family size defaults to 100, but is calculated from
/// family file if one is given.

/// priority: -p specified on command line  (poisson with specified value)
///           -f specified on command line  (specified root distribution to use)
///           -i + -p specified on command line  (poisson estimated from file)
///           -i Uniform distribution
void user_data::create_prior(const input_parameters& params)
{
    if (params.poisson_lambda > 0)
    {
        if (!rootdist.empty())
        {
            LOG(WARNING) << "\nBoth root distribution and Poisson distribution specified";
        }

        LOG(INFO) << "\nUsing Poisson root distribution with user provided lambda " << params.poisson_lambda;
        prior = root_equilibrium_distribution(params.poisson_lambda, DISCRETIZATION_RANGE);
    }
    else if (!rootdist.empty())
    {
        LOG(INFO) << "\nRoot distribution set by user provided file";
        prior = root_equilibrium_distribution(rootdist);
    }
    else if (!gene_families.empty() && params.use_poisson_dist_for_prior)
    {
        LOG(INFO) << "\nEstimating Poisson root distribution from gene families";
        prior = root_equilibrium_distribution(gene_families, DISCRETIZATION_RANGE);
    }
    else if (params.fixed_root_value > 0)
    {
        prior = root_equilibrium_distribution(params.fixed_root_value);
    }
    else 
    {
        LOG(WARNING) << "\nNo root family size distribution specified, using uniform distribution";
        prior = root_equilibrium_distribution(size_t(DISCRETIZATION_RANGE));
    }

    if (params.nsims > 0)
        prior.resize(params.nsims);
}

TEST_CASE("create_prior__creates__uniform_distribution")
{
    input_parameters params;
    user_data ud;
    ud.create_prior(params);
    gene_transcript gf;
    CHECK_EQ(doctest::Approx(0.005), ud.prior.compute(gf, 1));
    CHECK_EQ(doctest::Approx(0.005), ud.prior.compute(gf, 99));
    CHECK_EQ(0, ud.prior.compute(gf, DISCRETIZATION_RANGE+50));
}

TEST_CASE("create_prior__creates__specifed_distribution_if_given")
{
    input_parameters params;
    user_data ud;
    ud.rootdist[2] = 11;
    ud.rootdist[3] = 5;
    ud.rootdist[4] = 7;
    ud.rootdist[6] = 2;
    ud.max_root_family_size = 10;
    ud.create_prior(params);
    gene_transcript gf;
    CHECK_EQ(doctest::Approx(0.44), ud.prior.compute(gf, 2));
    CHECK_EQ(doctest::Approx(0.2), ud.prior.compute(gf, 3));
    CHECK_EQ(doctest::Approx(0.28), ud.prior.compute(gf, 4));
    CHECK_EQ(0, ud.prior.compute(gf, 5));
    CHECK_EQ(doctest::Approx(.08), ud.prior.compute(gf, 6));
}


TEST_CASE("create_prior__creates__poisson_distribution_if_given")
{
    input_parameters params;
    params.use_poisson_dist_for_prior = true;
    params.poisson_lambda = 0.75;
    user_data ud;
    ud.create_prior(params);
    gene_transcript gf;
    CHECK_EQ(doctest::Approx(0.47237f), ud.prior.compute(gf, 0));
    CHECK_EQ(doctest::Approx(0.35427f), ud.prior.compute(gf, 1));
    vector<int> r(100);
    int i = 0;
    generate(r.begin(), r.end(), [&ud, &i]() mutable { i++; return ud.prior.select_root_size(i);  });
    CHECK_EQ(94, count(r.begin(), r.end(), 1));
    CHECK_EQ(6, count(r.begin(), r.end(), 2));
    CHECK_EQ(0, count(r.begin(), r.end(), 3));
    CHECK_EQ(0, count(r.begin(), r.end(), 4));

    CHECK_EQ(0, ud.prior.compute(gf, 11));
}

TEST_CASE("create_prior__creates__poisson_distribution_if_given_distribution_and_poisson")
{
    input_parameters params;
    params.use_poisson_dist_for_prior = true;
    params.poisson_lambda = 0.75;
    user_data ud;
    ud.max_root_family_size = 100;
    ud.create_prior(params);
    gene_transcript gf;
    CHECK_EQ(doctest::Approx(0.35427f), ud.prior.compute(gf, 1));

}

TEST_CASE("create_prior__creates__poisson_distribution_from_families")
{
    randomizer_engine.seed(10);

    input_parameters params;
    params.use_poisson_dist_for_prior = true;
    user_data ud;
    ud.gene_families.resize(1);
    ud.create_prior(params);
    // bogus value that serves the purpose
    // depends on optimizer and poisson_scorer
    gene_transcript gf;
    CHECK_EQ(doctest::Approx(0.74174f), ud.prior.compute(gf, 0));
}

TEST_CASE("create_prior creates uniform distribution if poisson not specified")
{
    randomizer_engine.seed(10);

    input_parameters params;
    params.use_poisson_dist_for_prior = false;
    user_data ud;
    ud.gene_families.resize(1);
    ud.create_prior(params);
    // bogus value that serves the purpose
    // depends on optimizer and poisson_scorer
    gene_transcript gf;
    CHECK_EQ(doctest::Approx(0.005), ud.prior.compute(gf, 1));
    CHECK_EQ(doctest::Approx(0.005), ud.prior.compute(gf, 10));
    CHECK_EQ(doctest::Approx(0.005), ud.prior.compute(gf, 50));
    CHECK_EQ(doctest::Approx(0.005), ud.prior.compute(gf, 100));
}

TEST_CASE("create_prior__resizes_distribution_if_nsims_specified")
{
    randomizer_engine.seed(10);

    input_parameters params;
    params.nsims = 10;
    user_data ud;
    ud.max_root_family_size = 100;
    ud.create_prior(params);
    /// GCC's implementation of shuffle changed so the numbers that are
    /// returned are in a slightly different order, even with the same seed
#if __GNUC__ >= 7
    CHECK_EQ(152, ud.prior.select_root_size(9));
#else
    CHECK_EQ(80, ud.prior.select_root_size(9));
#endif
    CHECK_EQ(0, ud.prior.select_root_size(10));
}

