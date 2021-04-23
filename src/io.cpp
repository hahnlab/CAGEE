#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <sstream>
#include <cstring>
#include <set>
#include <numeric>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include <sys/stat.h>

#include "doctest.h"

#include "io.h"
#include "gene_transcript.h"
#include "error_model.h"
#include "clade.h"

using namespace std;

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "error_model", optional_argument, NULL, 'e' },
  { "output_prefix", required_argument, NULL, 'o'}, 
  { "tree", required_argument, NULL, 't' },
  { "fixed_sigma", required_argument, NULL, 'l' },
  { "fixed_multiple_sigmas", required_argument, NULL, 'm' },
  { "sigma_tree", required_argument, NULL, 'y' },
  { "n_gamma_cats", required_argument, NULL, 'k' },
  { "fixed_alpha", required_argument, NULL, 'a' },
  { "rootdist", required_argument, NULL, 'f'},
  { "fixed_root_value", required_argument, NULL, 'F'},
  { "poisson", optional_argument, NULL, 'p' },
  { "simulate", optional_argument, NULL, 's' },
  { "pvalue", required_argument, NULL, 'P' },
  { "zero_root", no_argument, NULL, 'z' },
  { "cores", required_argument, NULL, 'c' },
  { "lambda_per_family", no_argument, NULL, 'b' },
  { "log_config", required_argument, NULL, 'L' },
  { "optimizer_expansion", optional_argument, NULL, 'E' },
  { "optimizer_reflection", optional_argument, NULL, 'R' },
  { "optimizer_iterations", optional_argument, NULL, 'I' },
  { "help", no_argument, NULL, 'h'},
  { 0, 0, 0, 0 }
};

void input_parameters::check_input() {
    //! Options -l and -m cannot both specified.
    if (fixed_lambda > 0.0 && !fixed_multiple_lambdas.empty()) {
        throw runtime_error("Options -l and -m are mutually exclusive.");
    }

    //! Option -m requires a lambda tree (-y)
    if (!fixed_multiple_lambdas.empty() && lambda_tree_file_path.empty()) {
        throw runtime_error("Multiple lambda values (-m) specified with no lambda tree (-y)");
    }
    
    //! Options -l and -i have to be both specified (if estimating and not simulating).
    if (fixed_lambda > 0.0 && input_file_path.empty() && !is_simulating) {
        throw runtime_error("Options -l and -i must both be provided an argument.");
    }
    
    if (is_simulating)
    {
        // Must specify a lambda
        if (fixed_lambda <= 0.0 && fixed_multiple_lambdas.empty()) {
            throw runtime_error("Cannot simulate without initial sigma values");
        }

        if (fixed_alpha <= 0.0 && this->n_gamma_cats > 1) {
            throw runtime_error("Cannot simulate gamma clusters without an alpha value");
        }

        //! Options -i and -f cannot be both specified. Either one or the other is used to specify the root eq freq distr'n.
        if (!input_file_path.empty() && !rootdist.empty()) {
            throw runtime_error("Options -i and -f are mutually exclusive.");
        }
    }
    else
    {
        if (fixed_alpha >= 0.0 && n_gamma_cats == 1) {
            throw runtime_error("Alpha specified with 1 gamma category.");
    }


    if (lambda_per_family)
    {
        if (input_file_path.empty())
            throw runtime_error("No family file provided");
        if (tree_file_path.empty())
            throw runtime_error("No tree file provided");
    }

    if (n_gamma_cats > 1 && use_error_model && error_model_file_path.empty())
    {
        throw runtime_error("Estimating an error model with a gamma distribution is not supported at this time");
    }

    }
}

/* START: Reading in tree data */
//! Read tree from user-provided tree file
/*!
  This function is called by CAFExp's main function when "--tree"/"-t" is specified
*/
clade* read_tree(string tree_file_path, bool lambda_tree) {
    ifstream tree_file(tree_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
    if (!tree_file.is_open())
    {
        throw std::runtime_error("Failed to open " + tree_file_path);
    }

    string line;
    
    if (tree_file.good()) {
        getline(tree_file, line);
    }
    tree_file.close();
    
    clade *p_tree = parse_newick(line, lambda_tree);
    
    if (p_tree->is_leaf())
        throw std::runtime_error(tree_file_path + " does not seem to be a valid tree");

    return p_tree;
}
/* END: Reading in tree data */

//! Read gene family data from user-provided tab-delimited file
/*!
  This function is called by execute::read_gene_transcript_data, which is itself called by CAFExp's main function when "--infile"/"-i" is specified  
*/
void read_gene_families(std::istream& input_file, clade *p_tree, std::vector<gene_transcript> &gene_families) {
    map<int, std::string> sp_col_map; // For dealing with CAFE input format, {col_idx: sp_name} 
    std::string line;
    bool is_header = true;
    int first_gene_index = 2;

    while (getline(input_file, line)) {
        if (line.empty()) continue;

        std::vector<std::string> tokens = tokenize_str(line, '\t');

        // If still reading header
        if (is_header) {
            if (tokens[2] == "Treatment_Tissue")
                first_gene_index = 3;
            is_header = false;
                
            for (size_t i = first_gene_index; i < tokens.size(); ++i) {
                sp_col_map[i] = tokens[i];
            }
        }
        
        // Header has ended, reading gene family counts
        else {
            gene_transcript genfam(tokens[0], tokens[1], first_gene_index == 3 ? tokens[2] : "");
            
            for (size_t i = first_gene_index; i < tokens.size(); ++i) {
                std::string sp_name = sp_col_map[i];
                genfam.set_expression_value(sp_name, atof(tokens[i].c_str()));
            }
            
            gene_families.push_back(genfam);
        }
    }

    if (gene_families.empty())
        throw std::runtime_error("No families found");
}
/* END: Reading in gene family data */

double to_double(string s)
{
    return std::stod(s);
}
//! Read user-provided error model
/*!
  This function is called by execute::read_error_model, which is itself called by CAFExp's main function when "--error_model"/"-e" is specified  
*/
void read_error_model_file(std::istream& error_model_file, error_model *p_error_model) {
    std::string line;
    std::string max_header = "max";
    std::string cnt_diff_header = "cnt";
    
    while (getline(error_model_file, line)) {
        std::vector<std::string> tokens;
        
        // maxcnt line
        if (strncmp(line.c_str(), max_header.c_str(), max_header.size()) == 0) { 
            tokens = tokenize_str(line, ':');
            tokens[1].erase(remove_if(tokens[1].begin(), tokens[1].end(), ::isspace), tokens[1].end()); // removing whitespace
            int max_cnt = std::stoi(tokens[1]);
            p_error_model->set_max_family_size(max_cnt);
        }
        
        // cntdiff line
        else if (strncmp(line.c_str(), cnt_diff_header.c_str(), cnt_diff_header.size()) == 0) { 
            tokens = tokenize_str(line, ' ');
            
            if (tokens.size() % 2 != 0) { 
                throw std::runtime_error("Number of different count differences in the error model (including 0) is not an odd number. Exiting...");
            }
            
            std::vector<std::string> cnt_diffs(tokens.begin()+1, tokens.end());
            
//            cout << "Count diffs are" << endl;
//            for (auto it = cnt_diffs.begin(); it != cnt_diffs.end(); ++it) {
//                cout << *it << " ";
//            }
//            cout << endl;
            p_error_model->set_deviations(cnt_diffs);
        }
        else
        {
            tokens = tokenize_str(line, ' ');
            if (tokens.size() > 0)
            {
                int sz = std::stoi(tokens[0]);
                vector<double> values(tokens.size() - 1);
                transform(tokens.begin() + 1, tokens.end(), values.begin(), to_double);
                p_error_model->set_probabilities(sz, values);
            }
        }
    }
        // std::vector<std::string> tokens = tokenize_str(line, '\t'); 
}
/* END: Reading in error model data */

void write_error_model_file(std::ostream& ost, error_model& errormodel)
{
    ost << "maxcnt: " << errormodel.get_max_family_size()-1 << "\n";
    ost << "cntdiff:";
    for (int j : errormodel._deviations) {
        ost << " " << j;
    }
    ost << "\n";

    vector<double> last_probs;
    for (size_t j = 0; j < errormodel.get_max_family_size(); j++) {
        auto probs = errormodel.get_probs(j);
        if (probs == last_probs) continue;
        last_probs = probs;

        ost << j;
        for (auto p : probs)
            ost << " " << p;
        ost << endl;
    }
}

//! Split string into vector of strings given delimiter
std::vector<std::string> tokenize_str(std::string some_string, char some_delim) {
    std::istringstream ist(some_string);
    std::string token;
    std::vector<std::string> tokens;

    while (std::getline(ist, token, some_delim)) {
        tokens.push_back(token);
    }

    return tokens;
}

/// OS-specific. If mkdir succeeds it returns 0. If it returns -1
/// check to see if it failed because the directory already exists.
/// If it failed for some other reason, throw an exception.
void create_directory(std::string& dir)
{
    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if (errno != EEXIST)
            throw std::runtime_error("Failed to create directory");
    }
}

std::ostream& operator<<(std::ostream& ost, const gene_transcript& family)
{
    ost << family._id;
    if (!family._desc.empty())
        ost << "," << family._desc;
    ost << ":";
    for (auto& element : family._species_size_map) {
        ost << " " << element.first << ":" << element.second;
    }
    return ost;
}


input_parameters read_arguments(int argc, char* const argv[]);

struct option_test
{
    char* values[100];
    size_t argc;

    option_test(vector<string> arguments)
    {
        optind = 0;
        argc = arguments.size();
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            values[i] = strdup(arguments[i].c_str());
        }
    }

    ~option_test()
    {
        for (size_t i = 0; i < argc; ++i)
        {
            free(values[i]);
        }
    }
};

TEST_CASE("read_arguments translates short values ") {
    option_test c({ "cafe5", "-ifile" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}

TEST_CASE("read_arguments translates long values ") {
    option_test c({ "cafe5", "--infile", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}
TEST_CASE("Options, input_short_space_separated")
{
    option_test c({ "cafe5", "-i", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}
TEST_CASE("Options, simulate_long")
{
    option_test c({ "cafe5", "--simulate=1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, simulate_short")
{
    option_test c({ "cafe5", "-s1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, pvalue_long")
{
    option_test c({ "cafe5", "--pvalue=0.01" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.01, actual.pvalue);
}

TEST_CASE("Options, pvalue_short")
{
    option_test c({ "cafe5", "-P0.01" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.01, actual.pvalue);
}

TEST_CASE("Options, optimizer_long")
{
    option_test c({ "cafe5", "--optimizer_expansion=0.05", "--optimizer_reflection=3.2", "--optimizer_iterations=5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.neldermead_iterations);
}

TEST_CASE("Options, optimizer_short")
{
    option_test c({ "cafe5", "-E", "0.05", "-R", "3.2", "-I", "5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.neldermead_iterations);
}

TEST_CASE("Options, cores_long")
{
    option_test c({ "cafe5", "--cores=6" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(6, actual.cores);
}

TEST_CASE("Options, cores_short")
{
    option_test c({ "cafe5", "-c", "8" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(8, actual.cores);
}

TEST_CASE("Options: errormodel_accepts_argument")
{
    option_test c({ "cafe5", "-eerror.txt" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path == "error.txt");
}

TEST_CASE("Options: fixed_root_value")
{
    option_test c({ "cafe5", "--fixed_root_value", "12.7" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(12.7, actual.fixed_root_value);
}

TEST_CASE("Options, errormodel_accepts_no_argument")
{
    option_test c({ "cafe5", "-e" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path.empty());
}

TEST_CASE("Options: zero_root_familes")
{
    input_parameters by_default;
    CHECK_FALSE(by_default.exclude_zero_root_families);

    option_test c({ "cafe5", "-z" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.exclude_zero_root_families);
}

TEST_CASE("Options: cannot_have_space_before_optional_parameter")
{
    option_test c({ "cafe5", "-s", "1000" });

    CHECK_THROWS_WITH(read_arguments(c.argc, c.values), "Unrecognized parameter: '1000'");
}


TEST_CASE("Options: must_specify_sigma_for_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate without initial sigma values");
}

TEST_CASE("Options: must_specify_lambda_and_input_file_for_estimator")
{
    input_parameters params;
    params.fixed_lambda = 0.05;
    CHECK_THROWS_WITH(params.check_input(), "Options -l and -i must both be provided an argument.");
}

TEST_CASE("Options: must_specify_alpha_for_gamma_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_lambda = 0.05;
    params.n_gamma_cats = 3;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate gamma clusters without an alpha value");
}

TEST_CASE("Options: must_specify_alpha_and_k_for_gamma_inference")
{
    input_parameters params;
    params.fixed_alpha = 0.7;
    CHECK_THROWS_WITH(params.check_input(), "Alpha specified with 1 gamma category.");
}

TEST_CASE("Options: can_specify_alpha_without_k_for_gamma_simulation")
{
    input_parameters params;
    params.fixed_alpha = 0.7;
    params.fixed_lambda = 0.01;
    params.is_simulating = true;
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: check_input_does_not_throw_when_simulating_with_multiple_lambdas")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_multiple_lambdas = "0.01,0.05";
    params.lambda_tree_file_path = "./tree";
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: per_family_must_provide_families")
{
    input_parameters params;
    params.lambda_per_family = true;
    CHECK_THROWS_WITH_AS(params.check_input(), "No family file provided", runtime_error);
}

TEST_CASE("Options: per_family_must_provide_tree")
{
    input_parameters params;
    params.lambda_per_family = true;
    params.input_file_path = "/tmp/test";
    CHECK_THROWS_WITH_AS(params.check_input(), "No tree file provided", runtime_error);
}

TEST_CASE("Options: cannot_estimate_error_and_gamma_together")
{
    input_parameters params;
    params.n_gamma_cats = 3;
    params.use_error_model = true;
    params.error_model_file_path = "model.txt";
    params.check_input();
    CHECK(true);

    params.error_model_file_path.clear();
    CHECK_THROWS_WITH_AS(params.check_input(), "Estimating an error model with a gamma distribution is not supported at this time", runtime_error);

}

TEST_CASE("Options: Cannot specify rootdist for simulations with rootdist file and transcript file")
{
    input_parameters params;
    params.fixed_lambda = 10;
    params.input_file_path = "transcripts.txt";
    params.rootdist = "10 1";
    params.check_input();
    CHECK(true);

    params.is_simulating = true;
    CHECK_THROWS_WITH_AS(params.check_input(), "Options -i and -f are mutually exclusive.", runtime_error);

}

TEST_CASE("GeneFamilies: read_gene_families_throws_if_no_families_found")
{
    std::string empty;
    std::istringstream ist(empty);

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));
    std::vector<gene_transcript> families;
    CHECK_THROWS_WITH_AS(read_gene_families(ist, p_tree.get(), families), "No families found", runtime_error);
}

TEST_CASE("GeneFamilies: read_gene_families skips blank lines in input")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\n\n\n\t (null)1\t5\t10\t2\t6\n\n\n\n";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    CHECK_EQ(1, families.size());
}

TEST_CASE("GeneFamilies: read_gene_families_reads_cafe_files")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    CHECK_EQ(5, families.at(0).get_expression_value("A"));
    CHECK_EQ(10, families.at(0).get_expression_value("B"));
    CHECK_EQ(2, families.at(0).get_expression_value("C"));
    CHECK_EQ(6, families.at(0).get_expression_value("D"));
}

TEST_CASE("read_gene_families reads Treatment Tissue header")
{
    std::string str = "Desc\tFamily ID\tTreatment_Tissue\tA\tB\tC\tD\n\t (null)1\tovary\t5\t10\t2\t6\n\t (null)2\tlung\t5\t10\t2\t6\n\t (null)3\tbreast\t5\t10\t2\t6\n\t (null)4\tbrain\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    REQUIRE_EQ(4, families.size());

    CHECK(families[0].tissue() == "ovary");
    CHECK(families[1].tissue() == "lung");
    CHECK(families[2].tissue() == "breast");
    CHECK(families[3].tissue() == "brain");
}