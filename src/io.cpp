#include <string>
#include <fstream>
#include <iostream>
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
#include "proportional_variance.h"

using namespace std;
namespace pv = proportional_variance;

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
        if (line[0] == '#') continue;

        std::vector<std::string> tokens = tokenize_str(line, '\t');

        // If still reading header
        if (is_header) {
            if (tokens[2] == "SAMPLETYPE")
                first_gene_index = 3;
            is_header = false;
                
            for (size_t i = first_gene_index; i < tokens.size(); ++i) {
                sp_col_map[i] = tokens[i];
            }
        }
        
        // Header has ended, reading gene family counts
        else {
            gene_transcript genfam(tokens[1], tokens[0], first_gene_index == 3 ? tokens[2] : "");
            
            for (size_t i = first_gene_index; i < tokens.size(); ++i) {
                std::string sp_name = sp_col_map[i];
                double val = pv::to_computational_space(atof(tokens[i].c_str()));
                genfam.set_expression_value(sp_name, val);
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
#ifndef _WIN32
    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if (errno != EEXIST)
            throw std::runtime_error("Failed to create directory");
    }
#endif
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

TEST_CASE("GeneFamilies: read_gene_families skips comment lines in input")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\n#Comment1\n\n\t (null)1\t5\t10\t2\t6\n\n#Comment 2\n\n";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    CHECK_EQ(1, families.size());
}

TEST_CASE("GeneFamilies: read_gene_families_reads_cafe_files")
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\ndesc1\tid1\t5\t10\t2\t6\ndesc2\tid2\t5\t10\t2\t6\ndesc3\tid3\t5\t10\t2\t6\ndesc4\tid4\t\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    CHECK_EQ(string("id3"), families[2].id());
    CHECK_EQ(string("desc4"), families[3].description());
    CHECK_EQ(doctest::Approx(5.0), pv::to_user_space(families.at(0).get_expression_value("A")));
    CHECK_EQ(doctest::Approx(10.0), pv::to_user_space(families.at(0).get_expression_value("B")));
    CHECK_EQ(doctest::Approx(2.0), pv::to_user_space(families.at(0).get_expression_value("C")));
    CHECK_EQ(doctest::Approx(6.0), pv::to_user_space(families.at(0).get_expression_value("D")));
}

TEST_CASE("read_gene_families reads Treatment Tissue header")
{
    std::string str = "Desc\tFamily ID\tSAMPLETYPE\tA\tB\tC\tD\n\t (null)1\tovary\t5\t10\t2\t6\n\t (null)2\tlung\t5\t10\t2\t6\n\t (null)3\tbreast\t5\t10\t2\t6\n\t (null)4\tbrain\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_transcript> families;
    read_gene_families(ist, NULL, families);
    REQUIRE_EQ(4, families.size());

    CHECK(families[0].tissue() == "ovary");
    CHECK(families[1].tissue() == "lung");
    CHECK(families[2].tissue() == "breast");
    CHECK(families[3].tissue() == "brain");
}