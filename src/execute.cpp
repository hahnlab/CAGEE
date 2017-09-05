#include "io.h"
#include "execute.h"

//! Read user provided gene family data (whose path is stored in input_parameters instance)
std::vector<gene_family> * execute::read_gene_family_data(const input_parameters &my_input_parameters, int &max_family_size, int &max_root_family_size) {
    
    std::vector<gene_family> *p_gene_families;
    if (!my_input_parameters.input_file_path.empty()) {
        p_gene_families = read_gene_families(my_input_parameters.input_file_path);
            
        // Iterating over gene families to get max gene family size
        for (std::vector<gene_family>::iterator it = p_gene_families->begin(); it != p_gene_families->end(); ++it) {
            int this_family_max_size = it->get_parsed_max_size();
            
            if (max_family_size < this_family_max_size)
                max_family_size = this_family_max_size;
            }

            cout << max_family_size << endl;
        }
    
    return p_gene_families;
}

//! Read user provided phylogenetic tree (whose path is stored in input_parameters instance)
clade * execute::read_input_tree(const input_parameters &my_input_parameters, bool read_lambdas) {
    clade *p_tree = new clade();
    p_tree = read_tree(my_input_parameters.tree_file_path, read_lambdas);
    
    return p_tree;
}