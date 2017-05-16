#include <vector>
#include <map>
#include <random>
#include "clade.h"
#include "probability.h"
#include "family_generator.h"

static map<clade *, int> family_sizes; // key: clades, value: family size

std::default_random_engine gen(12);
std::uniform_real_distribution<> dis(0, 1); // draw random number from uniform distribution

/* Parameters to family size randomizer. TODO: Be more clever about passing these into set_node_familysize_random. Global variables are to be avoided */
int _max_family_size;
double _lambda;

//! Set the family size of a node to a random value, using parent's family size
/*!
  Starting from 0, the gene family size of the child (c) is increased until the cumulative probability of c (given the gene family size s of the parent) exceeds a random draw from a uniform distribution. When this happens, the last c becomes the child's gene family size.
 
  Note that the smaller the draw from the uniform, the higher the chance that c will be far away from s.
*/
//Before: void set_node_familysize_random(clade *node) {
void set_node_familysize_random(clade *node, trial *p_tth_trial) {

    if (node->is_root()) { return; } // if node is root, we do nothing

    /* Drawing random number from uniform */
    //  std::default_random_engine gen(static_cast<long unsigned int>(time(0)));
    // double rnd = dis(gen);
    double rnd = unifrnd(); // Ben: the smaller rnd is, the further away from the parent family size will c be; the larger rnd is, the closer c will be of the parent family size (?)
    // cout << "Max family size: " << _max_family_size << " and rnd = " << rnd << endl;
    // cout << "Branch length: " << node->get_branch_length() << endl;
    double cumul = 0;
    //Before: int parent_family_size = family_sizes[node->get_parent()];
    int parent_family_size = (*p_tth_trial)[node->get_parent()];
    int c = 0; // c is the family size we will go to
    
    if (parent_family_size > 0) {
        for (; c < _max_family_size - 1; c++) {
            double prob = the_probability_of_going_from_parent_fam_size_to_c(_lambda, node->get_branch_length(), parent_family_size, c);
            cumul += prob;
            // cout << "Probability value for " << c << ": " << prob << " (cumulative " << cumul << ")" << endl;
        
            if (cumul >= rnd) {
                // cout << "Stopping" << endl;
                break;
            }
        }

//    if (c == 349)
//    {
//      double to5 = the_probability_of_going_from_parent_fam_size_to_c(_lambda, node->get_branch_length(), parent_family_size, 5);
//      cout << "Setting " << node->get_taxon_name() << " to: " << c << endl;
//      cout << "Probability of going from " << parent_family_size << " to 5 is " << to5 << endl;
//    }
    }
    cout << "Setting " << node->get_taxon_name() << " to: " << c << endl;
    //Before: family_sizes[node] = c;
    (*p_tth_trial)[node] = c;
}

//! Simulate gene family evolution from a single root family size 
/*!
  Given a root gene family size and a lambda, simulates gene family counts for all nodes in the tree.
  Returns family_sizes map (= trial), key = pointer to node, value = gene family size
*/
//map<clade *, int> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda) {
vector<trial *> simulate_families_from_root_size(clade *tree, int num_trials, int root_family_size, int max_family_size, double lambda) {

    // family_sizes.clear(); // resets map {clade:fam size} (i.e., resets trial)
    trial *p_tth_trial = new trial;
    vector<trial *> result;
    _max_family_size = max_family_size; // we could pass these variables values in a different way
    _lambda = lambda;
    
    cout << "Number of simulations: " << num_trials << endl;
    
    for (int t = 0; t < num_trials; ++t) {
        
        cout << "Doing family " << t << endl;
        
        //Before: family_sizes[tree] = root_family_size; // set root family size
        (*p_tth_trial)[tree] = root_family_size;
        tree->apply_prefix_order(set_node_familysize_random); // Ben: here I'd like to pass set_node_familysize_random an instance of p_tth_trial
        //tree->apply_to_descendants(set_node_familysize_random); // setting the family size (random) of the descendants
        
        result.push_back(p_tth_trial);
    }

    //Before: return family_sizes;
    return result;
}

//! Simulate gene family evolution from a distribution of root family sizes
/*!
  Given a root family size distribution (provided as a map read using read_rootdist()), simulates gene family counts for all nodes in the tree, and the number of times specified by the distribution (e.g., 10 families with root size 2, 5 families with root size 3, etc.)
*/
/*
vector<vector<trial *> > simulate_families_from_distribution(clade *p_tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda) {
//vector<trial> simulate_families_from_distribution(clade *p_tree, int num_trials, const std::map<int, int>& root_dist, int max_family_size, double lambda) {


    //Before: vector<trial> result; // trial is a typedef for a map of kind {clade:fam size}, which is what simulate_families_from_root() returns
    vector<vector<trial *> > result;
    std::map<int, int>::const_iterator it = root_dist.begin();
    // Ben: while loop visits all keys in root_dist map (?)
    while (it != root_dist.end()) {
        int root_size = it->first;
        cout << "Root size: " << root_size << endl;
        int n_fams_this_size = it->second;
        result.push_back(simulate_families_from_root_size(p_tree, n_fams_this_size, root_size, max_family_size, lambda));
	
//        for (int i = 0; i<n_fams_this_size; ++i) {
//            cout << "Doing family " << i << endl; 
//            //result.push_back(simulate_families_from_root_size(p_tree, num_trials, root_size, max_family_size, lambda));
//            result.push_back(simulate_families_from_root_size(p_tree, 1, root_size, max_family_size, lambda));
//            //it++;
//        }
        
        it++;
    }

    cout << "The size of vector<vector<trial>> is " << result.size() << endl;
    cout << "When root = 1, the number of simulations should be 5: " << result[0].size() << endl;
    cout << "Then, there should be 7 nodes: " << result[0][0].size() << endl;
//    cout << "I should be able to access species A in map: " << result[0][0][0]->first << endl;
    return result;
}
*/

//void print_simulation(trial &sim, std::ostream& ost)
//{
//	trial::iterator it = sim.begin();
//	for (; it != sim.end(); ++it)
//	{
//		ost << "#" << it->first->get_taxon_name() << endl;
//	}
//	for (it = sim.begin(); it != sim.end(); ++it)
//	{
//		ost << it->second << "\t";
//	}
//	ost << endl;
//
//}

/*
//void print_simulation(std::vector<trial> &sim, std::ostream& ost) {
void print_simulation(std::vector<vector<trial *> > &sim, std::ostream& ost) {
    //Before: trial::iterator it = sim[0].begin();
    trial::iterator names_it = sim[0][0].begin();
    
    // printing header with '#' followed by species names, in the order they will appear below
    for (; names_it != sim[0][0].end(); ++names_it) {
	ost << "#" << names_it->first->get_taxon_name() << endl;
    }
    
    // now printing gene family sizes
//Before:    for (int i = 0; i < sim.size(); ++i) {
//        for (it = sim[i].begin(); it != sim[i].end(); ++it) {
//            ost << it->second << "\t";
//	}
//        
//        ost << endl;
//    }

    for (int i = 0; i < sim.size(); ++i) { // iterating over gene family sizes
        for (int j; j < sim[i].size(); ++j) { // iterating over trial of ith gene family size
            for (trial::iterator jth_trial_it = sim[i][j].begin(); jth_trial_it != sim[i][j].end(); ++jth_trial_it) { // accessing jth trial inside ith gene family
                ost << "porra" << endl;
                ost << jth_trial_it->second << "\t";
            }
            
            ost << endl;
        }
    }
}
*/