#ifndef FITCH_MARGOLIASH_H
#define FITCH_MARGOLIASH_H

#include <map>
#include <string>
#include <vector>

class clade;

class fitch_margoliash {
    // The fitch_margoliash algorithm builds a tree based on the distance matrix
public:
    using distance_matrix = std::map<std::pair<std::string, std::string>, double>;

    fitch_margoliash(distance_matrix distance_matrix, const std::vector<std::string>& species_names);

    double average_distance(std::string s) const;

    clade *build_tree();

    void compute_branch_lengths(std::string sp1, std::string sp2, clade *new_cluster);

    double distance(std::string sp1, std::string sp2) const;

    std::pair<std::string, std::string> closest_clusters();

private:
    distance_matrix _distance_matrix;
    std::map<std::string, clade*> _clusters;
};


#endif