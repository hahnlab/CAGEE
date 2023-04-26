#ifndef REPLICATE_MODEL_H
#define REPLICATE_MODEL_H

#include <map>
#include <string>
#include <Eigen/Dense>

class clade;
class gene_transcript;

class replicate_model
{

private:
    Eigen::IOFormat _fmt{Eigen::IOFormat(3, 0, ", ", "", "[ ", " ]", "", "", (char)32)};

public:
    // replicate id to species id
    std::map<std::string, std::string> _replicates;

    void apply(const clade* node, const gene_transcript& gene_transcript, int upper_bound, Eigen::VectorXd& result) const;

    void verify_replicates(const clade* p_tree, const gene_transcript& t) const;

    double get_average_expression_value(const gene_transcript& t, std::string taxon) const;

    void vlog_vector(std::string replicate_name, std::string taxon_name, Eigen::VectorXd likelihood_vector) const;
};


















#endif