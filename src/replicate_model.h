#ifndef REPLICATE_MODEL_H
#define REPLICATE_MODEL_H

#include <map>
#include <string>

class clade;
class gene_transcript;
class optional_probabilities;

using boundaries = std::pair<double, double>;


class replicate_model
{
public:
    // replicate id to species id
    std::map<std::string, std::string> _replicates;

    void apply(const clade* node, const gene_transcript& gene_transcript, boundaries bounds, optional_probabilities& result) const;

    void verify_replicates(const clade* p_tree, const gene_transcript& t) const;

    double get_average_expression_value(const gene_transcript& t, std::string taxon) const;
};


















#endif