#ifndef SPEARMAN_H
#define SPEARMAN_H

#include <map>
#include <vector>
#include <string>

class gene_transcript;

std::map<std::pair<std::string, std::string>, double> spearman_correlation_by_species(const std::vector<gene_transcript>& transcripts);


#endif