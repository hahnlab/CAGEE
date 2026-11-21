#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

#include "doctest.h"
#include "easylogging++.h"

#include "gene_transcript.h"

using namespace std;

// Helper: compute ranks for a vector
std::vector<double> compute_ranks(const std::vector<double>& values) {
    std::vector<size_t> idx(values.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::vector<double> ranks(values.size());
    std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) { return values[i] < values[j]; });

    for (size_t first = 0; first < idx.size();) {
        size_t last = first + 1;
        while (last < idx.size() && values[idx[last]] == values[idx[first]]) {
            ++last;
        }

        const double average_rank = (static_cast<double>(first + 1) + static_cast<double>(last)) / 2.0;
        for (size_t i = first; i < last; ++i) {
            ranks[idx[i]] = average_rank;
        }
        first = last;
    }
    return ranks;
}

// Main function: returns map from species to Spearman correlation
map<pair<string, string>, double> spearman_correlation_by_species(const std::vector<gene_transcript>& transcripts) {
    std::map<std::string, std::vector<double>> species_to_values;

    if (transcripts.size() < 2) {
        throw std::runtime_error("At least two transcripts are required to compute Spearman correlation.");
    }
    // Gather expression values for each species
    for (const auto& transcript : transcripts) {
        for (const auto& species : transcript.get_species()) {
            species_to_values[species].push_back(transcript.get_expression_value(species));
        }
    }

    std::vector<std::pair<std::string, std::string>> species_pairs;
    for (auto it1 = species_to_values.begin(); it1 != species_to_values.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != species_to_values.end(); ++it2) {
            species_pairs.emplace_back(it1->first, it2->first);
        }
    }

    map<pair<string, string>, double> spearman_correlations;
    for (auto species_pair : species_pairs) {
        const auto& sp1 = species_pair.first;
        const auto& sp2 = species_pair.second;

        auto ranks1 = compute_ranks(species_to_values[sp1]);
        auto ranks2 = compute_ranks(species_to_values[sp2]);
        if (ranks1.size() != ranks2.size()) continue;
        const double mean1 = accumulate(ranks1.begin(), ranks1.end(), 0.0) / ranks1.size();
        const double mean2 = accumulate(ranks2.begin(), ranks2.end(), 0.0) / ranks2.size();
        double covariance = 0.0;
        double variance1 = 0.0;
        double variance2 = 0.0;
        for (size_t i = 0; i < ranks1.size(); ++i) {
            const double centered1 = ranks1[i] - mean1;
            const double centered2 = ranks2[i] - mean2;
            covariance += centered1 * centered2;
            variance1 += centered1 * centered1;
            variance2 += centered2 * centered2;
        }
        const double denominator = sqrt(variance1 * variance2);
        const double spearman_corr = denominator == 0.0
            ? numeric_limits<double>::quiet_NaN()
            : covariance / denominator;
        
        spearman_correlations[species_pair] = spearman_corr;
    }

    return spearman_correlations;
}


TEST_CASE("spearman correlation_by_species computes correct correlations")
{
    // Work the example of spearman correlation here: 
    // https://www.statisticshowto.com/probability-and-statistics/correlation-coefficient-formula/spearman-rank-correlation-definition-calculate/
    std::vector<double> sp1 = {35, 23, 47, 17, 10, 43, 9, 6, 28};
    std::vector<double> sp2 = {30, 33, 45, 23, 8, 49, 12, 4, 31};

    std::vector<gene_transcript> transcripts(sp1.size());
    for (size_t i = 0; i < sp1.size(); i++)
    {
        transcripts[i].set_expression_value("species1", sp1[i]);
        transcripts[i].set_expression_value("species2", sp2[i]);    
    }
     
    auto result = spearman_correlation_by_species(transcripts);
    CHECK_EQ(1, result.size());
    CHECK_EQ(doctest::Approx(0.9), result[make_pair("species1","species2")]); 
}

TEST_CASE("spearman correlation_by_species handles tied values")
{
    std::vector<double> sp1 = {1, 1, 2, 3};
    std::vector<double> sp2 = {4, 3, 2, 1};

    std::vector<gene_transcript> transcripts(sp1.size());
    for (size_t i = 0; i < sp1.size(); i++)
    {
        transcripts[i].set_expression_value("species1", sp1[i]);
        transcripts[i].set_expression_value("species2", sp2[i]);
    }

    auto result = spearman_correlation_by_species(transcripts);
    CHECK_EQ(doctest::Approx(-0.9486832980505138), result[make_pair("species1", "species2")]);
}

TEST_CASE("spearman_correlation_by_species throws when less than two transcripts")
{
    string msg("At least two transcripts are required to compute Spearman correlation.");
    std::vector<gene_transcript> transcripts;
    CHECK_THROWS_WITH(spearman_correlation_by_species(transcripts), msg.c_str());

    transcripts.push_back(gene_transcript("transcript1", "", ""));
    CHECK_THROWS_WITH(spearman_correlation_by_species(transcripts), msg.c_str());
}

TEST_CASE("compute_ranks")
{
    std::vector<double> values = {10, 20, 30, 25, 15};
    auto ranks = compute_ranks(values);
    CHECK_EQ(ranks.size(), values.size());
    CHECK_EQ(ranks[0], 1); 
    CHECK_EQ(ranks[1], 3); 
    CHECK_EQ(ranks[2], 5); 
    CHECK_EQ(ranks[3], 4); 
    CHECK_EQ(ranks[4], 2); 
}

TEST_CASE("compute_ranks assigns average ranks to ties")
{
    std::vector<double> values = {10, 20, 20, 40};
    auto ranks = compute_ranks(values);
    CHECK_EQ(ranks[0], 1);
    CHECK_EQ(ranks[1], 2.5);
    CHECK_EQ(ranks[2], 2.5);
    CHECK_EQ(ranks[3], 4);
}

