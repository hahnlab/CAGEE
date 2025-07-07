#ifndef gene_transcript_H
#define gene_transcript_H

#include <string>
#include <map>
#include <vector>
#include <stdexcept>

#include "clade.h"

struct ci_less
{
    // case-independent (ci) compare_less binary function
    struct nocase_compare
    {
        bool operator() (const unsigned char& c1, const unsigned char& c2) const {
            return tolower(c1) < tolower(c2);
        }
    };
    bool operator() (const std::string& s1, const std::string& s2) const {
        return std::lexicographical_compare
        (s1.begin(), s1.end(),   // source range
            s2.begin(), s2.end(),   // dest range
            nocase_compare());  // comparison
    }
};

class missing_expression_value : public std::runtime_error
{
public:
    missing_expression_value(std::string transcript_id, std::string species) : std::runtime_error(species + " was not found in transcript " + transcript_id)
    {

    }


};

class gene_transcript {
private:
    std::string _id; 
    std::string _desc; 
    std::string _tissue;
    std::map<std::string, double, ci_less> _species_size_map; 

public:
    gene_transcript() { }
    gene_transcript(std::string id, std::string desc, std::string tissue) : _id(id), _desc(desc), _tissue(tissue)
    { 
    }

    gene_transcript(const gene_transcript& other) : _id(other._id), _desc(other._desc), _tissue(other._tissue), _species_size_map(other._species_size_map)
    {
    }
    gene_transcript(gene_transcript&& other) noexcept
    {
        *this = std::move(other);
    }

    void set_expression_value(std::string species, double val) {
        _species_size_map[species] = val;
    }

    std::vector<std::string> get_species() const;

    /// @brief Returns the minimum and maximum expression values for this transcript.
    /// @details If no values are set, returns (0.0, 0.0).
    std::pair<double, double> get_expression_boundaries() const;

    std::string id() const { return _id; }
    std::string tissue() const { return _tissue; }
    std::string description() const { return _desc; }

    double get_expression_value(std::string species) const;

    bool species_size_match(const gene_transcript& other) const
    {
        return _species_size_map == other._species_size_map;
    }

    bool exists_at_root(const clade *p_tree) const;

    //! Returns largest species size minus smallest species size
    double species_size_differential() const;

    // move assignment operator
    gene_transcript& operator=(gene_transcript&& other) noexcept
    {
        _species_size_map = std::move(other._species_size_map);
        _id = std::move(other._id);
        _desc = std::move(other._desc);
        _tissue = std::move(other._tissue);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& ost, const gene_transcript& transcript);

    static void remove_ungrouped_transcripts(const std::vector<std::string>& sample_groups, std::vector<gene_transcript>& transcripts);
};

typedef const std::vector<gene_transcript> transcript_vector;

#endif
