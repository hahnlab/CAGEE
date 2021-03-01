#ifndef gene_transcript_H
#define gene_transcript_H

#include <string>
#include <map>
#include <vector>

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


class gene_transcript {
private:
    std::string _id; //!< Gene family ID
    std::string _desc; //!< Gene family description

    std::map<std::string, double, ci_less> _species_size_map; //!< Map that stores each species gene family count: {sp1_name:count1, ...}

public:
    gene_transcript() { }
    gene_transcript(const gene_transcript& other) : _id(other._id), _desc(other._desc), _species_size_map(other._species_size_map)
    {
    }
    gene_transcript(gene_transcript&& other)
    {
        *this = std::move(other);
    }

    void set_desc(std::string desc) { _desc = desc; }

    void set_id(std::string id) { _id = id; }

    void set_species_size(std::string species, double gene_count) {
        _species_size_map[species] = gene_count;
    }

    std::vector<std::string> get_species() const;

    int get_max_size() const;

    std::string id() const { return _id; }

    double get_species_size(std::string species) const;

    //! Returns true if every species size for both gene families are identical
    bool species_size_match(const gene_transcript& other) const
    {
        return _species_size_map == other._species_size_map;
    }

    /// returns true if the family exists at the root of the given tree, according to their parsimony reconstruction.
    bool exists_at_root(const clade *p_tree) const;

    //! Returns largest species size minus smallest species size
    double species_size_differential() const;

    // move assignment operator
    gene_transcript& operator=(gene_transcript&& other)
    {
        _species_size_map = std::move(other._species_size_map);
        _id = std::move(other._id);
        _desc = std::move(other._desc);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& ost, const gene_transcript& family);
};
#endif
