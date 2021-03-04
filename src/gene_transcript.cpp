#include <algorithm>
#include <set>
#include <stdexcept>
#include <sstream>
#include <memory>

#include "doctest.h"

#include "gene_transcript.h"
#include "clade.h"
#include "io.h"

using namespace std;

// these functions are intended to work with maps (key, value pairs)
template <typename T, typename U>
bool max_key(const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
	return p1.first < p2.first;
}

template <typename T, typename U>
bool max_value(const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
	return p1.second < p2.second;
}

double gene_transcript::get_max_expression_value() const {
    // Max family size can only be found if there is data inside the object in the first place
    double max_family_size = 0.0;
    if (!_species_size_map.empty()) {
        max_family_size = (*max_element(_species_size_map.begin(), _species_size_map.end(), max_value<string, double>)).second;
    }
    return max_family_size;
}


double gene_transcript::get_expression_value(std::string species) const {
    // First checks if species data has been entered (i.e., is key in map?)
    if (_species_size_map.find(species) == _species_size_map.end()) {
        throw std::runtime_error(species + " was not found in gene family " + _id);
    }

    return _species_size_map.at(species);
}

//! Return first element of pair
template <typename type1, typename type2>
type1 do_get_species(const std::pair<type1, type2> & p1) {
    return p1.first;
}

//! Return vector of species names
vector<std::string> gene_transcript::get_species() const {
    vector<std::string> species_names(_species_size_map.size());
    transform(_species_size_map.begin(), _species_size_map.end(), species_names.begin(), do_get_species<string, int>); // Transform performs an operation on all elements of a container (here, the operation is the template)

    return species_names;
}

/// returns true if the family exists at the root, according to their parsimony reconstruction.
bool gene_transcript::exists_at_root(const clade *p_tree) const
{
    set<const clade *> exists;
    auto registered = [&exists](const clade *c) {
        return exists.find(c) != exists.end();
    };
    auto existence = [&](const clade *pc) {
        if (pc->is_leaf())
        {
            double sz = get_expression_value(pc->get_taxon_name());
            if (sz > 0)
            {
                auto p = pc;
                do
                {
                    exists.insert(p);
                    //cout << "Registering node " << p->get_taxon_name() << endl;
                    p = p->get_parent();
                } while (p && !registered(p));
            }
        }
    };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), existence);

    bool exists_at_all_children = true;
    auto does_child_exist = [&](const clade *child) { exists_at_all_children &= registered(child); };
    p_tree->apply_to_descendants(does_child_exist);

    return exists_at_all_children;
}

double gene_transcript::species_size_differential() const
{
    auto compare = [](const std::pair<string, int>& a, const std::pair<string, int>& b) { return a.second < b.second; };
    double max_species_size = max_element(_species_size_map.begin(), _species_size_map.end(), compare)->second;
    double min_species_size = min_element(_species_size_map.begin(), _species_size_map.end(), compare)->second;
    return max_species_size - min_species_size;
}


TEST_CASE("exists_at_root returns false if not all children exist")
{
    unique_ptr<clade> p_tree(parse_newick(" ((((cat:68.710687,horse:68.710687):4.566771,cow:73.277458):20.722542,(((((chimp:4.444178,human:4.444178):6.682660,orang:11.126837):2.285866,gibbon:13.412704):7.211528,(macaque:4.567239,baboon:4.567239):16.056993):16.060691,marmoset:36.684923):57.315077)mancat:38.738115,(rat:36.302467,mouse:36.302467):96.435648)"));

    istringstream ist(
        "Desc\tFamily ID\tcat\thorse\tcow\tchimp\thuman\torang\tgibbon\tmacaque\tbaboon\tmarmoset\trat\tmouse\n"
        "(null)\t1\t0\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0\n"
        "(null)\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\n");

    vector<gene_transcript> families;
    read_gene_families(ist, p_tree.get(), families);
    CHECK_FALSE(families[0].exists_at_root(p_tree.get()));
    CHECK_FALSE(families[1].exists_at_root(p_tree.get()));
}

TEST_CASE("exists_at_root_returns_true_if_all_children_exist")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    gene_transcript family;
    family.set_expression_value("A", 3);
    family.set_expression_value("B", 6);

    CHECK(family.exists_at_root(p_tree.get()));
}

TEST_CASE("exists_at_root handles values close to 0")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    gene_transcript family;
    family.set_expression_value("A", 0.2);
    family.set_expression_value("B", 0.3);

    CHECK(family.exists_at_root(p_tree.get()));
}

