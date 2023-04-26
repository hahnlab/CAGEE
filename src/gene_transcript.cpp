#include <algorithm>
#include <set>
#include <stdexcept>
#include <sstream>
#include <memory>

#include <boost/algorithm/string.hpp>

#include "doctest.h"
#include "easylogging++.h"

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

std::string gene_transcript::get_id() const {
    return this->_id;
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
    auto it = _species_size_map.find(species);
    if (it == _species_size_map.end()) {
        throw std::runtime_error(species + " was not found in transcript " + _id);
    }
    return it->second;
}

//! Return vector of species names
vector<std::string> gene_transcript::get_species() const {
    vector<std::string> species_names(_species_size_map.size());
    transform(_species_size_map.begin(), _species_size_map.end(), species_names.begin(), [](const std::pair<string, double>& p1) { return p1.first; });

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
    auto compare = [](const std::pair<string, double>& a, const std::pair<string, double>& b) { return a.second < b.second; };
    double max_species_size = max_element(_species_size_map.begin(), _species_size_map.end(), compare)->second;
    double min_species_size = min_element(_species_size_map.begin(), _species_size_map.end(), compare)->second;
    return max_species_size - min_species_size;
}

vector<string> expand_sample_groups(const vector<string>& sg)
{
    vector<string> result;
    for (auto s : sg)
    {
        vector<string> groups;
        boost::split(groups, s, boost::is_any_of(","));
        copy(groups.begin(), groups.end(), back_inserter(result));
    }
    return result;
}

void gene_transcript::remove_ungrouped_transcripts(const std::vector<std::string>& sample_groups, std::vector<gene_transcript>& transcripts)
{
    if (sample_groups.empty()) return;

    auto all_groups = expand_sample_groups(sample_groups);
    auto rem = std::remove_if(transcripts.begin(), transcripts.end(), [all_groups](const gene_transcript& t) {
        return find(all_groups.begin(), all_groups.end(), t.tissue()) == all_groups.end();
        });

    auto fmsize = transcripts.size();
    transcripts.erase(rem, transcripts.end());
    if (fmsize != transcripts.size())
    {
        LOG(WARNING) << "Some transcripts were not matched to a sample group";
    }
}

TEST_CASE("exists_at_root returns false if not all children exist")
{
    unique_ptr<clade> p_tree(parse_newick(" ((((cat:68.710687,horse:68.710687):4.566771,cow:73.277458):20.722542,(((((chimp:4.444178,human:4.444178):6.682660,orang:11.126837):2.285866,gibbon:13.412704):7.211528,(macaque:4.567239,baboon:4.567239):16.056993):16.060691,marmoset:36.684923):57.315077)mancat:38.738115,(rat:36.302467,mouse:36.302467):96.435648)"));

    istringstream ist(
        "Desc\tFamily ID\tcat\thorse\tcow\tchimp\thuman\torang\tgibbon\tmacaque\tbaboon\tmarmoset\trat\tmouse\n"
        "(null)\t1\t0\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0\n"
        "(null)\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\n");

    vector<gene_transcript> gt;
    read_gene_transcripts(ist, p_tree.get(), gt);
    CHECK_FALSE(gt[0].exists_at_root(p_tree.get()));
    CHECK_FALSE(gt[1].exists_at_root(p_tree.get()));
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

    gene_transcript gt;
    gt.set_expression_value("A", 0.2);
    gt.set_expression_value("B", 0.3);

    CHECK(gt.exists_at_root(p_tree.get()));
}

TEST_CASE("remove_ungrouped_transcripts removes nothing if no groups")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    istringstream ist(
        "Desc\tFamily ID\tSAMPLETYPE\tA\tB\n"
        "(null)\tt1\theart\t0\t0\n"
        "(null)\tt2\tlungs\t0\t0\n");

    vector<gene_transcript> transcripts;

    read_gene_transcripts(ist, p_tree.get(), transcripts);

    CHECK_EQ(2, transcripts.size());
    gene_transcript::remove_ungrouped_transcripts(vector<string>(), transcripts);

    CHECK_EQ(2, transcripts.size());
}

TEST_CASE("remove_ungrouped_transcripts removes specified single group")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    istringstream ist(
        "Desc\tFamily ID\tSAMPLETYPE\tA\tB\n"
        "(null)\tt1\theart\t0\t0\n"
        "(null)\tt2\tlungs\t0\t0\n");

    vector<gene_transcript> transcripts;

    read_gene_transcripts(ist, p_tree.get(), transcripts);

    CHECK_EQ(2, transcripts.size());
    gene_transcript::remove_ungrouped_transcripts(vector<string>({"lungs"}), transcripts);

    CHECK_EQ(1, transcripts.size());
}

TEST_CASE("remove_ungrouped_transcripts is aware of comma-specified groups")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    istringstream ist(
        "Desc\tFamily ID\tSAMPLETYPE\tA\tB\n"
        "(null)\tt1\tbrain\t0\t0\n"
        "(null)\tt2\tlungs\t0\t0\n");

    vector<gene_transcript> transcripts;

    read_gene_transcripts(ist, p_tree.get(), transcripts);

    CHECK_EQ(2, transcripts.size());
    gene_transcript::remove_ungrouped_transcripts(vector<string>({ "heart,lungs" }), transcripts);

    CHECK_EQ(1, transcripts.size());
}

TEST_CASE("expand_sample_groups")
{
    vector<string> sample_groups({ "heart,lungs", "brain", "kidneys,liver,spleen" });
    auto actual = expand_sample_groups(sample_groups);
    vector<string> expected({ "heart" , "lungs", "brain", "kidneys", "liver", "spleen" });
    CHECK_EQ(expected, actual);
}