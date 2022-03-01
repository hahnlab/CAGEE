#include <algorithm>
#include <iomanip>
#include <sstream>

#include "doctest.h"

#include "sigma.h"
#include "matrix_cache.h"
#include "clade.h"
#include "gene_transcript.h"

using namespace std;

void sigma_squared::update(const double* values)
{
    std::copy(values, values + _values.size(), _values.begin());
}

std::string sigma_squared::to_string() const
{
    ostringstream ost;
    ost << setw(15) << setprecision(14);
    for (size_t i = 0; i < _values.size(); ++i)
    {
        ost << _values[i];
        if (i != _values.size() - 1) ost << ", ";
    }
    return ost.str();
}

bool sigma_squared::is_valid() const
{
    return std::none_of(_values.begin(), _values.end(), [](double d) { return d < 0; });
}

double sigma_squared::get_value_for_clade(const clade *c) const {
    if (count() == 1)
        return _values[0];

    int index = _node_name_to_sigma_index.at(c->get_taxon_name());
    return _values[index];
}

double sigma_squared::get_named_value(const clade* c, const gene_transcript& t) const {
    string name;
    switch (_type)
    {
    case sigma_type::uniform:
        return _values[0];
    case sigma_type::lineage_specific:
        name = c->get_taxon_name();
        break;
    case sigma_type::sample_specific:
        name = t.tissue();
    }
    return _values[_node_name_to_sigma_index.at(name)];
}


TEST_CASE("lineage_specific sigma returns correct values")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    map<string, int> key;
    key["A"] = 5;
    key["B"] = 3;
    sigma_squared ml(key, { .03, .05, .07, .011, .013, .017 }, sigma_type::lineage_specific);
    CHECK_EQ(.017, ml.get_value_for_clade(p_tree->find_descendant("A")));
    CHECK_EQ(.011, ml.get_value_for_clade(p_tree->find_descendant("B")));

    gene_transcript t;
    CHECK_EQ(.017, ml.get_named_value(p_tree->find_descendant("A"), t));
    CHECK_EQ(.011, ml.get_named_value(p_tree->find_descendant("B"), t));

}

TEST_CASE("sample_specific sigma returns correct values")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    map<string, int> key;
    key["heart"] = 2;
    key["lungs"] = 4;
    sigma_squared ml(key, { .03, .05, .07, .011, .013, .017 }, sigma_type::sample_specific);
    gene_transcript t("C", "", "heart");
    gene_transcript t2("D", "", "lungs");
    CHECK_EQ(.07, ml.get_named_value(p_tree->find_descendant("A"), t));
    CHECK_EQ(.013, ml.get_named_value(p_tree->find_descendant("B"), t2));
}

TEST_CASE("is_valid returns false if any value is negative")
{
    map<string, int> key;
    sigma_squared ml(key, { .03, .05, .07, .011, .013, .017 }, sigma_type::sample_specific);

    CHECK(ml.is_valid());

    sigma_squared m2(key, { .03, .05, .07, .011, -.013, .017 }, sigma_type::sample_specific);
    CHECK_FALSE(m2.is_valid());
}

TEST_CASE("update")
{
    map<string, int> key;
    sigma_squared ml(key, { .03, .05, .07 }, sigma_type::sample_specific);
    double newvalues[3] = { .11, .15, .17 };
    ml.update(newvalues);
    CHECK_EQ(vector<double>({ .11, .15, .17 }), ml.get_values());
}