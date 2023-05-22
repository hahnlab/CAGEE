#include <algorithm>
#include <iomanip>
#include <sstream>
#include <memory>

#include <boost/algorithm/string.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "sigma.h"
#include "matrix_cache.h"
#include "clade.h"
#include "gene_transcript.h"

using namespace std;

class ss_resolver
{
public:
    virtual int get_value_index(const clade* c, const gene_transcript& t) const = 0;
    virtual size_t count() const = 0;
    virtual ss_resolver* clone() const = 0;
    virtual ~ss_resolver() {}
};

class ss_resolver_uniform : public ss_resolver
{
public:
    virtual int get_value_index(const clade* c, const gene_transcript& t) const { return 0; }

    virtual ss_resolver* clone() const { return new ss_resolver_uniform(); }

    virtual size_t count() const { return 1; }
};

class ss_resolver_lineage_specific : public ss_resolver
{
    std::map<std::string, int> _node_name_to_sigma_index;
public:
    ss_resolver_lineage_specific(std::map<std::string, int> nodename_index_map) :
        _node_name_to_sigma_index(nodename_index_map)
    {

    }

    virtual int get_value_index(const clade* c, const gene_transcript& t) const
    { 
        return _node_name_to_sigma_index.at(c->get_taxon_name());
    }

    virtual ss_resolver* clone() const { return new ss_resolver_lineage_specific(_node_name_to_sigma_index); }

    virtual size_t count() const 
    {
        set<size_t> s;
        for (auto a : _node_name_to_sigma_index)
            s.insert(a.second);
        return s.size(); 
    }
};

class ss_resolver_sample_specific : public ss_resolver
{
    std::map<std::string, int> _sample_to_sigma_index;
public:
    ss_resolver_sample_specific(std::map<std::string, int> nodename_index_map) :
        _sample_to_sigma_index(nodename_index_map)
    {

    }

    virtual int get_value_index(const clade* c, const gene_transcript& t) const
    {
        return _sample_to_sigma_index.at(t.tissue());
    }

    virtual ss_resolver* clone() const { return new ss_resolver_sample_specific(_sample_to_sigma_index); }

    virtual size_t count() const 
    { 
        set<size_t> s;
        for (auto a : _sample_to_sigma_index)
            s.insert(a.second);
        return s.size();
    }
};

class ss_resolver_composite : public ss_resolver
{
    using key = pair<size_t, size_t>;
    ss_resolver* _r1;
    ss_resolver* _r2;
    map<key,size_t> values;
public:
    ss_resolver_composite(ss_resolver *r1, ss_resolver* r2) :
        _r1(r1), _r2(r2)
    {
        int n = 0;
        for (size_t i = 0; i < r1->count(); ++i)
            for (size_t j = 0; j < r2->count(); ++j)
                values[key(i, j)] = n++;
    }

    ~ss_resolver_composite()
    {
        delete _r1;
        delete _r2;
    }
    virtual int get_value_index(const clade* c, const gene_transcript& t) const
    {
        return values.at(key(_r1->get_value_index(c, t), _r2->get_value_index(c, t)));
    }

    virtual ss_resolver* clone() const { return new ss_resolver_composite(_r1->clone(), _r2->clone()); }

    virtual size_t count() const { return values.size(); }
};

sigma_squared::sigma_squared(double val) : _p_resolver(new ss_resolver_uniform())
{
    _values.push_back(val);
}

sigma_squared::sigma_squared(ss_resolver* p_resolver, size_t sz) : 
    _p_resolver(p_resolver), _values(sz)
{
}

sigma_squared::sigma_squared(ss_resolver* p_resolver, const std::vector<double>& v) :
    _p_resolver(p_resolver), _values(v)
{

}

sigma_squared::~sigma_squared()
{
    delete _p_resolver;
}

sigma_squared* sigma_squared::multiply(double factor) const
{
    auto npi = _values;

    for (auto& i : npi)
        i *= factor;

    return new sigma_squared(_p_resolver->clone(), npi);
}

int sigma_squared::count() const {
    return _values.size();
}

std::vector<double> sigma_squared::get_values() const {
    return _values;
}

sigma_squared* sigma_squared::clone() const {
    return new sigma_squared(_p_resolver->clone(), _values);
}

void sigma_squared::update(const double* values)
{
    std::copy(values, values + _values.size(), _values.begin());
}

std::string sigma_squared::to_string() const
{
    ostringstream ost;
    ost << setprecision(14);
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

double sigma_squared::get_named_value(const clade* c, const gene_transcript& t) const {

    return _values[_p_resolver->get_value_index(c, t)];
}

map<string, int> get_sigma_index_map(const std::vector<string>& sample_groups)
{
    map<string, int> sample_to_sigma_index;
    for (size_t i = 0; i < sample_groups.size(); ++i)
    {
        vector<string> groups;
        boost::split(groups, sample_groups[i], boost::is_any_of(","));
        for (auto g : groups)
            sample_to_sigma_index[g] = i;
    }
    return sample_to_sigma_index;
}

sigma_squared* sigma_squared::create(clade* p_sigma_tree, const std::vector<string>& sample_groups)
{
    sigma_squared* result = nullptr;

    vector<ss_resolver*> resolvers;
    if (p_sigma_tree)
    {
        std::set<int> unique_sigmas;
        auto fn = [&unique_sigmas](const clade* p_node) { unique_sigmas.insert(p_node->get_sigma_index()); };
        p_sigma_tree->apply_prefix_order(fn);
        resolvers.push_back(new ss_resolver_lineage_specific(p_sigma_tree->get_sigma_index_map()));
    }

    if (!sample_groups.empty())
    {
       resolvers.push_back(new ss_resolver_sample_specific(get_sigma_index_map(sample_groups)));
    }

    switch (resolvers.size())
    {
    case 0:
        result = new sigma_squared(0.0);
        break;
    case 1:
        result = new sigma_squared(resolvers[0], resolvers[0]->count());
        break;
    case 2:
    {
        auto x = new ss_resolver_composite(resolvers[0], resolvers[1]);
        result = new sigma_squared(x, resolvers[0]->count() * resolvers[1]->count());
    }
        break;
    default:
        throw std::domain_error("Cannot handle more than two resolvers");
    }

    LOG(INFO) << "Searching for " << result->count() << " sigmas" << endl;
    return result;
}

sigma_squared* sigma_squared::create(clade* p_sigma_tree, const std::vector<double>& values)
{
    return new sigma_squared(new ss_resolver_lineage_specific(p_sigma_tree->get_sigma_index_map()), values);
}

TEST_CASE("lineage_specific sigma returns correct values")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    map<string, int> key;
    key["A"] = 5;
    key["B"] = 3;
    auto x = new ss_resolver_lineage_specific(key);

    sigma_squared ml(x, { .03, .05, .07, .011, .013, .017 });

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
    auto x = new ss_resolver_sample_specific(key);

    sigma_squared ml(x, { .03, .05, .07, .011, .013, .017 });
    gene_transcript t("C", "", "heart");
    gene_transcript t2("D", "", "lungs");
    CHECK_EQ(.07, ml.get_named_value(p_tree->find_descendant("A"), t));
    CHECK_EQ(.013, ml.get_named_value(p_tree->find_descendant("B"), t2));
}

TEST_CASE("is_valid returns false if any value is negative")
{
    map<string, int> key;
    sigma_squared ml(new ss_resolver_sample_specific(key), { .03, .05, .07, .011, .013, .017 });

    CHECK(ml.is_valid());

    sigma_squared m2(new ss_resolver_sample_specific(key), { .03, .05, .07, .011, -.013, .017 });
    CHECK_FALSE(m2.is_valid());
}

TEST_CASE("update")
{
    map<string, int> key;
    auto x = new ss_resolver_sample_specific(key);
    sigma_squared ml(x, { .03, .05, .07 });
    double newvalues[3] = { .11, .15, .17 };
    ml.update(newvalues);
    CHECK_EQ(vector<double>({ .11, .15, .17 }), ml.get_values());
}

TEST_CASE("sigma_squared::create returns single sigma if no arguments")
{
    unique_ptr<sigma_squared> sig(sigma_squared::create(nullptr, vector<string>()));
    CHECK_EQ(1, sig->get_values().size());
}

TEST_CASE("sigma_squared::create returns lineage-specific with a sigma tree")
{
    string s = "((((cat:1,horse:1):1,cow:1):1,(((((chimp:2,human:2):2,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:1,mouse:1):1);";
    unique_ptr<clade> p_tree(parse_newick(s, true));

    unique_ptr<sigma_squared> sig(sigma_squared::create(p_tree.get(), vector<string>()));
    REQUIRE_EQ(2, sig->get_values().size());

    double values[2] = { 7,14 };
    sig->update(values);
    CHECK_EQ(7, sig->get_named_value(p_tree->find_descendant("cow"), gene_transcript()));
    CHECK_EQ(14, sig->get_named_value(p_tree->find_descendant("chimp"), gene_transcript()));
}

TEST_CASE("sigma_squared::create returns sample-specific with sample groups")
{
    unique_ptr<sigma_squared> sig(sigma_squared::create(nullptr, vector<string>({ "heart", "lungs" })));
    REQUIRE_EQ(2, sig->get_values().size());

    double values[2] = { 7,14 };
    sig->update(values);
    CHECK_EQ(7, sig->get_named_value(nullptr, gene_transcript("A", "", "heart")));
    CHECK_EQ(14, sig->get_named_value(nullptr, gene_transcript("B", "", "lungs")));
}

TEST_CASE("get_sigma_index_map creates map")
{
    auto m = get_sigma_index_map(vector<string>({ "heart,lungs", "brain", "kidneys,liver,spleen" }));

    REQUIRE_EQ(6, m.size());

    CHECK_EQ(0, m["heart"]);
    CHECK_EQ(0, m["lungs"]);
    CHECK_EQ(1, m["brain"]);
    CHECK_EQ(2, m["kidneys"]);
    CHECK_EQ(2, m["liver"]);
    CHECK_EQ(2, m["spleen"]);
}

TEST_CASE("composite resolver count is all possible sigma pairs")
{
    string s = "((((cat:1,horse:1):1,cow:1):1,(((((chimp:2,human:2):2,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:3,mouse:3):3);";
    unique_ptr<clade> p_sigma_tree(parse_newick(s, true));

    auto ss = sigma_squared::create(p_sigma_tree.get(), vector<string>({ "heart", "lungs", "brain" }));
    CHECK_EQ(9, ss->count());

    gene_transcript gt("1", "desc", "heart");
    CHECK_EQ(0, ss->get_named_value(p_sigma_tree->find_descendant("cat"), gt));
}

TEST_CASE("Sample specific resolver only counts unique sample groups")
{
    ss_resolver_sample_specific ss({ {"heart",1}, {"lungs,brain",2} });

    CHECK_EQ(2, ss.count());
}
