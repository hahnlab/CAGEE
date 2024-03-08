#include <vector>
#include <numeric>

#include "replicate_model.h"
#include "gene_transcript.h"
#include "DiffMat.h"
#include "inference_pruner.h"

#include "doctest.h"
#include "easylogging++.h"

using namespace Eigen;
using namespace std;

void replicate_model::apply(const clade* node, const gene_transcript& gene_transcript, boundaries bounds, optional_probabilities& result) const
{
	int Npts = result.capacity();

	bool rep_found = false;
	VectorXd values(Npts);
	values.setZero();
	auto t = node->get_taxon_name();
	for (auto p : _replicates)
	{
		auto replicate = p.first;
		auto taxon = p.second;
		if (taxon == node->get_taxon_name())
		{
			rep_found = true;
			VectorXd rep_values(Npts);

			double expression_value = gene_transcript.get_expression_value(replicate);
			VectorPos_bounds(expression_value, bounds, rep_values);
			values += rep_values;
		}
	}
	if (rep_found)
	{
		result.set(values);
	}
	else
	{
		result.initialize(gene_transcript.get_expression_value(node->get_taxon_name()), bounds);
	}
}

void replicate_model::verify_replicates(const clade* p_tree, const gene_transcript& t) const
{
	set<string> taxa;
	p_tree->apply_prefix_order([&taxa](const clade* c) { if (c->is_leaf()) taxa.insert(c->get_taxon_name()); });
	for (auto p : _replicates)
	{
		if (find(taxa.begin(), taxa.end(), p.second) != taxa.end())
		{
			taxa.insert(p.first);
		}
	}
	for (auto p : _replicates)
	{
		if (find(taxa.begin(), taxa.end(), p.second) != taxa.end())
		{
			taxa.erase(p.second);
		}
	}
	for (auto taxon : taxa)
	{
		t.get_expression_value(taxon);
	}
}

double replicate_model::get_average_expression_value(const gene_transcript& t, string ttaxon) const
{
	vector<double> values;
	for (auto p : _replicates)
	{
		auto replicate = p.first;
		auto taxon = p.second;
		if (taxon == ttaxon)
		{
			values.push_back(t.get_expression_value(replicate));
		}
	}

	if (values.empty())
		return t.get_expression_value(ttaxon);
	else
		return accumulate(values.begin(), values.end(), 0) / double(values.size());
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("replicate_model replaces leaf values by adding together the replicates")
{
	gene_transcript g;
	g.set_expression_value("A-1", 5.0);
	g.set_expression_value("A-2", 7.0);

	replicate_model model;
	model._replicates["A-1"] = "A";
	model._replicates["A-2"] = "A";

	unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

	VectorXd v(10);
	optional_probabilities result;
	result.reserve(v);
	model.apply(p_tree->find_descendant("A"), g, boundaries(0,10), result);

	ostringstream ost;
	const IOFormat fmt(2, DontAlignCols, " ", " ", "", "", "", "");
	ost << result.probabilities().format(fmt) << endl;
	CHECK_STREAM_CONTAINS(ost, "0 0 0 0 0.5 0.5 0.7 0.3 0 0");
}

TEST_CASE("verify_replicates ")
{
	gene_transcript g;
	g.set_expression_value("A-1", 5.0);
	g.set_expression_value("A-2", 7.0);
	g.set_expression_value("B", 11.0);

	replicate_model model;
	model._replicates["A-1"] = "A";
	model._replicates["A-2"] = "A";

	unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

	model.verify_replicates(p_tree.get(), g);
}

TEST_CASE("verify_replicates throws if a node has no replicates")
{
	gene_transcript g;
	g.set_expression_value("A-1", 5.0);
	g.set_expression_value("A-2", 7.0);
	g.set_expression_value("B-1", 11.0);
	g.set_expression_value("B-2", 11.0);

	replicate_model model;
	model._replicates["A-1"] = "A";
	model._replicates["A-2"] = "A";
	model._replicates["B-1"] = "A";
	model._replicates["B-2"] = "A";

	unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

	CHECK_THROWS_WITH(model.verify_replicates(p_tree.get(), g), "B was not found in transcript ");
}

TEST_CASE("verify_replicates throws if a node has no replicates")
{
	replicate_model model;
	model._replicates["A-1"] = "A";
	model._replicates["A-2"] = "A";
	model._replicates["B-1"] = "B";
	model._replicates["B-2"] = "B";

	gene_transcript g;
	g.set_expression_value("A-1", 5.0);
	g.set_expression_value("A-2", 7.0);
	g.set_expression_value("B-1", 11.0);
	g.set_expression_value("B-2", 19.0);

	CHECK_EQ(6.0, model.get_average_expression_value(g, "A"));
	CHECK_EQ(15.0, model.get_average_expression_value(g, "B"));
}
