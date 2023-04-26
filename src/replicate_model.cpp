#include <vector>
#include <numeric>

#include "replicate_model.h"
#include "gene_transcript.h"
#include "DiffMat.h"

#include "doctest.h"
#include "easylogging++.h"

using namespace Eigen;
using namespace std;

void replicate_model::vlog_vector(std::string replicate_name, std::string taxon_name, VectorXd likelihood_vector) const {
	VLOG(9) << replicate_name << "\t" << taxon_name << " >>> " << likelihood_vector.transpose().format(_fmt) << endl;
	if( replicate_name == "SUM") VLOG(9) << endl;
}

void replicate_model::apply(const clade* node, const gene_transcript& gene_transcript, int upper_bound, VectorXd& result) const
{
	bool rep_found = false;
	result.setZero();
	auto t = node->get_taxon_name();
	// auto gene_id = gene_transript._id # can't access this because it's private
	std::string gene_id = gene_transcript.get_id();
	VLOG(9) << "gene_id = " << gene_id << endl;
	for (auto p : _replicates)
	{
		auto replicate = p.first;
		auto taxon = p.second;
		if (taxon == node->get_taxon_name())
		{
			rep_found = true;
			VectorXd values(result.size());

			double expression_value = gene_transcript.get_expression_value(replicate);
			VectorPos_bounds(expression_value, boundaries(0, upper_bound), values);
			
			// original code to be compatible with main branch:
			// this->vlog_vector(replicate, taxon, rep_values);
			// values += rep_values;
			
			// code to work with replicate_model branch
			this->vlog_vector(replicate, taxon, values);
			result += values;
		}
			
	}
	// original code for main branch: (within additional conditional block)
	// this->vlog_vector("sum", t, values);
	// result.set(values);

	// code to work with replicate_model branch
	this->vlog_vector("SUM", t, result);

	if (!rep_found)
	{
		VectorPos_bounds(gene_transcript.get_expression_value(node->get_taxon_name()), boundaries(0, upper_bound), result);
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

	VectorXd result(10);
	model.apply(p_tree->find_descendant("A"), g, 10, result);

	ostringstream ost;
	const IOFormat fmt(2, DontAlignCols, " ", " ", "", "", "", "");
	ost << result.format(fmt) << endl;
	CHECK_STREAM_CONTAINS(ost, "0 0 0 0 0.45 0.45 0.63 0.27 0 0");
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
