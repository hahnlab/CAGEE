#include <cmath>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <iterator>

#include "doctest.h"
#include "easylogging++.h"

#include "core.h"
#include "reconstruction.h"
#include "sigma.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "DiffMat.h"
#include "newick_ape_loader.h"
#include "proportional_variance.h"
#include "inference_pruner.h"
#include "io.h"

using namespace std;
using namespace Eigen;
namespace pv = proportional_variance;

bool credible_interval_intersects_parent(const clade* node, const gene_transcript& transcript, const reconstruction* r)
{
    if (node->is_root())
        return true;

    auto pci = r->get_internal_node_value(transcript, node->get_parent()).credible_interval;
    if (node->is_leaf())
    {
        auto val = transcript.get_expression_value(node->get_taxon_name());
        return val >= pci.first && val <= pci.second;
    }
    else
    {
        auto ci = r->get_internal_node_value(transcript, node).credible_interval;
        return !(ci.first >= pci.second || pci.first >= ci.second);
    }
}


void reconstruction::print_increases_decreases_by_clade(std::ostream& ost, const clade *p_tree, transcript_vector& gene_transcripts, bool count_all_changes)
{
    clademap<pair<int, int>> increase_decrease_map;

    for (auto& transcript : gene_transcripts)
    {
        cladevector nodes;
        // See 0015-credible_interval_calculation.md
        copy_if(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), back_inserter(nodes), [this, transcript, count_all_changes](const clade* node)
            {
                return count_all_changes || !credible_interval_intersects_parent(node, transcript, this);
            });

        for(auto node : nodes)
        {
            double diff = get_difference_from_parent(transcript, node);
            if (diff > 0)
                increase_decrease_map[node].first++;
            if (diff < 0)
                increase_decrease_map[node].second++;
        }
    }

    ost << "#Taxon_ID\tIncrease\tDecrease\n";
    for (auto& it : increase_decrease_map) {
        ost << clade_index_or_name(it.first) << "\t";
        ost << it.second.first << "\t";
        ost << it.second.second << endl;
    }
}

template<typename T>
void print_family_clade_table(std::ostream& ost, const cladevector& order, transcript_vector& gene_transcripts, const clade* p_tree, std::function<T(const gene_transcript& transcript, const clade *c)> get_family_clade_value)
{
    write_node_ordered<string>(ost, "TranscriptID", order, [](const clade* c) { return clade_index_or_name(c); });
    for (auto& transcript : gene_transcripts)
    {
        write_node_ordered<T>(ost, transcript.id(), order, [&transcript, get_family_clade_value](const clade* c) { return get_family_clade_value(transcript, c); });
    }
}

double node_value(const clade* node, const gene_transcript& transcript, const reconstruction* r)
{
    double val = node->is_leaf() ? transcript.get_expression_value(node->get_taxon_name()) : r->get_internal_node_value(transcript, node).most_likely_value;

    return pv::to_user_space(val);
}


void reconstruction::print_reconstructed_states(std::ostream& ost, transcript_vector& transcripts, const clade* p_tree)
{
    ost << "#nexus\nBEGIN TREES;\n";
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        auto& gene_transcript = transcripts[i];

        function<void(std::ostream& ost, const clade*)> text_func;
        text_func = [gene_transcript, this](std::ostream& ost, const clade* node) {
            ost << clade_index_or_name(node) << "_" << node_value(node, gene_transcript, this);

            if (!node->is_root())
                ost << ':' << node->get_branch_length();
        };

        ost << "  TREE " << gene_transcript.id() << " = ";
        p_tree->write_newick(ost, text_func);

        ost << ';' << endl;
    }
    ost << "\nEND;\n";
    write_nexus_extensions(ost);
    
}

string node_change(const clade* node, const gene_transcript& transcript, const reconstruction* r)
{
    ostringstream ost;
    ost << showpos << r->get_difference_from_parent(transcript, node);
    if (!credible_interval_intersects_parent(node, transcript, r))
    {
        ost << "*";
    }

    return ost.str();

}

string node_credible_interval(const clade* node, const gene_transcript& transcript, const reconstruction* r)
{
    double lower, upper;
    if (node->is_leaf())
    {
        auto val = pv::to_user_space(transcript.get_expression_value(node->get_taxon_name()));
        lower = val;
        upper = val;
    }
    else
    {
        auto val = r->get_internal_node_value(transcript, node);
        lower = pv::to_user_space(val.credible_interval.first);
        upper = pv::to_user_space(val.credible_interval.second);
    }

    ostringstream ost;
    ost << '[' << lower << "-" << upper << ']';
    return ost.str();

}

void reconstruction::write_results(std::string output_prefix, 
    const clade *p_tree, 
    transcript_vector& transcripts, 
    bool count_all_changes)
{
    cladevector order;
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&order](const clade* c)
        {
            order.push_back(c);
        });
    sort(order.begin(), order.end(), [](const clade* a, const clade* b) { return a->get_ape_index() < b->get_ape_index(); });

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing reconstructed states";
    std::ofstream ofst(filename("ancestral_states", output_prefix, "tre"));
    print_reconstructed_states(ofst, transcripts, p_tree);

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing node transcript values";
    std::ofstream counts(filename("ancestral_states", output_prefix, "tab"));
    print_family_clade_table<double>(counts, order, transcripts, p_tree, [this, transcripts](const gene_transcript& transcript, const clade* c) {
        return node_value(c, transcript, this);
        });

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing node changes";
    std::ofstream change(filename("change", output_prefix, "tab"));
    print_family_clade_table<string>(change, order, transcripts, p_tree, [this, &transcripts](const gene_transcript& transcript, const clade* c) {
        return node_change(c, transcript, this);
        });

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing node credible intervals";
    std::ofstream ci(filename("credible_intervals", output_prefix, "tab"));
    print_family_clade_table<string>(ci, order, transcripts, p_tree, [this, &transcripts](const gene_transcript& transcript, const clade* c) {
        return node_credible_interval(c, transcript, this);
        });

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing clade results";
    std::ofstream clade_results(filename("clade_results", output_prefix, "tab"));
    print_increases_decreases_by_clade(clade_results, p_tree, transcripts, count_all_changes);

    print_additional_data(transcripts, output_prefix);
}

double reconstruction::get_difference_from_parent(const gene_transcript & transcript, const clade * node) const
{
    if (node->is_root())
        return 0.0;

    double childval = node->is_leaf() ? transcript.get_expression_value(node->get_taxon_name()) : get_internal_node_value(transcript, node).most_likely_value;

    double parentval = get_internal_node_value(transcript, node->get_parent()).most_likely_value;

    return pv::to_user_space(childval) - pv::to_user_space(parentval);
}

class Reconstruction
{
public:
    gene_transcript fam;
    unique_ptr<clade> p_tree;

    Reconstruction() : fam("Family5", "", "")
    {
        p_tree.reset(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

        fam.set_expression_value("A", 11);
        fam.set_expression_value("B", 2);
        fam.set_expression_value("C", 5);
        fam.set_expression_value("D", 6);
    }
};

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript assigns actual values to leaves")
{
    sigma_squared sig(0.1);
    fam.set_expression_value("A", 3.7);
    fam.set_expression_value("D", 9.4);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_values(), p_tree->get_branch_lengths(), 20);

    inference_pruner tr(&sig, p_tree.get(), &calc);

    auto actual = tr.reconstruct(fam, 20);

    CHECK_EQ(3.7, actual[p_tree->find_descendant("A")].most_likely_value);
    CHECK_EQ(9.4, actual[p_tree->find_descendant("D")].most_likely_value);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript calculates parent node values correctly")
{
    int Npts = 200;
    sigma_squared sig(10.1);

    MatrixXd doubler = MatrixXd::Identity(Npts, Npts) * 2;
    matrix_cache calc;
    for (auto len : p_tree->get_branch_lengths())
    {
        calc.set_matrix(len, 10.1, 200, doubler);
    }
    inference_pruner tr(&sig, p_tree.get(), &calc);

    auto actual = tr.reconstruct(fam, 200);

    CHECK_EQ(2.5, actual[p_tree->find_descendant("AB")].most_likely_value);

}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript returns parent above 0 with 0s at the leafs")
{
    sigma_squared sig(10.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_values(), p_tree->get_branch_lengths(), 20);

    inference_pruner tr(&sig, p_tree.get(), &calc);
    auto actual = tr.reconstruct(fam, 20);
    CHECK_EQ(0.05, actual[p_tree->find_descendant("AB")].most_likely_value);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript returns correct value at root" * doctest::skip(true))
{
    // TODO: Create reasonable test case values
    sigma_squared sig(10.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_values(), p_tree->get_branch_lengths(), 20);

    inference_pruner tr(&sig, p_tree.get(), &calc);
    auto actual = tr.reconstruct(fam, 20);
    CHECK_EQ(50, actual[p_tree.get()].most_likely_value);
}

class mock_reconstruction : public reconstruction
{
    node_reconstruction get_internal_node_value(const gene_transcript& gf, const clade* c) const override
    { 
        node_reconstruction nr;
        if (_reconstructions.find(gf.id()) == _reconstructions.end())
            nr.most_likely_value = 3;
        else
            nr = _reconstructions.at(gf.id()).at(c);

        return nr;
    }
public:
    std::map<std::string, clademap<node_reconstruction>> _reconstructions;

};

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE_FIXTURE(Reconstruction, "node_value returns expression value for leaves")
{
    auto leaf = p_tree->find_descendant("D");
    mock_reconstruction r;

    CHECK_EQ( doctest::Approx(6), node_value(leaf, fam, &r));

}

TEST_CASE_FIXTURE(Reconstruction, "node_value returns reconstruction value for internal nodes")
{
    auto node = p_tree->find_descendant("CD");
    mock_reconstruction r;
    r._reconstructions["Family5"][node].most_likely_value = pv::to_computational_space(7);

    CHECK_EQ(doctest::Approx(7), node_value(node, fam, &r));

}

TEST_CASE("Reconstruction: print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;

    mock_reconstruction bmr;

    bmr.print_increases_decreases_by_clade(empty, p_tree.get(), {}, true);
    CHECK_EQ(string("#Taxon_ID\tIncrease\tDecrease\n"), empty.str());

    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].most_likely_value = 17.6;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].credible_interval.first = 15;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].credible_interval.second = 20;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 22.11);
    gf.set_expression_value("B", 9.3);

    ostringstream ost;
    bmr.print_increases_decreases_by_clade(ost, p_tree.get(), {gf}, true);
    CHECK_STREAM_CONTAINS(ost, "#Taxon_ID\tIncrease\tDecrease");
    CHECK_STREAM_CONTAINS(ost, "A<1>\t1\t0");
    CHECK_STREAM_CONTAINS(ost, "B<2>\t0\t1");
}

TEST_CASE("Reconstruction: print_increases_decreases_by_clade skips clades inside the credible interval if count_all_changes is false")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;

    mock_reconstruction bmr;

    bmr.print_increases_decreases_by_clade(empty, p_tree.get(), {}, false);
    CHECK_EQ(string("#Taxon_ID\tIncrease\tDecrease\n"), empty.str());

    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].most_likely_value = 17.6;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].credible_interval.first = 8.4;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].credible_interval.second = 19.9;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 22.11);
    gf.set_expression_value("B", 9.3);

    ostringstream ost;
    bmr.print_increases_decreases_by_clade(ost, p_tree.get(), { gf }, false);
    CHECK_STREAM_CONTAINS(ost, "#Taxon_ID\tIncrease\tDecrease");
    CHECK_STREAM_CONTAINS(ost, "A<1>\t1\t0");
    CHECK_MESSAGE(ost.str().find("B<2>\t0\t1") == std::string::npos, ost.str());
}

TEST_CASE("Reconstruction: print_increases_decreases_by_clade includes clades inside the credible interval if count_all_changes is true")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    mock_reconstruction bmr;

    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].most_likely_value = 17.6;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].credible_interval.first = 8.4;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")].credible_interval.second = 19.9;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 22.11);
    gf.set_expression_value("B", 9.3);

    ostringstream ost;
    bmr.print_increases_decreases_by_clade(ost, p_tree.get(), { gf }, true);
    CHECK_STREAM_CONTAINS(ost, "#Taxon_ID\tIncrease\tDecrease");
    CHECK_STREAM_CONTAINS(ost, "A<1>\t1\t0");
    CHECK_STREAM_CONTAINS(ost, "B<2>\t0\t1");
}

TEST_CASE_FIXTURE(Reconstruction, "get_difference_from_parent")
{
    mock_reconstruction r;
    r._reconstructions["Family5"][p_tree->find_descendant("AB")].most_likely_value = pv::to_computational_space(17.6);

    fam.set_expression_value("A", pv::to_computational_space(22.11));
    CHECK_EQ(doctest::Approx(4.51), r.get_difference_from_parent(fam, p_tree->find_descendant("A")));
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_index_in_angle_brackets_for_non_leaf")
{
    CHECK_EQ(string("<5>"), clade_index_or_name(p_tree.get()));
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_name_plus_index_in_angle_brackets_for_leaf")
{
    auto a = p_tree->find_descendant("A");
    CHECK_EQ(string("A<1>"), clade_index_or_name(a));
}

TEST_CASE_FIXTURE(Reconstruction, "credible_interval_intersects_parent is correct for leaf nodes")
{
    mock_reconstruction r;

    r._reconstructions["Family5"][p_tree->find_descendant("AB")].credible_interval = pair<double, double>(pv::to_computational_space(10), pv::to_computational_space(12));
    CHECK(credible_interval_intersects_parent(p_tree->find_descendant("A"), fam, &r));
    CHECK(!credible_interval_intersects_parent(p_tree->find_descendant("B"), fam, &r));
}

TEST_CASE_FIXTURE(Reconstruction, "credible_interval_intersects_parent is correct for internal nodes")
{
    mock_reconstruction r;

    r._reconstructions["Family5"][p_tree->find_descendant("AB")].credible_interval = pair<double, double>(pv::to_computational_space(10), pv::to_computational_space(12));
    r._reconstructions["Family5"][p_tree->find_descendant("CD")].credible_interval = pair<double, double>(pv::to_computational_space(3), pv::to_computational_space(7));
    r._reconstructions["Family5"][p_tree.get()].credible_interval = pair<double, double>(pv::to_computational_space(5), pv::to_computational_space(9));
    CHECK(!credible_interval_intersects_parent(p_tree->find_descendant("AB"), fam, &r));
    CHECK(credible_interval_intersects_parent(p_tree->find_descendant("CD"), fam, &r));
}

TEST_CASE_FIXTURE(Reconstruction, "credible_interval_intersects_parent returns false at root")
{
    mock_reconstruction r;

    CHECK(credible_interval_intersects_parent(p_tree.get(), fam, &r));
}
