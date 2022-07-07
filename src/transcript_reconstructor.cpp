#include <cmath>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <iterator>

#include "doctest.h"
#include "easylogging++.h"

#include "transcript_reconstructor.h"
#include "sigma.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "DiffMat.h"
#include "newick_ape_loader.h"
#include "proportional_variance.h"
#include "probability.h"

using namespace std;
using namespace Eigen;
namespace pv = proportional_variance;

transcript_reconstructor::transcript_reconstructor(const sigma_squared* p_sigma, const clade* p_tree, const matrix_cache* p_cache)
    : _p_sigma(p_sigma),
    _p_tree(p_tree),
    _p_cache(p_cache)
{
    for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
    {
        all_node_Ls[*it] = VectorXd::Zero(200);
    }
}

clademap<double> transcript_reconstructor::reconstruct_gene_transcript(const gene_transcript& t, int upper_bound)
{
    inference_pruner pruner(*_p_cache, _p_sigma, nullptr, _p_tree, 1.0);
    return pruner.reconstruct(t, upper_bound);
}

string newick_node(const clade *node, const cladevector& order, bool significant, std::function<std::string(const clade *c)> textwriter)
{
    ostringstream ost;
    ost << clade_index_or_name(node, order) << (significant ? "*" : "") << "_" << textwriter(node);

    if (!node->is_root())
        ost << ':' << node->get_branch_length();

    return ost.str();
}

void reconstruction::print_node_change(std::ostream& ost, const cladevector& order, transcript_vector& gene_families, const clade* p_tree)
{
    print_family_clade_table(ost, order, gene_families, p_tree, [this, &gene_families](int family_index, const clade* c) {
        ostringstream ost;
        ost << showpos << get_difference_from_parent(gene_families[family_index], c);
        return ost.str();
        });
}

void reconstruction::print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order, transcript_vector& gene_transcripts) {
    clademap<pair<int, int>> increase_decrease_map;

    for (auto& transcript : gene_transcripts)
    {
        for (auto node : order)
        {
            if (node)
            {
                double val = get_difference_from_parent(transcript, node);
                if (val > 0)
                    increase_decrease_map[node].first++;
                if (val < 0)
                    increase_decrease_map[node].second++;
            }
        }
    }

    ost << "#Taxon_ID\tIncrease\tDecrease\n";
    for (auto& it : increase_decrease_map) {
        ost << clade_index_or_name(it.first, order) << "\t";
        ost << it.second.first << "\t";
        ost << it.second.second << endl;
    }
}

void write_node_ordered(std::ostream& ost, std::string title, const cladevector& order, std::function<string(const clade* c)> f)
{
    ost << title;
    for (auto node : order)
    {
        if (node)
        {
            ost << "\t" << f(node);
        }
    }
    ost << endl;
}

void reconstruction::print_family_clade_table(std::ostream& ost, const cladevector& order, transcript_vector& gene_families, const clade* p_tree, std::function<string(int family_index, const clade *c)> get_family_clade_value)
{
    write_node_ordered(ost, "FamilyID", order, [order](const clade* c) { return clade_index_or_name(c, order); });
    for (size_t i = 0; i < gene_families.size(); ++i)
    {
        write_node_ordered(ost, gene_families[i].id(), order, [i, get_family_clade_value](const clade* c) { return get_family_clade_value(i, c); });
    }
}

void reconstruction::print_reconstructed_states(std::ostream& ost, const cladevector& order, transcript_vector& transcripts, const clade* p_tree, double test_pvalue)
{
    ost << "#nexus\nBEGIN TREES;\n";
    for (size_t i = 0; i < transcripts.size(); ++i)
    {
        auto& gene_transcript = transcripts[i];
        auto g = [gene_transcript, this](const clade* node) {
            return std::to_string(pv::to_user_space(get_node_value(gene_transcript, node)));
        };

        function<string(const clade*)> text_func;
        text_func = [g, order](const clade* node) {
            return newick_node(node, order, false, g);
        };

        ost << "  TREE " << gene_transcript.id() << " = ";
        p_tree->write_newick(ost, text_func);

        ost << ';' << endl;
    }
    ost << "\nEND;\n";
    write_nexus_extensions(ost);
    
}

void reconstruction::print_node_values(std::ostream& ost, const cladevector& order, transcript_vector& transcripts, const clade* p_tree)
{
    print_family_clade_table(ost, order, transcripts, p_tree, [this, transcripts](int family_index, const clade* c) {
        auto& gf = transcripts[family_index];
        if (c->is_leaf())
            return to_string(pv::to_user_space(gf.get_expression_value(c->get_taxon_name())));
        else
            return to_string(pv::to_user_space(get_node_value(gf, c)));
        });
}

void reconstruction::write_results(std::string model_identifier, 
    std::string output_prefix, 
    const clade *p_tree, 
    transcript_vector& families, 
    double test_pvalue)
{
    //cladevector order;
    //for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&order](const clade* c) { order.push_back(c); });
    auto order = get_ape_order(p_tree);

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing reconstructed states";
    std::ofstream ofst(filename(model_identifier + "_asr", output_prefix, "tre"));
    print_reconstructed_states(ofst, order, families, p_tree, test_pvalue);

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing node transcript values";
    std::ofstream counts(filename(model_identifier + "_count", output_prefix, "tab"));
    print_node_values(counts, order, families, p_tree);

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing node changes";
    std::ofstream change(filename(model_identifier + "_change", output_prefix, "tab"));
    print_node_change(change, order, families, p_tree);

    VLOG(TRANSCRIPT_RECONSTRUCTION) << "writing clade results";
    std::ofstream clade_results(filename(model_identifier + "_clade_results", output_prefix));
    print_increases_decreases_by_clade(clade_results, order, families);

    print_additional_data(order, families, output_prefix);
}

double reconstruction::get_difference_from_parent(const gene_transcript& gf, const clade* c)
{
    if (c->is_root())
        return 0;

    return pv::to_user_space(get_node_value(gf, c)) - pv::to_user_space(get_node_value(gf, c->get_parent()));
}

class Reconstruction
{
public:
    gene_transcript fam;
    unique_ptr<clade> p_tree;
    cladevector order;

    Reconstruction() : fam("Family5", "", "")
    {
        p_tree.reset(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

        fam.set_expression_value("A", 11);
        fam.set_expression_value("B", 2);
        fam.set_expression_value("C", 5);
        fam.set_expression_value("D", 6);

        vector<string> nodes{ "A", "B", "C", "D", "AB", "CD", "ABCD" };
        order.resize(nodes.size());
        const clade* t = p_tree.get();
        transform(nodes.begin(), nodes.end(), order.begin(), [t](string s) { return t->find_descendant(s); });

    }
};

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript assigns actual values to leaves")
{
    sigma_squared sig(0.1);
    fam.set_expression_value("A", 3.7);
    fam.set_expression_value("D", 9.4);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_values(), p_tree->get_branch_lengths(), 20);

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);

    auto actual = tr.reconstruct_gene_transcript(fam, 20);

    CHECK_EQ(3.7, actual[p_tree->find_descendant("A")]);
    CHECK_EQ(9.4, actual[p_tree->find_descendant("D")]);
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
    transcript_reconstructor tr(&sig, p_tree.get(), &calc);

    auto actual = tr.reconstruct_gene_transcript(fam, 200);

    CHECK_EQ(2.5, actual[p_tree->find_descendant("AB")]);

}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript returns parent above 0 with 0s at the leafs")
{
    sigma_squared sig(10.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_values(), p_tree->get_branch_lengths(), 20);

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);
    auto actual = tr.reconstruct_gene_transcript(fam, 20);
    CHECK_EQ(0.05, actual[p_tree->find_descendant("AB")]);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript returns correct value at root" * doctest::skip(true))
{
    // TODO: Create reasonable test case values
    sigma_squared sig(10.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_values(), p_tree->get_branch_lengths(), 20);

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);
    auto actual = tr.reconstruct_gene_transcript(fam, 20);
    CHECK_EQ(50, actual[p_tree.get()]);
}

class mock_reconstruction : public reconstruction
{
    double get_node_value(const gene_transcript& gf, const clade* c) const override
    { 
        if (_reconstructions.find(gf.id()) == _reconstructions.end())
            return 3;

        return _reconstructions.at(gf.id()).at(c);
    }
public:
    std::map<std::string, clademap<double>> _reconstructions;

};

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE_FIXTURE(Reconstruction, "print_node_counts")
{
    mock_reconstruction r;

    ostringstream ost;
    r.print_node_values(ost, order, std::vector<gene_transcript>({ fam }), p_tree.get());
    CHECK_STREAM_CONTAINS(ost, "FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>");
}

TEST_CASE_FIXTURE(Reconstruction, "print_node_counts handles no zero id")
{
    mock_reconstruction r;

    ostringstream ost;
    order.push_back(order[0]);
    order[0] = nullptr;
    r.print_node_values(ost, order, std::vector<gene_transcript>({ fam }), p_tree.get());
    CHECK_STREAM_CONTAINS(ost, "FamilyID\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>\tA<7>");
}

TEST_CASE("Reconstruction: print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    cladevector order{ 
        nullptr,
        p_tree->find_descendant("A"),
        p_tree->find_descendant("B"),
        p_tree->find_descendant("AB") };

    ostringstream empty;

    mock_reconstruction bmr;

    bmr.print_increases_decreases_by_clade(empty, order, {});
    CHECK_EQ(string("#Taxon_ID\tIncrease\tDecrease\n"), empty.str());

    bmr._reconstructions["myid"][p_tree->find_descendant("A")] = 22.11;
    bmr._reconstructions["myid"][p_tree->find_descendant("B")] = 9.3;
    bmr._reconstructions["myid"][p_tree->find_descendant("AB")] = 17.6;

    gene_transcript gf("myid", "", "");

    ostringstream ost;
    bmr.print_increases_decreases_by_clade(ost, order, { gf });
    CHECK_STREAM_CONTAINS(ost, "#Taxon_ID\tIncrease\tDecrease");
    CHECK_STREAM_CONTAINS(ost, "A<1>\t1\t0");
    CHECK_STREAM_CONTAINS(ost, "B<2>\t0\t1");
}

TEST_CASE_FIXTURE(Reconstruction, "get_difference_from_parent")
{
    mock_reconstruction r;
    r._reconstructions["Family5"][p_tree->find_descendant("A")] = pv::to_computational_space(22.11);
    r._reconstructions["Family5"][p_tree->find_descendant("AB")] = pv::to_computational_space(17.6);

    CHECK_EQ(doctest::Approx(4.51), r.get_difference_from_parent(fam, p_tree->find_descendant("A")));
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_index_in_angle_brackets_for_non_leaf")
{
    CHECK_EQ(string("<0>"), clade_index_or_name(p_tree.get(), { p_tree.get() }).c_str());
}

TEST_CASE_FIXTURE(Reconstruction, "clade_index_or_name__returns_node_name_plus_index_in_angle_brackets_for_leaf")
{
    auto a = p_tree->find_descendant("A");
    CHECK_EQ(string("A<1>"), clade_index_or_name(a, { p_tree.get(), a }).c_str());
}

#if 0
TEST_CASE_FIXTURE(Reconstruction, "print_branch_probabilities__shows_NA_for_invalids")
{
    gene_transcript gf("Family5", "", "");
    std::ostringstream ost;
    branch_probabilities probs;
    for (auto c : order)
        probs.set(gf, c, 0.05);
    probs.set(gf, p_tree->find_descendant("B"), branch_probabilities::invalid());
    probs.set(gf, p_tree.get(), branch_probabilities::invalid());

    print_branch_probabilities(ost, order, { fam }, probs);
    CHECK_STREAM_CONTAINS(ost, "FamilyID\tA<0>\tB<1>\tC<2>\tD<3>\t<4>\t<5>\t<6>");
    CHECK_STREAM_CONTAINS(ost, "Family5\t0.05\tN/A\t0.05\t0.05\t0.05\t0.05\tN/A\n");
}

TEST_CASE_FIXTURE(Reconstruction, "print_branch_probabilities__skips_families_without_reconstructions")
{
    std::ostringstream ost;
    branch_probabilities probs;

    print_branch_probabilities(ost, order, { fam }, probs);
    CHECK(ost.str().find("Family5") == string::npos);
}

#endif
