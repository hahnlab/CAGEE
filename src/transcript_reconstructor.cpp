#include <cmath>
#include <sstream>
#include <algorithm>
#include <fstream>

#include "doctest.h"
#include "easylogging++.h"

#include "transcript_reconstructor.h"
#include "sigma.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "DiffMat.h"

using namespace std;
using namespace Eigen;

transcript_reconstructor::transcript_reconstructor(const sigma* p_sigma, const clade* p_tree, const matrix_cache* p_cache)
    : _p_sigma(p_sigma),
    _p_tree(p_tree),
    _p_cache(p_cache)
{
    for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
    {
        all_node_Ls[*it] = VectorXd::Zero(200);
    }
}

double get_value(const gene_transcript& t, const VectorXd& likelihood)
{
    auto bound = get_upper_bound(t);
    Index i;
    auto m = likelihood.maxCoeff(&i);
    return (float(i) / 200.0) * bound;
}

clademap<double> transcript_reconstructor::reconstruct_gene_transcript(const gene_transcript& t)
{
    IOFormat CleanFmt(4, DontAlignCols, " ", " ", "", "");
    for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
    {
        const clade* c = *it;
        if (c->is_leaf())
        {
            std::pair<double, double> bounds(0, get_upper_bound(t));

            all_node_Ls[c] = VectorPos_bounds(t.get_expression_value(c->get_taxon_name()), DISCRETIZATION_RANGE, bounds);
            VLOG(TRANSCRIPT_RECONSTRUCTION) << c->get_taxon_name() << " " << all_node_Ls[c].format(CleanFmt) << endl;
        }
        else
        {
            // Each child should be propagated separately, and then take the mean of the resulting vectors
            //    So mean(child1 * m, child2 * m)
            
            vector<VectorXd> child_likelihoods(2);
            transform(c->descendant_begin(), c->descendant_end(), child_likelihoods.begin(), [this, t, &CleanFmt](const clade* child) {
                auto matrix = _p_cache->get_matrix(child->get_branch_length(), _p_sigma->get_named_value(child, t), get_upper_bound(t));
                VectorXd result = matrix * all_node_Ls[child];
                VLOG(TRANSCRIPT_RECONSTRUCTION) << "ConvPropBounds * L[" << child->get_taxon_name() << "] " << result.format(CleanFmt) << endl;

                return result;
                });
            for (int i = 0; i < DISCRETIZATION_RANGE; ++i)
            {
                all_node_Ls[c][i] = (child_likelihoods[0][i] + child_likelihoods[1][i]) / 2.0;
            }
            VLOG(TRANSCRIPT_RECONSTRUCTION) << c->get_taxon_name() << " " << all_node_Ls[c].format(CleanFmt) << endl;
        }

    }

    clademap<double> result;
    for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
    {
        result[*it] = (*it)->is_leaf() ? t.get_expression_value((*it)->get_taxon_name()) : get_value(t, all_node_Ls[*it]);
    }
    return result;
}

string newick_node(const clade *node, const cladevector& order, bool significant, std::function<std::string(const clade *c)> textwriter)
{
    ostringstream ost;
    ost << clade_index_or_name(node, order) << (significant ? "*" : "") << "_" << textwriter(node);

    if (!node->is_root())
        ost << ':' << node->get_branch_length();

    return ost.str();
}

void reconstruction::print_node_change(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree)
{
    print_family_clade_table(ost, order, gene_families, p_tree, [this, &gene_families](int family_index, const clade* c) {
        ostringstream ost;
        ost << showpos << get_difference_from_parent(gene_families[family_index], c);
        return ost.str();
        });
}


void reconstruction::print_increases_decreases_by_clade(std::ostream& ost, const cladevector& order, familyvector& gene_families) {
    clademap<pair<int, int>> increase_decrease_map;

    for (size_t j = 0; j < gene_families.size(); ++j) {
        for (size_t i = 0; i < order.size(); ++i)
        {
            int val = get_difference_from_parent(gene_families[j], order[i]);
            if (val > 0)
                increase_decrease_map[order[i]].first++;
            if (val < 0)
                increase_decrease_map[order[i]].second++;
        }
    }

    ost << "#Taxon_ID\tIncrease\tDecrease\n";
    for (auto& it : increase_decrease_map) {
        ost << clade_index_or_name(it.first, order) << "\t";
        ost << it.second.first << "\t";
        ost << it.second.second << endl;
    }
}

void reconstruction::print_family_clade_table(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree, std::function<string(int family_index, const clade *c)> get_family_clade_value)
{
    ost << "FamilyID";
    for (auto c : order)
    {
        ost << "\t" << clade_index_or_name(c, order);
    }
    ost << endl;
    for (size_t i = 0; i < gene_families.size(); ++i)
    {
        ost << gene_families[i].id();
        for (auto node : order)
        {
            ost << "\t";
            ost << get_family_clade_value(i, node);
        }
        ost << endl;
    }
}

void print_branch_probabilities(std::ostream& ost, const cladevector& order, const vector<gene_transcript>& gene_families, const branch_probabilities& branch_probabilities)
{
    ost << "#FamilyID\t";
    for (auto& it : order) {
        ost << clade_index_or_name(it, order) << "\t";
    }
    ost << endl;

    for (auto& gf : gene_families) 
    {
        if (branch_probabilities.contains(gf))
        {
            ost << gf.id();
            for (auto c : order)
            {
                ost << '\t';
                if (branch_probabilities.at(gf, c)._is_valid)
                    ost << branch_probabilities.at(gf, c)._value;
                else
                    ost << "N/A";
            }
            ost << endl;
        }
    }

}

void reconstruction::print_reconstructed_states(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree, double test_pvalue, const branch_probabilities& branch_probabilities)
{
    ost << "#nexus\nBEGIN TREES;\n";
    for (size_t i = 0; i < gene_families.size(); ++i)
    {
        auto& gene_transcript = gene_families[i];
        auto g = [gene_transcript, this](const clade* node) {
            return std::to_string(get_node_count(gene_transcript, node));
        };

        function<string(const clade*)> text_func;
        if (branch_probabilities.contains(gene_transcript))
        {
            auto is_significant = [&branch_probabilities, test_pvalue, gene_transcript](const clade* node) {
                const auto& p = branch_probabilities.at(gene_transcript, node);
                return p._is_valid ? p._value < test_pvalue : false;
            };

            text_func = [g, order, is_significant, &gene_transcript](const clade* node) {
                return newick_node(node, order, is_significant(node), g);
            };
        }
        else
        {
            text_func = [g, order](const clade* node) {
                return newick_node(node, order, false, g);
            };
        }

        ost << "  TREE " << gene_transcript.id() << " = ";
        p_tree->write_newick(ost, text_func);

        ost << ';' << endl;
    }
    ost << "\nEND;\n";
    write_nexus_extensions(ost);
    
}

void reconstruction::print_node_counts(std::ostream& ost, const cladevector& order, familyvector& gene_families, const clade* p_tree)
{
    print_family_clade_table(ost, order, gene_families, p_tree, [this, gene_families](int family_index, const clade* c) {
        auto& gf = gene_families[family_index];
        if (c->is_leaf())
            return to_string(gf.get_expression_value(c->get_taxon_name()));
        else
            return to_string(get_node_count(gf, c));
        });
}


void reconstruction::write_results(std::string model_identifier, 
    std::string output_prefix, 
    const clade *p_tree, 
    familyvector& families, 
    double test_pvalue, 
    const branch_probabilities& branch_probabilities)
{
    cladevector order;
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&order](const clade* c) { order.push_back(c); });

    std::ofstream ofst(filename(model_identifier + "_asr", output_prefix, "tre"));
    print_reconstructed_states(ofst, order, families, p_tree, test_pvalue, branch_probabilities);

    std::ofstream counts(filename(model_identifier + "_count", output_prefix, "tab"));
    print_node_counts(counts, order, families, p_tree);

    std::ofstream change(filename(model_identifier + "_change", output_prefix, "tab"));
    print_node_change(change, order, families, p_tree);

    std::ofstream clade_results(filename(model_identifier + "_clade_results", output_prefix));
    print_increases_decreases_by_clade(clade_results, order, families);

    std::ofstream branch_probabilities_file(filename(model_identifier + "_branch_probabilities", output_prefix, "tab"));
    print_branch_probabilities(branch_probabilities_file, order, families, branch_probabilities);

    print_additional_data(order, families, output_prefix);
}

int reconstruction::get_difference_from_parent(const gene_transcript& gf, const clade* c)
{
    if (c->is_root())
        return 0;

    return get_node_count(gf, c) - get_node_count(gf, c->get_parent()); 
}

branch_probabilities::branch_probability compute_viterbi_sum(const clade* c, 
    const gene_transcript& transcript, 
    const reconstruction* rec, 
    const matrix_cache& cache, 
    const sigma* p_lambda)
{
    if (c->is_root())
    {
        return branch_probabilities::invalid();
    }

    auto probs = cache.get_matrix(c->get_branch_length(), p_lambda->get_named_value(c, transcript), get_upper_bound(transcript));

    int parent_size = rec->get_node_count(transcript, c->get_parent());
    int child_size = rec->get_node_count(transcript, c);
    double result = 0;
    double calculated_probability = probs(parent_size, child_size);
    for (int m = 0; m < 200; m++)
    {
        double probability_to_m = probs(parent_size, m);
        if (probability_to_m == calculated_probability)
        {
            result += probability_to_m / 2.0;
        }
        else if (probability_to_m < calculated_probability)
        {
            result += probability_to_m;
        }
    }
    if (result < 0.05)
    {
        VLOG(1) << transcript.id() << ":" << c->get_taxon_name() << " probability " << result << " calculated for parent: " << parent_size << ", child: " << child_size;
    }
    else
    {
        VLOG(2) << transcript.id() << ":" << c->get_taxon_name() << " probability " << result << " calculated for parent: " << parent_size << ", child: " << child_size;
    }
    return branch_probabilities::branch_probability(result);
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
    sigma sig(0.1);
    fam.set_expression_value("A", 3.7);
    fam.set_expression_value("D", 9.4);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_lambdas(), set<int>({ get_upper_bound(fam) }), p_tree->get_branch_lengths());

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);

    auto actual = tr.reconstruct_gene_transcript(fam);

    CHECK_EQ(3.7, actual[p_tree->find_descendant("A")]);
    CHECK_EQ(9.4, actual[p_tree->find_descendant("D")]);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript calculates parent node values correctly")
{
    sigma sig(10.1);
    fam.set_expression_value("A", 45.2);
    fam.set_expression_value("B", 61.8);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_lambdas(), set<int>({ get_upper_bound(fam) }), p_tree->get_branch_lengths());

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);

    auto actual = tr.reconstruct_gene_transcript(fam);

    CHECK_EQ(47.0, actual[p_tree->find_descendant("AB")]);

}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript returns parent 0 with 0s at the leafs")
{
    sigma sig(10.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_lambdas(), set<int>({ get_upper_bound(fam) }), p_tree->get_branch_lengths());

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);
    auto actual = tr.reconstruct_gene_transcript(fam);
    CHECK_EQ(0, actual[p_tree->find_descendant("AB")]);
}

TEST_CASE_FIXTURE(Reconstruction, "reconstruct_gene_transcript returns correct value at root" * doctest::skip(true))
{
    // TODO: Create reasonable test case values
    sigma sig(10.1);
    fam.set_expression_value("A", 0);
    fam.set_expression_value("B", 0);

    matrix_cache calc;
    calc.precalculate_matrices(sig.get_lambdas(), set<int>({ get_upper_bound(fam) }), p_tree->get_branch_lengths());

    transcript_reconstructor tr(&sig, p_tree.get(), &calc);
    auto actual = tr.reconstruct_gene_transcript(fam);
    CHECK_EQ(50, actual[p_tree.get()]);
}

