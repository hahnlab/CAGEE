#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

#include <boost/algorithm/string.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "core.h"
#include "user_data.h"
#include "matrix_cache.h"
#include "gamma_core.h"
#include "base_model.h"
#include "error_model.h"
#include "sigma.h"
#include "arguments.h"
#include "reconstruction.h"

using namespace std;

std::vector<model *> build_models(const input_parameters& user_input, user_data& user_data) {

    model *p_model = NULL;

    std::vector<gene_transcript> *p_gene_families = &user_data.gene_families;

    if (user_input.is_simulating) {
        p_gene_families = NULL;
    }

    if (user_input.fixed_alpha > 0 || user_input.n_gamma_cats > 1)
    {
        auto gmodel = new gamma_model(user_data.p_lambda, &user_data.gene_families, 
            user_input.n_gamma_cats, user_input.fixed_alpha, user_data.p_error_model);
#ifndef SILENT
        if (user_input.fixed_alpha >= 0 && !user_input.is_simulating)
            gmodel->write_probabilities(cout);
#endif
        p_model = gmodel;
    }
    else
    {
        error_model* p_error_model = user_data.p_error_model;
        if (user_input.use_error_model && !p_error_model)
        {
            p_error_model = new error_model();
            p_error_model->set_probabilities(0, { 0, .95, 0.05 });
            p_error_model->set_probabilities(200, { 0.05, .9, 0.05 });
        }

        p_model = new base_model(user_data.p_lambda, p_gene_families, p_error_model);
    }

    return std::vector<model *>{p_model};
}

int upper_bound_from_transcript_values(const vector<gene_transcript>& transcripts)
{
    vector<int> bounds(transcripts.size());
    transform(transcripts.begin(), transcripts.end(), bounds.begin(), [](const gene_transcript& gt) {

        return max(1.0, ceil(gt.get_max_expression_value() + MATRIX_SIZE_MULTIPLIER));
        });

    return *max_element(bounds.begin(), bounds.end());
}


model::model(sigma_squared* p_lambda,
    const vector<gene_transcript> *p_gene_families,
    error_model *p_error_model) :
    _ost(cout), _p_sigma(p_lambda),  _p_error_model(p_error_model) 
{
    if (p_gene_families)
        references = build_reference_list(*p_gene_families);
}

void model::write_vital_statistics(std::ostream& ost, const clade *p_tree, double final_likelihood, const input_parameters& p)
{
    ost << PROJECT_NAME " " PROJECT_VER << endl;
    ost << "Command line: " << p.command_line << endl;
    ost << "Final Likelihood (-lnL): " << final_likelihood << endl;
    ost << "Sigma2: " << *get_sigma() << endl;
    
    ost << "IDs of Nodes: ";
    p_tree->write_newick(ost, [](const clade*c) {
        return "<" + to_string(c->get_ape_index()) + ">";
        });
    ost << endl;

    if (_p_error_model)
        ost << "Epsilon: " << _p_error_model->get_epsilons()[0] << endl;

    get_monitor().log(ost);

    write_extra_vital_statistics(ost);
}

sigma_squared* model::get_simulation_lambda()
{
    return _p_sigma->clone();
}

void model::write_error_model(int max_family_size, std::ostream& ost) const
{
    auto em = _p_error_model;
    if (!em)
    {
        em = new error_model();
        em->set_probabilities(max_family_size, { 0, 1, 0 });
    }
    write_error_model_file(ost, *em);
}

void event_monitor::Event_InferenceAttempt_Started() 
{ 
    attempts++;
}

void event_monitor::log(el::base::type::ostream_t& ost) const
{
    if (attempts == 0)
    {
        ost << "No attempts made\n";
        return;
    }
    ost << this->attempts << " values were attempted (" << round(double(rejects) / double(attempts) * 100) << "% rejected)\n";
    if (!failure_count.empty())
    {
        auto failures = [](const pair<string, int>& a, const pair<string, int>& b) { return a.second < b.second; };
        auto worst_performing_family = std::max_element(failure_count.begin(), failure_count.end(), failures);
        if (worst_performing_family->second * 5 > (attempts - rejects))    // at least one family had 20% rejections
        {
            ost << "The following families had failure rates >20% of the time:\n";
            for (auto& a : this->failure_count)
            {
                if (a.second * 5 > (attempts - rejects))
                    ost << a.first << " had " << a.second << " failures\n";
            }
        }
    }
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

//! Create a sigma based on the sigma tree model the user passed.
/// Called when the user has provided no sigma value and one must
/// be estimated. If the p_sigma_tree is NULL, uses a single
/// sigma; otherwise uses the number of unique sigmas in the provided
/// tree
sigma_squared* initialize_search_sigma(clade* p_sigma_tree, const std::vector<string>& sample_groups)
{
    sigma_squared* p_sigma = NULL;
    if (p_sigma_tree != NULL)
    {
        std::set<int> unique_sigmas;
        auto fn = [&unique_sigmas](const clade* p_node) { unique_sigmas.insert(p_node->get_sigma_index()); };
        p_sigma_tree->apply_prefix_order(fn);
        auto node_name_to_sigma_index = p_sigma_tree->get_sigma_index_map();
        p_sigma = new sigma_squared(node_name_to_sigma_index, std::vector<double>(unique_sigmas.size()), sigma_type::lineage_specific);
        LOG(INFO) << "Searching for " << unique_sigmas.size() << " sigmas" << endl;
    }
    else if (!sample_groups.empty())
    {
        auto sample_name_to_sigma_index = get_sigma_index_map(sample_groups);
        p_sigma = new sigma_squared(sample_name_to_sigma_index, std::vector<double>(sample_groups.size()), sigma_type::sample_specific);
        LOG(INFO) << "Searching for " << sample_groups.size() << " sigmas" << endl;
    }
    else
    {
        p_sigma = new sigma_squared(0.0);
    }

    return p_sigma;
}

TEST_CASE("initialize_search_sigma returns single sigma if no arguments")
{
    unique_ptr<sigma_squared> sig(initialize_search_sigma(nullptr, vector<string>()));
    CHECK_EQ(1, sig->get_values().size());
}

TEST_CASE("initialize_search_sigma returns lineage-specific with a sigma tree")
{
    string s = "((((cat:1,horse:1):1,cow:1):1,(((((chimp:2,human:2):2,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:1,mouse:1):1);";
    unique_ptr<clade> p_tree(parse_newick(s, true));

    unique_ptr<sigma_squared> sig(initialize_search_sigma(p_tree.get(), vector<string>()));
    REQUIRE_EQ(2, sig->get_values().size());

    double values[2] = { 7,14 };
    sig->update(values);
    CHECK_EQ(7, sig->get_named_value(p_tree->find_descendant("cow"), gene_transcript()));
    CHECK_EQ(14, sig->get_named_value(p_tree->find_descendant("chimp"), gene_transcript()));
}

TEST_CASE("initialize_search_sigma returns sample-specific with sample groups")
{
    unique_ptr<sigma_squared> sig(initialize_search_sigma(nullptr, vector<string>({"heart", "lungs"})));
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

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

class mock_model : public model {
    // Inherited via model
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override { return nullptr; }
    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior) override { return nullptr; }
    bool _invalid_likelihood = false;
public:
    mock_model(sigma_squared*s) : model(s, NULL, NULL) {}
    virtual double infer_family_likelihoods(const user_data& ud, const sigma_squared* p_lambda, const std::gamma_distribution<double>& prior) override { return 0;  }
};

TEST_CASE("Model: write_vital_statistics")
{
    sigma_squared s(75.5);
    mock_model model(&s);
    std::ostringstream ost;
    input_parameters p;
    p.command_line = "cagee -t tree.txt -i transcript.txt";
    model.write_vital_statistics(ost, new clade("A", 5), 0.01, p);
    CHECK_STREAM_CONTAINS(ost, PROJECT_NAME " " PROJECT_VER);
    CHECK_STREAM_CONTAINS(ost, "Final Likelihood (-lnL): 0.01");
    CHECK_STREAM_CONTAINS(ost, "Sigma2: 75.5");
    CHECK_STREAM_CONTAINS(ost, "No attempts made");
    CHECK_STREAM_CONTAINS(ost, "Command line: cagee -t tree.txt -i transcript.txt");
}

TEST_CASE("write_vital_statistics writes node IDs")
{
    string s = "((((cat:1,horse:1):1,cow:1):1,(((((chimp:2,human:2):2,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:1,mouse:1):1);";
    unique_ptr<clade> p_tree(parse_newick(s, false));

    sigma_squared ss(75.5);
    mock_model model(&ss);
    std::ostringstream ost;
    input_parameters p;
    model.write_vital_statistics(ost, p_tree.get(), 0.01, p);
    CHECK_STREAM_CONTAINS(ost, "IDs of Nodes: ((((<1>,<2>)<16>,<3>)<15>,(((((<4>,<5>)<21>,<6>)<20>,<7>)<19>,(<8>,<9>)<22>)<18>,<10>)<17>)<14>,(<11>,<12>)<23>)<13>\n");
}

/// Bounds are next largest integer if in log space.
TEST_CASE("upper_bound_from_transcript_values returns next multiple of 5")
{
    gene_transcript gt;
    gt.set_expression_value("A", .3);
    gt.set_expression_value("B", .8);
    CHECK_EQ(4, upper_bound_from_transcript_values({ gt }));

    gt.set_expression_value("A", 1.1);
    CHECK_EQ(5, upper_bound_from_transcript_values({ gt }));

    gt.set_expression_value("B", 7.3);
    CHECK_EQ(11, upper_bound_from_transcript_values({ gt }));
}

TEST_CASE("upper_bound_from_transcript_values never returns less than MATRIX_SIZE_MULTIPLIER")
{
    gene_transcript gt;
    gt.set_expression_value("A", 0.00005);
    gt.set_expression_value("B", 0.00004);
    CHECK_EQ(MATRIX_SIZE_MULTIPLIER+1, upper_bound_from_transcript_values({ gt }));
}

TEST_CASE("upper_bound_from_transcript_values never returns less than MATRIX_SIZE_MULTIPLIER even if all values are very small")
{
    gene_transcript gt;
    gt.set_expression_value("A", 0.0000000002);
    gt.set_expression_value("B", 0.0000000005);
    CHECK_EQ(MATRIX_SIZE_MULTIPLIER+1, upper_bound_from_transcript_values({ gt }));
}

TEST_CASE("upper_bound_from_transcript_values returns MATRIX_SIZE_MULTIPLIER if all values are zero")
{
    gene_transcript gt;
    gt.set_expression_value("A", 0);
    gt.set_expression_value("B", 0);
    CHECK_EQ(MATRIX_SIZE_MULTIPLIER, upper_bound_from_transcript_values({ gt }));
}

