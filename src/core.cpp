#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

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
#include "prior.h"
#include "freerate_paired_model.h"
#include "freerate_global_model.h"

using namespace std;

std::vector<model *> build_models(const input_parameters& user_input, user_data& user_data) {

    model *p_model = NULL;

    std::vector<gene_transcript> *p_gene_transcripts = &user_data.gene_transcripts;

    if (user_input.is_simulating) {
        p_gene_transcripts = NULL;
    }

    if (user_input.fixed_alpha > 0 || user_input.n_gamma_cats > 1)
    {
        auto gmodel = new gamma_model(user_data.p_sigma, &user_data.gene_transcripts, 
            user_input.n_gamma_cats, user_input.fixed_alpha, user_data.p_error_model);
#ifndef SILENT
        if (user_input.fixed_alpha >= 0 && !user_input.is_simulating)
            gmodel->write_probabilities(cout);
#endif
        p_model = gmodel;
    }
    else if (!user_input.free_rate.empty())
    {
        if (user_input.free_rate == "global")
            p_model = new freerate_global_model(user_input.input_file_has_ratios);
        else
            p_model = new freerate_paired_model(user_input.input_file_has_ratios);
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

        p_model = new base_model(user_data.p_sigma, p_gene_transcripts, p_error_model);
    }

    return std::vector<model *>{p_model};
}

model::model(sigma_squared* p_sigma,
    const vector<gene_transcript> *p_gene_transcripts,
    error_model *p_error_model) :
    _ost(cout), _p_sigma(p_sigma),  _p_error_model(p_error_model) 
{
    if (p_gene_transcripts)
        references = build_reference_list(*p_gene_transcripts);
}

void model::write_vital_statistics(std::ostream& ost, const clade *p_tree, double final_likelihood, const input_parameters& p)
{
    ost << PROJECT_NAME " " PROJECT_VER << endl;
    ost << "Command line: " << p.command_line << endl;
    ost << "Final Likelihood (-lnL): " << final_likelihood << endl;

    auto p_sigma = get_sigma();
    if (p_sigma)
        ost << "Sigma2: " << *p_sigma << endl;
    
    ost << "IDs of Nodes: ";
    p_tree->write_newick(ost, [](std::ostream& ost, const clade*c) {
        ost << "<" << c->get_ape_index() << ">";
        });
    ost << endl;

    if (_p_error_model)
        ost << "Epsilon: " << _p_error_model->get_epsilons()[0] << endl;

    get_monitor().log(ost);

    write_extra_vital_statistics(ost);
}

sigma_squared* model::get_simulation_sigma()
{
    return new sigma_squared(*_p_sigma);
}

void model::write_error_model(int max_transcript_size, std::ostream& ost) const
{
    auto em = _p_error_model;
    if (!em)
    {
        em = new error_model();
        em->set_probabilities(max_transcript_size, { 0, 1, 0 });
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
        auto worst_performing_transcript = std::max_element(failure_count.begin(), failure_count.end(), failures);
        if (worst_performing_transcript->second * 5 > (attempts - rejects))   
        {
            ost << "The following transcripts had failure rates >20% of the time:\n";
            for (auto& a : this->failure_count)
            {
                if (a.second * 5 > (attempts - rejects))
                    ost << a.first << " had " << a.second << " failures\n";
            }
        }
    }
}

bool verify_results(const std::vector<size_t>& references, const std::vector<std::vector<double>>& partial_likelihoods, int& error)
{
    for (size_t i = 0; i < references.size(); ++i)
    {
        if (references[i] == i && partial_likelihoods[i].empty())
        {
            error = i;
            return false;
        }
    }
    return true;
}

double compute_distribution_mean(const clade* p_tree, const transcript_vector& transcripts)
{
    if (transcripts.empty()) throw runtime_error("No gene transcripts provided");
    if (!p_tree) throw runtime_error("No tree provided");

    vector<double> variances;
    for (auto& tt : transcripts)
    {
        auto species = tt.get_species();
        if (species.size() < 2) continue;

        vector<double> v(species.size());
        transform(species.begin(), species.end(), v.begin(), [&tt](string s) {return tt.get_expression_value(s);  });
        double sz = v.size();
        
        auto mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
        variances.push_back(std::accumulate(v.begin(), v.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
            return accumulator + ((val - mean) * (val - mean) / (sz - 1));
            }));
    }

    double species_variance = std::accumulate(variances.begin(), variances.end(), 0.0) / double(variances.size());
    return species_variance / p_tree->distance_from_root_to_tip();
}


#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

class mock_model : public model {
    // Inherited via model
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override { return nullptr; }
    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups) override { return nullptr; }
    bool _invalid_likelihood = false;
public:
    mock_model(sigma_squared*s) : model(s, NULL, NULL) {}
    virtual double infer_transcript_likelihoods(const user_data& ud, const sigma_squared* p_ss) override { return 0;  }

    virtual std::string get_name() const override { return "Mock"; }
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

TEST_CASE("Inference, event_monitor_shows_poor_performing_families")
{
    event_monitor evm;

    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Started();
    evm.Event_InferenceAttempt_Saturation("test");
    ostringstream ost;

    evm.log(ost);
    CHECK_STREAM_CONTAINS(ost, "2 values were attempted (0% rejected)");
    CHECK_STREAM_CONTAINS(ost, "The following transcripts had failure rates >20% of the time:");
    CHECK_STREAM_CONTAINS(ost, "test had 1 failures");
}

TEST_CASE("Inference: create_gamma_model_if__n_gamma_cats__provided")
{
    input_parameters params;
    params.input_file_path = "foo";
    params.n_gamma_cats = 3;
    params.fixed_alpha = 1.5;
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<gamma_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("build_models creates a freerate model")
{
    input_parameters params;
    params.free_rate = "paired";
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<freerate_paired_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("build_models creates a global freerate model")
{
    input_parameters params;
    params.free_rate = "global";
    user_data data;
    auto models = build_models(params, data);
    CHECK_EQ(1, models.size());
    CHECK(dynamic_cast<freerate_global_model*>(models[0]));
    for (auto m : models)
        delete m;
}

TEST_CASE("verify_results fails if empty partial likelihoods")
{
    std::vector<size_t> references = {0, 1, 2, 3};
    std::vector<std::vector<double>> partial_likelihoods = {{1.0, 2.0}, {}, {4.0, 5.0, 6.0}, {}};

    int error = -1;
    bool result = verify_results(references, partial_likelihoods, error);

    CHECK_FALSE(result);
    CHECK_EQ(error, 1);
}

TEST_CASE("verify_results succeeds if empty partial likelihoods are referenced")
{
    std::vector<size_t> references = {0, 0, 2, 3};
    std::vector<std::vector<double>> partial_likelihoods = {{1.0, 2.0}, {}, {3.0}, {4.0, 5.0, 6.0}};

    int error = -1;
    bool result = verify_results(references, partial_likelihoods, error);

    CHECK(result);
}

TEST_CASE("verify_results succeeds if no empty partial likelihoods")
{
    std::vector<size_t> references = {0, 1, 2, 3};
    std::vector<std::vector<double>> partial_likelihoods = {{1.0, 2.0}, {7.0}, {3.0}, {4.0, 5.0, 6.0}};

    int error = -1;
    bool result = verify_results(references, partial_likelihoods, error);

    CHECK(result);
}

TEST_CASE("compute_distribution_mean")
{
    vector<gene_transcript> transcripts;
    transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    transcripts[0].set_expression_value("A", 1);
    transcripts[0].set_expression_value("B", 2);

    unique_ptr<clade> tree(parse_newick("(A:1,B:1);"));

    CHECK_EQ(0.5, compute_distribution_mean(tree.get(), transcripts));
}

TEST_CASE("compute_distribution_mean skips transcripts with less than two values")
{
    vector<gene_transcript> transcripts;
    transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    transcripts.push_back(gene_transcript("TestFamily2", "", ""));
    transcripts[0].set_expression_value("A", 5);
    transcripts[1].set_expression_value("A", 1);
    transcripts[1].set_expression_value("B", 2);

    unique_ptr<clade> tree(parse_newick("(A:1,B:1);"));

    CHECK_EQ(0.5, compute_distribution_mean(tree.get(), transcripts));
}