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
#include "reconstruction.h"
#include "prior.h"

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
    ost << "Sigma2: " << *get_sigma() << endl;
    
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

double computational_space_prior(double val, const prior *p_prior)
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    return exp(val) * p_prior->pdf(exp(val));
#else
    return p_prior->pdf(val);
#endif

}

double compute_prior_likelihood(const vector<double>& partial_likelihood, const vector<double>& priors)
{
    std::vector<double> full(partial_likelihood.size());
    std::transform(partial_likelihood.begin(), partial_likelihood.end(), priors.begin(), full.begin(), std::multiplies<double>());
    std::transform(full.begin(), full.end(), full.begin(), [](double d) {
        return isnan(d) ? -numeric_limits<double>::infinity() : d;
        });

#ifdef USE_MAX_PROBABILITY
    double likelihood = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
#else
    double likelihood = accumulate(full.begin(), full.end(), 0.0, [](double a, double b) { return isinf(b) ? a : a+b; }); // sum over all sizes (Felsenstein's approach)
#endif
    return log(likelihood);
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



#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

class mock_model : public model {
    // Inherited via model
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override { return nullptr; }
    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups) override { return nullptr; }
    bool _invalid_likelihood = false;
public:
    mock_model(sigma_squared*s) : model(s, NULL, NULL) {}
    virtual double infer_transcript_likelihoods(const user_data& ud, const sigma_squared* p_ss) override { return 0;  }
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

