#include <assert.h>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <random>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>

#include "doctest.h"
#include "easylogging++.h"

#include "gamma_core.h"
#include "gamma.h"
#include "reconstruction.h"
#include "matrix_cache.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "optimizer_scorer.h"
#include "simulator.h"
#include "sigma.h"
#include "DiffMat.h"
#include "proportional_variance.h"
#include "inference_pruner.h"
#include "prior.h"

using namespace std;
namespace pv = proportional_variance;

extern mt19937 randomizer_engine;

//! @brief Holds data for reconstructing a tree based on the Gamma model
//! \ingroup gamma
class gamma_model_reconstruction : public reconstruction
{
    discretized_gamma _gamma;
    virtual void write_nexus_extensions(std::ostream& ost) override;

public:
    gamma_model_reconstruction(transcript_vector& transcripts, discretized_gamma gamma) :
        reconstruction(transcripts),
        _gamma(gamma)
    {

    }

    gamma_model_reconstruction(transcript_vector& transcripts, replicate_model* p_model, discretized_gamma gamma) :
        reconstruction(transcripts, p_model),
        _gamma(gamma)
    {
    }

    void print_additional_data(std::string output_prefix) override;

    void print_category_likelihoods(std::ostream& ost);

    node_reconstruction get_internal_node_value(const gene_transcript& transcript, const clade* c) const;

    struct gamma_reconstruction {
        std::vector<clademap<node_reconstruction>> category_reconstruction;
        clademap<double> reconstruction;
        std::vector<double> _category_likelihoods;
    };

    std::map<std::string, gamma_reconstruction> _reconstructions;
};

gamma_model::gamma_model(sigma_squared* p_sigma, std::vector<gene_transcript>* p_gene_transcripts, int n_gamma_cats, double fixed_alpha, error_model* p_error_model) :
    model(p_sigma, p_gene_transcripts, p_error_model) {

    _gamma_cat_probs.resize(n_gamma_cats);
    _sigma_multipliers.resize(n_gamma_cats);
    if (p_gene_transcripts)
        _category_likelihoods.resize(p_gene_transcripts->size());
    set_alpha(fixed_alpha);
}

void gamma_model::write_extra_vital_statistics(std::ostream& ost)
{
    ost << "Alpha: " << _alpha << endl;
}

string comma_separated(const std::vector<double>& items)
{
    string s;
    for (auto i : items)
        s += (s.empty() ? "" : ",") + to_string(i);
    return s;
}

//! Set alpha for gamma distribution
void gamma_model::set_alpha(double alpha) {

    _alpha = alpha;
    if (_gamma_cat_probs.size() > 1)
        get_gamma(_gamma_cat_probs, _sigma_multipliers, alpha); // passing vectors by reference
}

void gamma_model::write_probabilities(ostream& ost)
{
    ost << "Alpha: " << _alpha << endl;
    ost << "Gamma cat probs are: " << comma_separated(_gamma_cat_probs) << endl;
    ost << "Sigma multipliers are: " << comma_separated(_sigma_multipliers) << endl;
}

sigma_squared* gamma_model::get_simulation_sigma()
{
    discrete_distribution<int> dist(_gamma_cat_probs.begin(), _gamma_cat_probs.end());
    return new sigma_squared(*_p_sigma, _sigma_multipliers[dist(randomizer_engine)]);
}

std::vector<double> gamma_model::get_posterior_probabilities(std::vector<double> cat_likelihoods)
{
    size_t process_count = cat_likelihoods.size();

    vector<double> numerators(process_count);
    transform(cat_likelihoods.begin(), cat_likelihoods.end(), _gamma_cat_probs.begin(), numerators.begin(), multiplies<double>());

    double denominator = accumulate(numerators.begin(), numerators.end(), 0.0);
    vector<double> posterior_probabilities(process_count);
    transform(numerators.begin(), numerators.end(), posterior_probabilities.begin(), [denominator](double d) { return d/denominator; });

    return posterior_probabilities;
}

bool gamma_model::can_infer() const
{
    if (!_p_sigma->is_valid())
        return false;

    if (_alpha < 0)
        return false;

    return true;
}

bool gamma_model::prune(const gene_transcript& transcript, const prior *p_prior, const matrix_cache& diff_mat, const sigma_squared*p_sigma,
    const clade *p_tree, std::vector<double>& category_likelihoods, boundaries bounds) 
{
    #if 0
    category_likelihoods.clear();

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        inference_pruner pruner(diff_mat, p_sigma, _p_error_model, nullptr, p_tree, bounds);
        auto partial_likelihood = pruner.prune(transcript);
        if (accumulate(partial_likelihood.begin(), partial_likelihood.end(), 0.0) == 0.0)
            return false;   // saturation

#ifdef MODEL_GENE_EXPRESSION_LOGS
        throw std::runtime_error("Log values not implemented for gamma yet");
#endif
        std::vector<double> full(partial_likelihood.size());
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = p_prior->pdf(j);
            full[j] = partial_likelihood[j] * eq_freq;
        }

#ifdef USE_MAX_PROBABILITY
        category_likelihoods.push_back(*max_element(full.begin(), full.end()) * _gamma_cat_probs[k]); // get max (CAFE's approach)
#else
        category_likelihoods.push_back(accumulate(full.begin(), full.end(), 0.0) * _gamma_cat_probs[k]); // sum over all sizes (Felsenstein's approach)
#endif
    }
#endif
    return true;
}

void flatten_vector(const vector<vector<double>>& v, vector<double>& result)
{
    for (auto& vv : v)
    {
        result.insert(result.end(), vv.begin(), vv.end());
    }
}

//! Infer bundle
double gamma_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma) {

    _monitor.Event_InferenceAttempt_Started();

    if (!can_infer())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    using namespace std;

    vector<double> all_bundles_likelihood(ud.gene_transcripts.size());

    vector<bool> failure(ud.gene_transcripts.size());

    vector<sigma_squared> sigmas;
    for (auto m : _sigma_multipliers)
    {
        sigmas.emplace_back(*p_sigma, m);
    }

    double final_likelihood = 0;
    for (auto s: sigmas)
    {
    
        matrix_cache calc;
            auto v = calc.create_vector();
    vector<double> priors(v.size());
    copy(v.begin(), v.end(), priors.begin());
    try
    {
        for (size_t j = 0; j < priors.size(); ++j) {
            double x = (double(j) + 0.5) * double(ud.bounds.second) / (priors.size() - 1);

            priors[j] = computational_space_prior(x, ud.p_prior);
        }
    }
    catch (std::domain_error& e)
    {
        LOG(DEBUG) << e.what();
        LOG(WARNING) << "Prior not valid for this sigma and data set";
        return -log(0);
    }
        std::vector<double> all_transcripts_likelihood(ud.gene_transcripts.size());

        calc.precalculate_matrices(s.get_values(),  ud.p_tree->get_branch_lengths(), ud.bounds);

        vector<vector<double>> partial_likelihoods(ud.gene_transcripts.size());

        vector<inference_pruner> pruners;
        pruners.reserve(ud.gene_transcripts.size());
        std::generate_n(std::back_inserter(pruners), ud.gene_transcripts.size(), [&]() {return inference_pruner(calc, &s, _p_error_model, ud.p_replicate_model, ud.p_tree, ud.bounds); });
    #pragma omp parallel for
        for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i) {
            if ((int)references[i] == i)
                partial_likelihoods[i] = pruners[i].prune(ud.gene_transcripts.at(i));
        }

    #pragma omp parallel for
        for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i) {

            all_transcripts_likelihood[i] = compute_prior_likelihood(partial_likelihoods[references[i]], priors);
            
        }
        double val = -std::accumulate(all_transcripts_likelihood.begin(), all_transcripts_likelihood.end(), 0.0);        
        LOG(DEBUG) << "  Sub-sigma " << s << " -lnL: " << val;
        final_likelihood += val;
    }

    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;

    return final_likelihood;

}

sigma_optimizer_scorer* gamma_model::get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups)
{
    bool estimate_sigma = data.p_sigma == NULL;
    bool estimate_alpha = _alpha <= 0.0;

    if (estimate_sigma && estimate_alpha)
    {
        _p_sigma = sigma_squared::create(data.p_sigma_tree, sample_groups);
        return new sigma_optimizer_scorer(this, data, _p_sigma);
    }
    else if (estimate_sigma && !estimate_alpha)
    {
        _p_sigma = sigma_squared::create(data.p_sigma_tree, sample_groups);
        return new sigma_optimizer_scorer(dynamic_cast<model *>(this), data, _p_sigma);
    }
    else if (!estimate_sigma && estimate_alpha)
    {
        _p_sigma = new sigma_squared(*data.p_sigma);
        return new sigma_optimizer_scorer(this, data);
    }
    else
    {
        return nullptr;
    }
}

clademap<double> get_weighted_averages(const std::vector<clademap<node_reconstruction>>& m, const vector<double>& probabilities)
{
    cladevector nodes(m[0].size());
    std::transform(m[0].begin(), m[0].end(), nodes.begin(), [](std::pair<const clade *, node_reconstruction> v) { return v.first;  });

    clademap<double> result;
    for (auto node : nodes)
    {
        double val = 0.0;
        for (size_t i = 0; i<probabilities.size(); ++i)
        {
            val += probabilities[i] * double(m[i].at(node).most_likely_value);
        }
        result[node] = val;
    }

    return result;
}

reconstruction* gamma_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *calc)
{
    LOG(INFO) << "Starting reconstruction processes for Gamma model";

    auto values = _p_sigma->get_values();
    vector<double> all;
    for (double multiplier : _sigma_multipliers)
    {
        for (double sigma : values)
        {
            all.push_back(sigma*multiplier);
        }
    }

    calc->precalculate_matrices(_p_sigma->get_values(), ud.p_tree->get_branch_lengths(), ud.bounds);

    gamma_model_reconstruction* result = new gamma_model_reconstruction(ud.gene_transcripts, ud.p_replicate_model, discretized_gamma(_alpha, _sigma_multipliers.size()));
    vector<gamma_model_reconstruction::gamma_reconstruction *> recs(ud.gene_transcripts.size());
    for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
    {
        recs[i] = &result->_reconstructions[ud.gene_transcripts[i].id()];
        result->_reconstructions[ud.gene_transcripts[i].id()]._category_likelihoods = _category_likelihoods[i];
        result->_reconstructions[ud.gene_transcripts[i].id()].category_reconstruction.resize(_sigma_multipliers.size());
    }

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        VLOG(1) << "Reconstructing for multiplier " << _sigma_multipliers[k];
        sigma_squared ml(*_p_sigma, _sigma_multipliers[k]);

        inference_pruner tr(&ml, ud.p_tree, ud.p_replicate_model, calc, ud.bounds);

        for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
        {
            recs[i]->category_reconstruction[k] = tr.reconstruct(ud.gene_transcripts[i]);
        }
    }

    for (auto reconstruction : recs)
    {
        // multiply every reconstruction by gamma_cat_prob
        reconstruction->reconstruction = get_weighted_averages(reconstruction->category_reconstruction, _gamma_cat_probs);
    }

    LOG(INFO) << "Done!\n";

    return result;
}

void gamma_model_reconstruction::write_nexus_extensions(std::ostream& ost)
{
    ost << "\nBEGIN SIGMA_MULTIPLIERS;\n";
    _gamma.write_multipliers(ost, false);
    ost << "END;\n\n";
}

node_reconstruction gamma_model_reconstruction::get_internal_node_value(const gene_transcript& transcript, const clade* c) const
{
    node_reconstruction nr;
    nr.most_likely_value = _reconstructions.at(transcript.id()).reconstruction.at(c);
    nr.credible_interval.first = nr.most_likely_value;
    nr.credible_interval.second = nr.most_likely_value;
    return nr;

}

void gamma_model_reconstruction::print_category_likelihoods(std::ostream& ost)
{
    ost << "Transcript ID\t";
    _gamma.write_multipliers(ost, true);
    ost << endl;

    for (auto& gf : _transcripts)
    {
        ost << gf.id() << '\t';
        auto rc = _reconstructions[gf.id()];
        ostream_iterator<double> ct(ost, "\t");
        copy(rc._category_likelihoods.begin(), rc._category_likelihoods.end(), ct);
        ost << endl;
    }
}

void gamma_model_reconstruction::print_additional_data(std::string output_prefix)
{
    std::ofstream cat_likelihoods(filename("category_likelihoods", output_prefix));
    print_category_likelihoods(cat_likelihoods);

}

discretized_gamma::discretized_gamma(double alpha, int bins) : _alpha(alpha) 
{
    if (bins > 1)
    {
        _sigma_multipliers.resize(bins);
        _gamma_cat_probs.resize(bins);
        get_gamma(_gamma_cat_probs, _sigma_multipliers, alpha); // passing vectors by reference
    }
}

sigma_squared* discretized_gamma::get_random_sigma(const sigma_squared& ss)
{
    discrete_distribution<int> dist(_gamma_cat_probs.begin(), _gamma_cat_probs.end());
    return new sigma_squared(ss, _sigma_multipliers[dist(randomizer_engine)]);       
}

vector<sigma_squared> discretized_gamma::get_discrete_sigmas(const sigma_squared &ss)
{
    vector<sigma_squared> sigmas;
    for (auto m : _sigma_multipliers)
    {
        sigmas.emplace_back(ss, m);
    }

    return sigmas;
}

void discretized_gamma::write_multipliers(std::ostream& ost, bool single_line) const
{
    for (auto& lm : _sigma_multipliers)
    {
        ost << (single_line ? "" : "  ") << lm << (single_line ? "\t" : ";\n");
    }
}

TEST_CASE("Inference: gamma_model__creates sigma optimizer__if_alpha_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, NULL, 4, 0.25, NULL);

    user_data data;
    data.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    data.gene_transcripts[0].set_expression_value("A", 1);
    data.gene_transcripts[0].set_expression_value("B", 2);
    data.p_tree = p_tree.get();

    auto opt = model.get_sigma_optimizer(data, vector<string>());

    REQUIRE(opt);
    CHECK_EQ("Optimizing Sigma ", opt->description());
    delete model.get_sigma();
}

TEST_CASE("Inference: gamma_model__creates__gamma_optimizer__if_sigma_provided")
{
    gamma_model model(NULL, NULL, 4, -1, NULL);

    user_data data;

    sigma_squared sl(0.05);
    data.p_sigma = &sl;

    auto opt = model.get_sigma_optimizer(data, vector<string>());

    REQUIRE(opt);
    CHECK_EQ("Optimizing Alpha ", opt->description());

    delete model.get_sigma();
}

TEST_CASE("Inference: gamma_model_creates_nothing_if_sigma_and_alpha_provided")
{
    gamma_model model(NULL, NULL, 4, .25, NULL);

    user_data data;

    sigma_squared sl(0.05);
    data.p_sigma = &sl;

    CHECK(model.get_sigma_optimizer(data, vector<string>()) == nullptr);
}


TEST_CASE("get_weighted_averages")
{
    clade c1;
    clade c2;

    clademap<node_reconstruction> rc1;
    node_reconstruction nr;
    nr.most_likely_value = 10;
    rc1[&c1] = nr;

    nr.most_likely_value = 2;
    rc1[&c2] = nr;

    clademap<node_reconstruction> rc2;
    nr.most_likely_value = 20;
    rc2[&c1] = nr;
    nr.most_likely_value = 8;
    rc2[&c2] = nr;

    auto avg = get_weighted_averages({ rc1, rc2 }, { .25, .75 });
    CHECK_EQ(17.5, avg[&c1]);
    CHECK_EQ(6.5, avg[&c2]);
}

class Reconstruction
{
public:
    unique_ptr<transcript_vector> p_transcripts;
    unique_ptr<clade> p_tree;

    Reconstruction()
    {
        gene_transcript fam("Family5", "", "");
        p_tree.reset(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

        fam.set_expression_value("A", pv::to_computational_space(11));
        fam.set_expression_value("B", pv::to_computational_space(2));
        fam.set_expression_value("C", pv::to_computational_space(5));
        fam.set_expression_value("D", pv::to_computational_space(6));

        p_transcripts.reset(new transcript_vector{fam});

    }
};

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction print_reconstructed_states")
{
    gamma_model_reconstruction gmr(*p_transcripts, discretized_gamma(0.5, 1));

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()].most_likely_value = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")].most_likely_value = 0;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")].most_likely_value = 0;

    rec.reconstruction[p_tree.get()] = pv::to_computational_space(7);
    rec.reconstruction[p_tree->find_descendant("AB")] = pv::to_computational_space(8);
    rec.reconstruction[p_tree->find_descendant("CD")] = pv::to_computational_space(6);

    ostringstream ost;
    gmr.print_reconstructed_states(ost, p_tree.get());
    CHECK_STREAM_CONTAINS(ost, "  TREE Family5 = ((A<1>_11:1,B<2>_2:3)<6>_8:7,(C<3>_5:11,D<4>_6:17)<7>_6:23)<5>_7;");
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__print_additional_data__prints_likelihoods")
{
    gamma_model_reconstruction gmr(*p_transcripts, discretized_gamma(0.5, 3));
    gmr._reconstructions["Family5"]._category_likelihoods = { 0.01, 0.03, 0.09, 0.07 };
    ostringstream ost;
    gmr.print_category_likelihoods(ost);
    CHECK_STREAM_CONTAINS(ost, "Transcript ID\t0.0550773\t0.565867\t2.37906\t\n");
    CHECK_STREAM_CONTAINS(ost, "Family5\t0.01\t0.03\0.09\t0.07");
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__prints_sigma_multipiers")
{
    gamma_model_reconstruction gmr(*p_transcripts, discretized_gamma(0.5, 2));

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()].most_likely_value = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")].most_likely_value = 8;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")].most_likely_value = 6;

    rec.reconstruction[p_tree.get()] = 7;
    rec.reconstruction[p_tree->find_descendant("AB")] = 8;
    rec.reconstruction[p_tree->find_descendant("CD")] = 6;

    std::ostringstream ost;
    gmr.print_reconstructed_states(ost, p_tree.get());

    CHECK_STREAM_CONTAINS(ost, "BEGIN SIGMA_MULTIPLIERS;");
    CHECK_STREAM_CONTAINS(ost, "  0.142516;");
    CHECK_STREAM_CONTAINS(ost, "  1.85748;");
    CHECK_STREAM_CONTAINS(ost, "END;");
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction get_internal_node_value returns reconstruction value for internal nodes")
{
    auto node = p_tree->find_descendant("CD");
    gamma_model_reconstruction gmr(*p_transcripts, discretized_gamma(0.5, 1));
    gmr._reconstructions["Family5"].reconstruction[node] = 7;

    CHECK_EQ(7, gmr.get_internal_node_value(p_transcripts->at(0), node).most_likely_value);

}

TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    vector<double> multipliers({ .2, .75 });
    vector<double> em;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    gf.set_expression_value("B", 2);
    transcript_vector transcripts{ gf };
    gamma_model_reconstruction gmr(transcripts, discretized_gamma(0.5, 2));

    gmr._reconstructions["myid"].reconstruction[p_tree->find_descendant("AB")] = 5;

    ostringstream ost;
    gmr.print_increases_decreases_by_clade(ost, p_tree.get(), true);
    CHECK_STREAM_CONTAINS(ost, "#Taxon_ID\tIncrease\tDecrease");
    CHECK_STREAM_CONTAINS(ost, "A<1>\t1\t0");
    CHECK_STREAM_CONTAINS(ost, "B<2>\t0\t1");
}

TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_clade empty")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;

    vector<double> multipliers({ .2, .75 });
    vector<double> em;

    transcript_vector transcripts;
    //gamma_model_reconstruction gmr(transcripts, em);
    gamma_model_reconstruction gmr(transcripts, discretized_gamma(0.5, 2));

    gmr.print_increases_decreases_by_clade(empty, p_tree.get(), true);
    CHECK_EQ(empty.str(), "#Taxon_ID\tIncrease\tDecrease\n");
}

#if 0
TEST_CASE("gamma_model_prune_returns_false_if_saturated" * doctest::skip(true))
{
    vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    sigma_squared lambda(0.9);
    matrix_cache cache;

    vector<double> cat_likelihoods;

    gamma_model model(&lambda, &families, { 1.0,1.0 }, { 0.1, 0.5 }, NULL);

    CHECK(!model.prune(families[0], new prior("gamma", 0.0, 1600), cache, &lambda, p_tree.get(), cat_likelihoods, boundaries(0, 20)));
}
#endif

inline void CHECK_SIGMA_VALUE(double val, const sigma_squared& sigma)
{
    CHECK_EQ(doctest::Approx(val), sigma.get_values()[0]);
}

TEST_CASE("Check sigma multipliers for a given alpha")
{
    vector<double> probs(3);
    vector<double> multipliers(3);
    double alpha = 0.635735;
    sigma_squared sigma((double)0.613693);
    get_gamma(probs, multipliers, alpha);
    CHECK_EQ(doctest::Approx(0.0976623), multipliers[0]);
    CHECK_EQ(doctest::Approx(0.653525), multipliers[1]);
    CHECK_EQ(doctest::Approx(2.24881), multipliers[2]);
    vector<sigma_squared> sigmas;
    for (auto m : multipliers)
    {
        sigmas.emplace_back(sigma, m);
    }
    REQUIRE_EQ(3, sigmas.size());
    CHECK_SIGMA_VALUE(0.0599347, sigmas[0]);
    CHECK_SIGMA_VALUE(0.401064, sigmas[1]);
    CHECK_SIGMA_VALUE(1.38008, sigmas[2]);
}

TEST_CASE("gamma_model_prune" * doctest::skip(true))
{
    vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 6);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    sigma_squared lambda(0.005);
    matrix_cache cache;

    gamma_model model(&lambda, &families, 1, 3, NULL);

    vector<double> cat_likelihoods;
    CHECK(model.prune(families[0], new prior("gamma", 0.0, 1600), cache, &lambda, p_tree.get(), cat_likelihoods, boundaries(0, 20)));

    CHECK_EQ(2, cat_likelihoods.size());
    CHECK_EQ(doctest::Approx(-23.04433), log(cat_likelihoods[0]));
    CHECK_EQ(doctest::Approx(-16.68005), log(cat_likelihoods[1]));
}

TEST_CASE("Simulation: get_simulation_sigma uses multiplier based on category probability")
{
    vector<double> gamma_categories{ 0.3, 0.7 };
    vector<double> multipliers{ 0.5, 1.5 };
    sigma_squared lam(0.05);
    gamma_model m(&lam, NULL, 3, 1.0, NULL);
    vector<double> results(100);
    generate(results.begin(), results.end(), [&m]() {
        unique_ptr<sigma_squared> new_lam(dynamic_cast<sigma_squared*>(m.get_simulation_sigma()));
        return new_lam->get_values()[0];
        });

    CHECK_EQ(doctest::Approx(0.0426), accumulate(results.begin(), results.end(), 0.0) / 100.0);

}
