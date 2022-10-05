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
#include "rootdist_estimator.h"
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

using namespace std;
namespace pv = proportional_variance;

extern mt19937 randomizer_engine;

//! @brief Holds data for reconstructing a tree based on the Gamma model
//! \ingroup gamma
class gamma_model_reconstruction : public reconstruction
{
    const std::vector<double> _lambda_multipliers;
    virtual void write_nexus_extensions(std::ostream& ost) override;

public:
    gamma_model_reconstruction(const std::vector<double>& lambda_multipliers) :
        _lambda_multipliers(lambda_multipliers)
    {
    }

    void print_additional_data(transcript_vector& gene_families, std::string output_prefix) override;

    void print_category_likelihoods(std::ostream& ost, transcript_vector& gene_families);

    node_reconstruction get_internal_node_value(const gene_transcript& transcript, const clade* c) const;

    struct gamma_reconstruction {
        std::vector<clademap<node_reconstruction>> category_reconstruction;
        clademap<double> reconstruction;
        std::vector<double> _category_likelihoods;
    };

    std::map<std::string, gamma_reconstruction> _reconstructions;
};

gamma_model::gamma_model(sigma_squared* p_lambda, std::vector<gene_transcript>* p_gene_families, int n_gamma_cats, double fixed_alpha, error_model* p_error_model) :
    model(p_lambda, p_gene_families, p_error_model) {

    _gamma_cat_probs.resize(n_gamma_cats);
    _lambda_multipliers.resize(n_gamma_cats);
    if (p_gene_families)
        _category_likelihoods.resize(p_gene_families->size());
    set_alpha(fixed_alpha);
}

gamma_model::gamma_model(sigma_squared* p_lambda, std::vector<gene_transcript>* p_gene_families, std::vector<double> gamma_categories, std::vector<double> multipliers, error_model *p_error_model) :
    model(p_lambda, p_gene_families,  p_error_model)
{
    _gamma_cat_probs = gamma_categories;
    _lambda_multipliers = multipliers;
    if (p_gene_families)
        _category_likelihoods.resize(p_gene_families->size());
}

void gamma_model::write_extra_vital_statistics(std::ostream& ost)
{
    ost << "Alpha: " << _alpha << endl;
}

//! Set alpha for gamma distribution
void gamma_model::set_alpha(double alpha) {

    _alpha = alpha;
    if (_gamma_cat_probs.size() > 1)
        get_gamma(_gamma_cat_probs, _lambda_multipliers, alpha); // passing vectors by reference

}

string comma_separated(const std::vector<double>& items)
{
    string s;
    for (auto i : items)
        s += (s.empty() ? "" : ",") + to_string(i);
    return s;
}

void gamma_model::write_probabilities(ostream& ost)
{
    ost << "Gamma cat probs are: " << comma_separated(_gamma_cat_probs) << endl;
    ost << "Lambda multipliers are: " << comma_separated(_lambda_multipliers) << endl;
}

sigma_squared* gamma_model::get_simulation_lambda()
{
    discrete_distribution<int> dist(_gamma_cat_probs.begin(), _gamma_cat_probs.end());
    return _p_sigma->multiply(_lambda_multipliers[dist(randomizer_engine)]);
}

std::vector<double> gamma_model::get_posterior_probabilities(std::vector<double> cat_likelihoods)
{
    size_t process_count = cat_likelihoods.size();

    vector<double> numerators(process_count);
    transform(cat_likelihoods.begin(), cat_likelihoods.end(), _gamma_cat_probs.begin(), numerators.begin(), multiplies<double>());

    double denominator = accumulate(numerators.begin(), numerators.end(), 0.0);
    vector<double> posterior_probabilities(process_count);
    transform(numerators.begin(), numerators.end(), posterior_probabilities.begin(), bind2nd(divides<double>(), denominator));

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

bool gamma_model::prune(const gene_transcript& family, const std::gamma_distribution<double>& prior, const matrix_cache& diff_mat, const sigma_squared*p_lambda,
    const clade *p_tree, std::vector<double>& category_likelihoods, int upper_bound) 
{
    category_likelihoods.clear();

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        inference_pruner pruner(diff_mat, p_lambda, _p_error_model, p_tree, _lambda_multipliers[k]);
        auto partial_likelihood = pruner.prune(family, upper_bound);
        if (accumulate(partial_likelihood.begin(), partial_likelihood.end(), 0.0) == 0.0)
            return false;   // saturation

#ifdef MODEL_GENE_EXPRESSION_LOGS
        throw std::runtime_error("Log values not implemented for gamma yet");
#endif
        std::vector<double> full(partial_likelihood.size());
        for (size_t j = 0; j < partial_likelihood.size(); ++j) {
            double eq_freq = gammapdf(j, prior);
            full[j] = partial_likelihood[j] * eq_freq;
        }

#ifdef USE_MAX_PROBABILITY
        category_likelihoods.push_back(*max_element(full.begin(), full.end()) * _gamma_cat_probs[k]); // get max (CAFE's approach)
#else
        category_likelihoods.push_back(accumulate(full.begin(), full.end(), 0.0) * _gamma_cat_probs[k]); // sum over all sizes (Felsenstein's approach)
#endif
    }

    return true;
}

//! Infer bundle
double gamma_model::infer_family_likelihoods(const user_data& ud, const sigma_squared*p_sigma, const std::gamma_distribution<double>& prior) {

    _monitor.Event_InferenceAttempt_Started();

    int upper_bound = upper_bound_from_transcript_values(ud.gene_families);

    results.clear();

    if (!can_infer())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    using namespace std;

    vector<double> all_bundles_likelihood(ud.gene_families.size());

    vector<bool> failure(ud.gene_families.size());

    vector<vector<family_info_stash>> pruning_results(ud.gene_families.size());
    matrix_cache cache;
    vector<double> multipliers;
    for (auto multiplier : _lambda_multipliers)
    {
        unique_ptr<sigma_squared> mult(p_sigma->multiply(multiplier));
        auto values = mult->get_values();
        multipliers.insert(multipliers.end(), values.begin(), values.end());
    }
    cache.precalculate_matrices(multipliers, ud.p_tree->get_branch_lengths(), upper_bound);

#pragma omp parallel for
    for (size_t i = 0; i < ud.gene_families.size(); i++) {
        auto& cat_likelihoods = _category_likelihoods[i];

        if (prune(ud.gene_families.at(i), prior, cache, p_sigma, ud.p_tree, cat_likelihoods, upper_bound))
        {
            double family_likelihood = accumulate(cat_likelihoods.begin(), cat_likelihoods.end(), 0.0);

            vector<double> posterior_probabilities = get_posterior_probabilities(cat_likelihoods);

            pruning_results[i].resize(cat_likelihoods.size());
            for (size_t k = 0; k < cat_likelihoods.size(); ++k)
            {
                pruning_results[i][k] = family_info_stash(ud.gene_families.at(i).id(),_lambda_multipliers[k], cat_likelihoods[k],
                    family_likelihood, posterior_probabilities[k], posterior_probabilities[k] > 0.95);
                //            cout << "Bundle " << i << " Process " << k << " family likelihood = " << family_likelihood << endl;
            }
            all_bundles_likelihood[i] = std::log(family_likelihood);
        }
        else
        {
            // we got here because one of the gamma categories was saturated - reject this 
            failure[i] = true;
        }
    }

    if (find(failure.begin(), failure.end(), true) != failure.end())
    {
        for (size_t i = 0; i < ud.gene_families.size(); i++) {
            if (failure[i])
            {
                _monitor.Event_InferenceAttempt_Saturation(ud.gene_families.at(i).id());
            }
        }
        return -log(0);
    }
    for (auto& stashes : pruning_results)
    {
        for (auto& stash : stashes)
        {
            results.push_back(stash);
        }
    }
    double final_likelihood = -accumulate(all_bundles_likelihood.begin(), all_bundles_likelihood.end(), 0.0);

    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;
    return final_likelihood;
}

sigma_optimizer_scorer* gamma_model::get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior)
{
    bool estimate_sigma = data.p_lambda == NULL;
    bool estimate_alpha = _alpha <= 0.0;

    if (estimate_sigma && estimate_alpha)
    {
        _p_sigma = initialize_search_sigma(data.p_lambda_tree, sample_groups);
        return new sigma_optimizer_scorer(this, data, prior, _p_sigma);
    }
    else if (estimate_sigma && !estimate_alpha)
    {
        _p_sigma = initialize_search_sigma(data.p_lambda_tree, sample_groups);
        return new sigma_optimizer_scorer(dynamic_cast<model *>(this), data, prior, _p_sigma);
    }
    else if (!estimate_sigma && estimate_alpha)
    {
        _p_sigma = data.p_lambda->clone();
        return new sigma_optimizer_scorer(this, data, prior);
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

    int upper_bound = upper_bound_from_transcript_values(ud.gene_families);

    auto values = _p_sigma->get_values();
    vector<double> all;
    for (double multiplier : _lambda_multipliers)
    {
        for (double lambda : values)
        {
            all.push_back(lambda*multiplier);
        }
    }

    calc->precalculate_matrices(_p_sigma->get_values(), ud.p_tree->get_branch_lengths(), upper_bound);

    gamma_model_reconstruction* result = new gamma_model_reconstruction(_lambda_multipliers);
    vector<gamma_model_reconstruction::gamma_reconstruction *> recs(ud.gene_families.size());
    for (size_t i = 0; i < ud.gene_families.size(); ++i)
    {
        recs[i] = &result->_reconstructions[ud.gene_families[i].id()];
        result->_reconstructions[ud.gene_families[i].id()]._category_likelihoods = _category_likelihoods[i];
        result->_reconstructions[ud.gene_families[i].id()].category_reconstruction.resize(_lambda_multipliers.size());
    }

    for (size_t k = 0; k < _gamma_cat_probs.size(); ++k)
    {
        VLOG(1) << "Reconstructing for multiplier " << _lambda_multipliers[k];
        unique_ptr<sigma_squared> ml(_p_sigma->multiply(_lambda_multipliers[k]));

        inference_pruner tr(ml.get(), ud.p_tree, calc);

        for (size_t i = 0; i < ud.gene_families.size(); ++i)
        {
            recs[i]->category_reconstruction[k] = tr.reconstruct(ud.gene_families[i], upper_bound);
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
    ost << "\nBEGIN LAMBDA_MULTIPLIERS;\n";
    for (auto& lm : _lambda_multipliers)
    {
        ost << "  " << lm << ";\n";
    }
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

void gamma_model_reconstruction::print_category_likelihoods(std::ostream& ost, transcript_vector& transcripts)
{
    ost << "Family ID\t";
    ostream_iterator<double> lm(ost, "\t");
    copy(_lambda_multipliers.begin(), _lambda_multipliers.end(), lm);
    ost << endl;

    for (auto& gf : transcripts)
    {
        ost << gf.id() << '\t';
        auto rc = _reconstructions[gf.id()];
        ostream_iterator<double> ct(ost, "\t");
        copy(rc._category_likelihoods.begin(), rc._category_likelihoods.end(), ct);
        ost << endl;
    }
}

void gamma_model_reconstruction::print_additional_data(transcript_vector& gene_families, std::string output_prefix)
{
    std::ofstream cat_likelihoods(filename("category_likelihoods", output_prefix));
    print_category_likelihoods(cat_likelihoods, gene_families);

}

TEST_CASE("Inference: gamma_model__creates__lambda_optimizer__if_alpha_provided")
{
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:1):1,(C:1,D:1):1);"));

    gamma_model model(NULL, NULL, 4, 0.25, NULL);

    user_data data;
    data.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    data.gene_families[0].set_expression_value("A", 1);
    data.gene_families[0].set_expression_value("B", 2);
    data.p_tree = p_tree.get();

    auto opt = model.get_sigma_optimizer(data, vector<string>(), std::gamma_distribution<double>(1, 2));

    REQUIRE(opt);
    CHECK_EQ("Optimizing Sigma ", opt->description());
    delete model.get_sigma();
}

TEST_CASE("Inference: gamma_model__creates__gamma_optimizer__if_lambda_provided")
{
    gamma_model model(NULL, NULL, 4, -1, NULL);

    user_data data;

    sigma_squared sl(0.05);
    data.p_lambda = &sl;

    auto opt = model.get_sigma_optimizer(data, vector<string>(), std::gamma_distribution<double>(1, 2));

    REQUIRE(opt);
    CHECK_EQ("Optimizing Alpha ", opt->description());

    delete model.get_sigma();
}

TEST_CASE("Inference: gamma_model_creates_nothing_if_lambda_and_alpha_provided")
{
    gamma_model model(NULL, NULL, 4, .25, NULL);

    user_data data;

    sigma_squared sl(0.05);
    data.p_lambda = &sl;

    CHECK(model.get_sigma_optimizer(data, vector<string>(), std::gamma_distribution<double>(1, 2)) == nullptr);
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
    gene_transcript fam;
    unique_ptr<clade> p_tree;

    Reconstruction() : fam("Family5", "", "")
    {
        p_tree.reset(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

        fam.set_expression_value("A", pv::to_computational_space(11));
        fam.set_expression_value("B", pv::to_computational_space(2));
        fam.set_expression_value("C", pv::to_computational_space(5));
        fam.set_expression_value("D", pv::to_computational_space(6));
    }
};

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction print_reconstructed_states")
{
    gamma_model_reconstruction gmr(vector<double>({ 1.0 }));

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()].most_likely_value = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")].most_likely_value = 0;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")].most_likely_value = 0;

    rec.reconstruction[p_tree.get()] = pv::to_computational_space(7);
    rec.reconstruction[p_tree->find_descendant("AB")] = pv::to_computational_space(8);
    rec.reconstruction[p_tree->find_descendant("CD")] = pv::to_computational_space(6);

    ostringstream ost;
    gmr.print_reconstructed_states(ost, { fam }, p_tree.get(), 0.05);
    CHECK_STREAM_CONTAINS(ost, "  TREE Family5 = ((A<1>_11.000000:1,B<2>_2.000000:3)<6>_8.000000:7,(C<3>_5.000000:11,D<4>_6.000000:17)<7>_6.000000:23)<5>_7.000000;");
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__print_additional_data__prints_likelihoods")
{
    gamma_model_reconstruction gmr(vector<double>({ 0.3, 0.9, 1.4, 2.0 }));
    gmr._reconstructions["Family5"]._category_likelihoods = { 0.01, 0.03, 0.09, 0.07 };
    ostringstream ost;
    gmr.print_category_likelihoods(ost, { fam });
    CHECK_STREAM_CONTAINS(ost, "Family ID\t0.3\t0.9\t1.4\t2\t\n");
    CHECK_STREAM_CONTAINS(ost, "Family5\t0.01\t0.03\0.09\t0.07");
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction__prints_lambda_multipiers")
{
    vector<double> multipliers{ 0.13, 1.4 };
    gamma_model_reconstruction gmr(multipliers);

    auto& rec = gmr._reconstructions["Family5"];
    rec.category_reconstruction.resize(1);
    rec.category_reconstruction[0][p_tree.get()].most_likely_value = 7;
    rec.category_reconstruction[0][p_tree->find_descendant("AB")].most_likely_value = 8;
    rec.category_reconstruction[0][p_tree->find_descendant("CD")].most_likely_value = 6;

    rec.reconstruction[p_tree.get()] = 7;
    rec.reconstruction[p_tree->find_descendant("AB")] = 8;
    rec.reconstruction[p_tree->find_descendant("CD")] = 6;

    std::ostringstream ost;
    gmr.print_reconstructed_states(ost, { fam }, p_tree.get(), 0.05);

    CHECK_STREAM_CONTAINS(ost, "BEGIN LAMBDA_MULTIPLIERS;");
    CHECK_STREAM_CONTAINS(ost, "  0.13;");
    CHECK_STREAM_CONTAINS(ost, "  1.4;");
    CHECK_STREAM_CONTAINS(ost, "END;");
}

TEST_CASE_FIXTURE(Reconstruction, "gamma_model_reconstruction get_internal_node_value returns reconstruction value for internal nodes")
{
    auto node = p_tree->find_descendant("CD");
    gamma_model_reconstruction gmr({ .5 });
    gmr._reconstructions["Family5"].reconstruction[node] = 7;

    CHECK_EQ(7, gmr.get_internal_node_value(fam, node).most_likely_value);

}
TEST_CASE("Reconstruction: gamma_model_print_increases_decreases_by_clade")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    ostringstream empty;

    vector<double> multipliers({ .2, .75 });
    vector<double> em;
    gamma_model_reconstruction gmr(em);

    gmr.print_increases_decreases_by_clade(empty, p_tree.get(), {});
    CHECK_EQ(empty.str(), "#Taxon_ID\tIncrease\tDecrease\n");

    gmr._reconstructions["myid"].reconstruction[p_tree->find_descendant("AB")] = 5;

    gene_transcript gf("myid", "", "");
    gf.set_expression_value("A", 7);
    gf.set_expression_value("B", 2);

    ostringstream ost;
    gmr.print_increases_decreases_by_clade(ost, p_tree.get(), {gf});
    CHECK_STREAM_CONTAINS(ost, "#Taxon_ID\tIncrease\tDecrease");
    CHECK_STREAM_CONTAINS(ost, "A<1>\t1\t0");
    CHECK_STREAM_CONTAINS(ost, "B<2>\t0\t1");
}


