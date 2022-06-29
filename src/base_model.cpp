#include <cmath>
#include <numeric>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>
#include <iomanip>

#include "doctest.h"

#include "base_model.h"
#include "rootdist_estimator.h"
#include "matrix_cache.h"
#include "gene_transcript.h"
#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "optimizer_scorer.h"
#include "simulator.h"
#include "sigma.h"
#include "io.h"
#include "DiffMat.h"
#include "transcript_reconstructor.h"
#include "proportional_variance.h"

using namespace std;
namespace pv = proportional_variance;

extern mt19937 randomizer_engine;

base_model::base_model(sigma_squared* p_lambda, const vector<gene_transcript>* p_gene_families,
    error_model *p_error_model) :
    model(p_lambda, p_gene_families, p_error_model)
{

}

vector<size_t> build_reference_list(const vector<gene_transcript>& families)
{
    const size_t invalid_index = std::numeric_limits<size_t>::max();

    vector<size_t> reff;
    const int num_families = families.size();
    reff.resize(num_families, invalid_index);
    for (int i = 0; i < num_families; ++i) {
        if (reff[i] != invalid_index) continue;

        reff[i] = i;

        for (int j = i + 1; j < num_families; ++j) {
            if (reff[j] == invalid_index)
            {
                if (families[i].species_size_match(families[j]))
                {
                    reff[j] = i;
                }
            }
        }
    }

    return reff;
}


inline double computational_space_prior(double val, const gamma_distribution<double>& prior)
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    return exp(val) * gammapdf(exp(val), prior);
#else
    return gammapdf(val, prior);
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

double base_model::infer_family_likelihoods(const user_data& ud, const sigma_squared *p_sigma, const gamma_distribution<double>& prior) {
    //TIMED_FUNC(timerObj);
    _monitor.Event_InferenceAttempt_Started();

    if (!p_sigma->is_valid())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    int upper_bound = upper_bound_from_transcript_values(ud.gene_families);
    matrix_cache calc;
    auto v = calc.create_vector();
    vector<double> priors(v.size());
    copy(v.begin(), v.end(), priors.begin());
    for (size_t j = 0; j < priors.size(); ++j) {
        double x = (double(j) + 0.5) * double(upper_bound) / (priors.size() - 1);

        priors[j] = computational_space_prior(x, prior);
    }

    if (all_of(priors.begin(), priors.end(), [](double prior) { return prior <= 0 || isnan(prior); }))
    {
        LOG(WARNING) << "Prior not valid for this sigma and data set";
        return -log(0);
    }
    results.resize(ud.gene_families.size());
    std::vector<double> all_families_likelihood(ud.gene_families.size());

    calc.precalculate_matrices(p_sigma->get_values(),  ud.p_tree->get_branch_lengths(), upper_bound);

    vector<vector<double>> partial_likelihoods(ud.gene_families.size());

    vector<inference_pruner> pruners;
    pruners.reserve(ud.gene_families.size());
    std::generate_n(std::back_inserter(pruners), ud.gene_families.size(), [&]() {return inference_pruner(calc, p_sigma, _p_error_model, ud.p_tree, 1.0); });
#pragma omp parallel for
    for (int i = 0; i < ud.gene_families.size(); ++i) {
        if (references[i] == i)
            partial_likelihoods[i] = pruners[i].prune(ud.gene_families.at(i), upper_bound);
            // probabilities of various family sizes
    }

    // prune all the families with the same lambda
#pragma omp parallel for
    for (int i = 0; i < ud.gene_families.size(); ++i) {

        all_families_likelihood[i] = compute_prior_likelihood(partial_likelihoods[references[i]], priors);
        
        results[i] = family_info_stash(ud.gene_families.at(i).id(), 0.0, 0.0, 0.0, all_families_likelihood[i], false);
    }
    double final_likelihood = -std::accumulate(all_families_likelihood.begin(), all_families_likelihood.end(), 0.0); // sum over all families

    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;

    return final_likelihood;
}

void base_model::write_family_likelihoods(std::ostream& ost)
{
    ost << "#FamilyID\tLikelihood of Family" << endl;
    for (const auto& r : results)
    {
        ost << r.family_id << "\t" << r.posterior_probability << endl;
    }
}

sigma_optimizer_scorer* base_model::get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior)
{
    if (data.p_lambda != NULL)  // already have a lambda, nothing we want to optimize
        return nullptr;

    _p_sigma = initialize_search_sigma(data.p_lambda_tree, sample_groups);

    if (_p_error_model && !data.p_error_model)
    {
        return new sigma_optimizer_scorer(this, data, prior, _p_sigma, _p_error_model);
    }
    else
    {
        return new sigma_optimizer_scorer(this, data, prior, _p_sigma);
    }
}

#define EPSILON_RANGES

reconstruction* base_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc)
{
    LOG(INFO) << "Starting reconstruction processes for Base model";

    auto result = new base_model_reconstruction();

    int upper_bound = upper_bound_from_transcript_values(ud.gene_families);

    p_calc->precalculate_matrices(_p_sigma->get_values(), ud.p_tree->get_branch_lengths(), upper_bound);

    for (size_t i = 0; i < ud.gene_families.size(); ++i)
    {
        auto &rc = result->_reconstructions[ud.gene_families[i].id()];
        ud.p_tree->apply_prefix_order([&rc](const clade* c) {
            rc[c] = 0;
            });
    }

    transcript_reconstructor tr(_p_sigma, ud.p_tree, p_calc);

    for (int i = 0; i< ud.gene_families.size(); ++i)
    {
        result->_reconstructions[ud.gene_families[i].id()] = tr.reconstruct_gene_transcript(ud.gene_families[i], upper_bound);
    }

    LOG(INFO) << "Done!\n";

    return result;
}

sigma_squared* base_model::get_simulation_lambda()
{
    return _p_sigma->multiply(simulation_lambda_multiplier);
}

double base_model_reconstruction::get_node_value(const gene_transcript& family, const clade *c) const
{
    if (c->is_leaf())
        return family.get_expression_value(c->get_taxon_name());

    if (_reconstructions.find(family.id()) == _reconstructions.end())
        throw runtime_error("Family " + family.id() + " not found in reconstruction");

    return _reconstructions.at(family.id()).at(c);
}

TEST_CASE("compute_prior_likelihood combines prior and inference correctly")
{
    gene_transcript gt;
    gt.set_expression_value("A", 12);
    gt.set_expression_value("B", 24);

    vector<double> inf{ 0.1, 0.2, 0.3};

    vector<double> priors({ 1.43078e-15,    2.5363e-23,  5.65526e-35 });
    double actual = compute_prior_likelihood(inf, priors);

#ifdef USE_MAX_PROBABILITY
    CHECK_EQ(doctest::Approx(-35.7683), actual);
#else
    CHECK_EQ(doctest::Approx(-36.4831), actual);
#endif
}

TEST_CASE("base optimizer guesses sigma only")
{
    user_data ud;
    ud.p_lambda = NULL;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    base_model model(ud.p_lambda, NULL, NULL);

    unique_ptr<sigma_optimizer_scorer> opt(model.get_sigma_optimizer(ud, vector<string>(), std::gamma_distribution<double>(1, 2)));
    auto guesses = opt->initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(0.5, guesses[0]);
    delete model.get_sigma();
}

TEST_CASE("infer_family_likelihoods")
{
    user_data ud;
    unique_ptr<clade> tree(parse_newick("(A:1,B:1);"));
    ud.p_tree = tree.get();

    auto& families = ud.gene_families;
    gene_transcript fam;
    fam.set_expression_value("A", 1);
    fam.set_expression_value("B", 2);
    families.push_back(fam);
    gene_transcript fam2;
    fam2.set_expression_value("A", 2);
    fam2.set_expression_value("B", 1);
    families.push_back(fam2);
    gene_transcript fam3;
    fam3.set_expression_value("A", 3);
    fam3.set_expression_value("B", 6);
    families.push_back(fam3);
    gene_transcript fam4;
    fam4.set_expression_value("A", 6);
    fam4.set_expression_value("B", 3);
    families.push_back(fam4);

    sigma_squared ss(0.01);

    base_model core(&ss, &families, NULL);

    double multi = core.infer_family_likelihoods(ud, &ss, std::gamma_distribution<double>(1, 2));

    CHECK_EQ(doctest::Approx(114.7321), multi);
}


TEST_CASE("infer_family_likelihoods returns invalid on invalid prior")
{
    sigma_squared ss(5.0);
    user_data ud;
    ud.p_lambda = NULL;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_families.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_families[0].set_expression_value("A", 1);
    ud.gene_families[0].set_expression_value("B", 2);

    base_model model(&ss, &ud.gene_families, NULL);
    CHECK_EQ(-log(0),  model.infer_family_likelihoods(ud, &ss, gamma_distribution<double>(0.0, 1600)));
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

        fam.set_expression_value("A", pv::to_computational_space(11));
        fam.set_expression_value("B", pv::to_computational_space(2));
        fam.set_expression_value("C", pv::to_computational_space(5));
        fam.set_expression_value("D", pv::to_computational_space(6));

        vector<string> nodes{ "A", "B", "C", "D", "AB", "CD", "ABCD" };
        order.resize(nodes.size());
        const clade* t = p_tree.get();
        transform(nodes.begin(), nodes.end(), order.begin(), [t](string s) { return t->find_descendant(s); });

    }
};

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE_FIXTURE(Reconstruction, "base_model_reconstruction__print_reconstructed_states")
{
    base_model_reconstruction bmr;
    auto& values = bmr._reconstructions["Family5"];

    values[p_tree.get()] = pv::to_computational_space(7);
    values[p_tree->find_descendant("AB")] = pv::to_computational_space(8);
    values[p_tree->find_descendant("CD")] = pv::to_computational_space(6);

    branch_probabilities branch_probs;

    ostringstream ost;

    bmr.print_reconstructed_states(ost, order, { fam }, p_tree.get(), 0.05, branch_probs);
    CHECK_STREAM_CONTAINS(ost, "#nexus");
    CHECK_STREAM_CONTAINS(ost, "BEGIN TREES;");
    CHECK_STREAM_CONTAINS(ost, "  TREE Family5 = ((A<0>_11.000000:1,B<1>_2.000000:3)<4>_8.000000:7,(C<2>_5.000000:11,D<3>_6.000000:17)<5>_6.000000:23)<6>_7.000000;");
    CHECK_STREAM_CONTAINS(ost, "END;");

}

TEST_CASE("increase_decrease")
{
    base_model_reconstruction bmr;
    gene_transcript gf("myid", "", "");

    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

    auto a = p_tree->find_descendant("A");
    auto b = p_tree->find_descendant("B");
    auto ab = p_tree->find_descendant("AB");
    auto abcd = p_tree->find_descendant("ABCD");

    gf.set_expression_value("A", pv::to_computational_space(4));
    gf.set_expression_value("B", pv::to_computational_space(2));
    bmr._reconstructions["myid"][ab] = pv::to_computational_space(3);
    bmr._reconstructions["myid"][abcd] = pv::to_computational_space(3);

    CHECK_EQ(doctest::Approx(1.0), bmr.get_difference_from_parent(gf, a));
    CHECK_EQ(doctest::Approx(-1.0), bmr.get_difference_from_parent(gf, b));
    CHECK_EQ(0, bmr.get_difference_from_parent(gf, ab));
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities" * doctest::skip(true))
{
    sigma_squared lm(0.05);
    matrix_cache cache;
    cache.precalculate_matrices(lm.get_values(), { 1,3,7 }, 20);
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree->find_descendant("AB")] = 10;
    rec._reconstructions[fam.id()][p_tree->find_descendant("ABCD")] = 12;
    CHECK_EQ(doctest::Approx(0.537681), compute_viterbi_sum(p_tree->find_descendant("AB"), fam, &rec, cache, &lm, 20)._value);
}

TEST_CASE_FIXTURE(Reconstruction, "viterbi_sum_probabilities_returns_invalid_if_root")
{
    sigma_squared lm(0.05);
    matrix_cache cache;
    cache.precalculate_matrices(lm.get_values(), { 1,3,7 }, 20);
    base_model_reconstruction rec;
    rec._reconstructions[fam.id()][p_tree.get()] = 11;
    CHECK_FALSE(compute_viterbi_sum(p_tree.get(), fam, &rec, cache, &lm, 20)._is_valid);
}







