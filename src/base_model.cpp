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
#include "reconstruction.h"
#include "proportional_variance.h"
#include "inference_pruner.h"

using namespace std;
namespace pv = proportional_variance;

extern mt19937 randomizer_engine;

//! \ingroup base Base Model
class base_model_reconstruction : public reconstruction
{
public:

    base_model_reconstruction()
    {

    }

    std::map<std::string, clademap<node_reconstruction>> _reconstructions;

    node_reconstruction get_internal_node_value(const gene_transcript& transcript, const clade* c) const
    {
        if (_reconstructions.find(transcript.id()) == _reconstructions.end())
            throw runtime_error("Transcript " + transcript.id() + " not found in reconstruction");

        return _reconstructions.at(transcript.id()).at(c);
    }

};

base_model::base_model(sigma_squared* p_sigma, const vector<gene_transcript>* p_gene_transcripts,
    error_model *p_error_model) :
    model(p_sigma, p_gene_transcripts, p_error_model)
{

}

vector<size_t> build_reference_list(const vector<gene_transcript>& transcripts)
{
    const size_t invalid_index = std::numeric_limits<size_t>::max();

    vector<size_t> reff;
    const int num_transcripts = transcripts.size();
    reff.resize(num_transcripts, invalid_index);
    for (int i = 0; i < num_transcripts; ++i) {
        if (reff[i] != invalid_index) continue;

        reff[i] = i;

        for (int j = i + 1; j < num_transcripts; ++j) {
            if (reff[j] == invalid_index)
            {
                if (transcripts[i].species_size_match(transcripts[j]))
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

double base_model::infer_transcript_likelihoods(const user_data& ud, const sigma_squared *p_sigma, const gamma_distribution<double>& prior) {
    //TIMED_FUNC(timerObj);
    _monitor.Event_InferenceAttempt_Started();

    if (!p_sigma->is_valid())
    {
        _monitor.Event_InferenceAttempt_InvalidValues();
        return -log(0);
    }

    int upper_bound = upper_bound_from_transcript_values(ud.gene_transcripts);
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
    std::vector<double> all_transcripts_likelihood(ud.gene_transcripts.size());

    calc.precalculate_matrices(p_sigma->get_values(),  ud.p_tree->get_branch_lengths(), upper_bound);

    vector<vector<double>> partial_likelihoods(ud.gene_transcripts.size());

    vector<inference_pruner> pruners;
    pruners.reserve(ud.gene_transcripts.size());
    std::generate_n(std::back_inserter(pruners), ud.gene_transcripts.size(), [&]() {return inference_pruner(calc, p_sigma, _p_error_model, ud.p_tree, 1.0); });
#pragma omp parallel for
    for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i) {
        if ((int)references[i] == i)
            partial_likelihoods[i] = pruners[i].prune(ud.gene_transcripts.at(i), upper_bound);
    }

#pragma omp parallel for
    for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i) {

        all_transcripts_likelihood[i] = compute_prior_likelihood(partial_likelihoods[references[i]], priors);
        
    }
    double final_likelihood = -std::accumulate(all_transcripts_likelihood.begin(), all_transcripts_likelihood.end(), 0.0);

    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;

    return final_likelihood;
}

sigma_optimizer_scorer* base_model::get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior)
{
    if (data.p_sigma != NULL)  // already have a sigma, nothing we want to optimize
        return nullptr;

    _p_sigma = initialize_search_sigma(data.p_sigma_tree, sample_groups);

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

    int upper_bound = upper_bound_from_transcript_values(ud.gene_transcripts);

    p_calc->precalculate_matrices(_p_sigma->get_values(), ud.p_tree->get_branch_lengths(), upper_bound);

    for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
    {
        auto &rc = result->_reconstructions[ud.gene_transcripts[i].id()];
        ud.p_tree->apply_prefix_order([&rc](const clade* c) {
            rc[c].most_likely_value = 0;
            });
    }

    inference_pruner tr(_p_sigma, ud.p_tree, p_calc);

    for (size_t i = 0; i< ud.gene_transcripts.size(); ++i)
    {
        result->_reconstructions[ud.gene_transcripts[i].id()] = tr.reconstruct(ud.gene_transcripts[i], upper_bound);
    }

    LOG(INFO) << "Done!\n";

    return result;
}

sigma_squared* base_model::get_simulation_sigma()
{
    return _p_sigma->multiply(simulation_sigma_multiplier);
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
    ud.p_sigma = NULL;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    base_model model(ud.p_sigma, NULL, NULL);

    unique_ptr<sigma_optimizer_scorer> opt(model.get_sigma_optimizer(ud, vector<string>(), std::gamma_distribution<double>(1, 2)));
    auto guesses = opt->initial_guesses();
    CHECK_EQ(1, guesses.size());
    CHECK_EQ(0.5, guesses[0]);
    delete model.get_sigma();
}

TEST_CASE("infer_transcript_likelihoods")
{
    user_data ud;
    unique_ptr<clade> tree(parse_newick("(A:1,B:1);"));
    ud.p_tree = tree.get();

    auto& transcripts = ud.gene_transcripts;
    gene_transcript fam;
    fam.set_expression_value("A", 1);
    fam.set_expression_value("B", 2);
    transcripts.push_back(fam);
    gene_transcript fam2;
    fam2.set_expression_value("A", 2);
    fam2.set_expression_value("B", 1);
    transcripts.push_back(fam2);
    gene_transcript fam3;
    fam3.set_expression_value("A", 3);
    fam3.set_expression_value("B", 6);
    transcripts.push_back(fam3);
    gene_transcript fam4;
    fam4.set_expression_value("A", 6);
    fam4.set_expression_value("B", 3);
    transcripts.push_back(fam4);

    sigma_squared ss(0.01);

    base_model core(&ss, &transcripts, NULL);

    double multi = core.infer_transcript_likelihoods(ud, &ss, std::gamma_distribution<double>(1, 2));

#ifdef GSL_FOUND
    CHECK_EQ(doctest::Approx(114.708), multi);
#else
    CHECK_EQ(doctest::Approx(114.7321), multi);
#endif
}


TEST_CASE("infer_transcript_likelihoods returns invalid on invalid prior")
{
    sigma_squared ss(5.0);
    user_data ud;
    ud.p_sigma = NULL;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    base_model model(&ss, &ud.gene_transcripts, NULL);
    CHECK_EQ(-log(0),  model.infer_transcript_likelihoods(ud, &ss, gamma_distribution<double>(0.0, 1600)));
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

    values[p_tree.get()].most_likely_value = pv::to_computational_space(7);
    values[p_tree->find_descendant("AB")].most_likely_value = pv::to_computational_space(8);
    values[p_tree->find_descendant("CD")].most_likely_value = pv::to_computational_space(6);

    ostringstream ost;

    bmr.print_reconstructed_states(ost, { fam }, p_tree.get());
    CHECK_STREAM_CONTAINS(ost, "#nexus");
    CHECK_STREAM_CONTAINS(ost, "BEGIN TREES;");
    CHECK_STREAM_CONTAINS(ost, "  TREE Family5 = ((A<1>_11.000000:1,B<2>_2.000000:3)<6>_8.000000:7,(C<3>_5.000000:11,D<4>_6.000000:17)<7>_6.000000:23)<5>_7.000000;");
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
    bmr._reconstructions["myid"][ab].most_likely_value = pv::to_computational_space(3);
    bmr._reconstructions["myid"][abcd].most_likely_value = pv::to_computational_space(3);

    CHECK_EQ(doctest::Approx(1.0), bmr.get_difference_from_parent(gf, a));
    CHECK_EQ(doctest::Approx(-1.0), bmr.get_difference_from_parent(gf, b));
    CHECK_EQ(0, bmr.get_difference_from_parent(gf, ab));
}


TEST_CASE("base_model_reconstruction")
{
    user_data ud;
    ud.p_tree = parse_newick("(A:1,B:1);");
    ud.gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    ud.gene_transcripts[0].set_expression_value("A", 1);
    ud.gene_transcripts[0].set_expression_value("B", 2);

    sigma_squared sl(0.05);
    ud.p_sigma = &sl;

    std::vector<gene_transcript> families(1);
    families[0].set_expression_value("A", 3);
    families[0].set_expression_value("B", 4);
    boundaries b(0, 20);
    base_model model(&sl, &families, NULL);

    matrix_cache calc;
    root_distribution_uniform dist(size_t(10));

    std::unique_ptr<base_model_reconstruction> rec(dynamic_cast<base_model_reconstruction*>(model.reconstruct_ancestral_states(ud, &calc)));

    CHECK_EQ(1, rec->_reconstructions.size());

}

TEST_CASE("build_reference_list")
{
    std::string str = "Desc\tFamily ID\tA\tB\n"
        "\t (null)1\t5\t10\n"
        "\t (null)2\t5\t7\n"
        "\t (null)3\t5\t10\n"
        "\t (null)4\t5\t7\n";
    std::istringstream ist(str);
    std::vector<gene_transcript> gt;
    read_gene_transcripts(ist, NULL, gt);
    auto actual = build_reference_list(gt);
    vector<int> expected({ 0, 1, 0, 1 });
    CHECK_EQ(expected.size(), actual.size());

}