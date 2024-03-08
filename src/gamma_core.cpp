#include <assert.h>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>

#include "doctest.h"
#include "easylogging++.h"

#include "gamma_core.h"
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
    model(p_sigma, p_gene_transcripts, p_error_model),
    _n_gamma_cats(n_gamma_cats) {

    _gamma = discretized_gamma(fixed_alpha, n_gamma_cats);
}

void gamma_model::write_extra_vital_statistics(std::ostream& ost)
{
    ost << "Alpha: " << get_alpha() << endl;
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
    _gamma = discretized_gamma(alpha, _n_gamma_cats);
}

void gamma_model::write_probabilities(ostream& ost)
{
    ost << "Alpha: " << get_alpha() << endl;
    ost << "Gamma cat probs are: ";
    _gamma.write_probabilities(ost);
    ost << endl << "Sigma multipliers are: ";
    _gamma.write_multipliers(ost, true);
    ost << endl;
}

sigma_squared* gamma_model::get_simulation_sigma()
{
    return _gamma.get_random_sigma(*_p_sigma);
}

bool gamma_model::can_infer() const
{
    if (!_p_sigma->is_valid())
        return false;

    if (get_alpha() < 0)
        return false;

    return true;
}

void flatten_vector(const vector<vector<double>>& v, vector<double>& result)
{
    for (auto& vv : v)
    {
        result.insert(result.end(), vv.begin(), vv.end());
    }
}

class likelihood_combiner
{
    typedef std::vector<double> v;
    typedef const sigma_squared *k;

    std::map<k,v> _likelihoods;
    int _num_transcripts;
public:
    likelihood_combiner(int num_transcripts) : _num_transcripts(num_transcripts){
    }

    void add_sigma(k key)
    {
        if (_likelihoods.find(key) == _likelihoods.end())
        {
            _likelihoods[key] = vector<double>(_num_transcripts);
        }
    }
    void set_value(k key, int transcript_index, double value)
    {
        _likelihoods[key][transcript_index] = value;
    }

    vector<k> get_sigmas() const {
        
        vector<k> keys(_likelihoods.size());
        std::transform(_likelihoods.begin(), _likelihoods.end(), keys.begin(),
               [](const std::pair<k,v>& kv) { return kv.first; });

        sort(keys.begin(), keys.end(), [](k a, k b) { return a->get_values()[0] < b->get_values()[0]; });
        return keys;
    }

    double get_value(k key, int transcript_index) const
    {
        return _likelihoods.at(key)[transcript_index];
    }

    double final_value()
    {
        // For each transcript, we add the likelihoods of all the categories, then take the log of that and add them all together
        auto keys = get_sigmas();
        double final_likelihood = 0.0;
        for (int i = 0; i<_num_transcripts; ++i)
        {
            final_likelihood += log(accumulate(keys.begin(), keys.end(), 0.0, [&](double a, k b) { return a + _likelihoods[b][i]; }));
        }
        return -final_likelihood;
    }
};

vector<double> get_priors(const matrix_cache& calc, const user_data& ud)
{
    vector<double> priors(calc.create_vector().size());
    for (size_t j = 0; j < priors.size(); ++j) {
        double x = (double(j) + 0.5) * double(ud.bounds.second) / (priors.size() - 1);

        priors[j] = computational_space_prior(x, ud.p_prior);
    }

    return priors;
}

gamma_model::~gamma_model()
{

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

    auto sigmas = _gamma.get_discrete_sigmas(*p_sigma);

    _p_all_transcripts_likelihood.reset(new likelihood_combiner(ud.gene_transcripts.size()));

    for (auto& s: sigmas)
    {    
        matrix_cache calc;
        vector<double> priors;
        try
        {
            priors = get_priors(calc, ud);
        }
        catch (std::domain_error& e)
        {
            LOG(DEBUG) << e.what();
            LOG(WARNING) << "Prior not valid for this sigma and data set";
            return -log(0);
        }

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

        _p_all_transcripts_likelihood->add_sigma(&s);
    #pragma omp parallel for
        for (int i = 0; i < (int)ud.gene_transcripts.size(); ++i) {
            _p_all_transcripts_likelihood->set_value(&s, i, compute_prior_likelihood(partial_likelihoods[references[i]], priors));            
        }
    }

    double final_likelihood = _p_all_transcripts_likelihood->final_value();
    LOG(INFO) << "Score (-lnL): " << std::setw(15) << std::setprecision(14) << final_likelihood;

    return final_likelihood;

}

sigma_optimizer_scorer* gamma_model::get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups)
{
    bool estimate_sigma = data.p_sigma == NULL;
    bool estimate_alpha = get_alpha() <= 0.0;

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

clademap<double> get_weighted_averages(const std::vector<clademap<node_reconstruction>>& m, const discretized_gamma& gamma)
{
    cladevector nodes(m[0].size());
    std::transform(m[0].begin(), m[0].end(), nodes.begin(), [](std::pair<const clade *, node_reconstruction> v) { return v.first;  });

    clademap<double> result;
    for (auto node : nodes)
    {
        vector<double> probabilities(m.size());
        for (size_t i = 0; i<probabilities.size(); ++i)
        {
            probabilities[i] = double(m[i].at(node).most_likely_value);
        }
        auto weighted = gamma.weight(probabilities);
        result[node] = accumulate(weighted.begin(), weighted.end(), 0.0);
    }

    return result;
}

reconstruction* gamma_model::reconstruct_ancestral_states(const user_data& ud, matrix_cache *calc)
{
    LOG(INFO) << "Starting reconstruction processes for Gamma model";

    vector<double> all;
    auto ss = _gamma.get_discrete_sigmas(*_p_sigma);
    for (auto& s : ss)
    {
        for (auto t : s.get_values())
        {
            all.push_back(t);
        }
    }

    calc->precalculate_matrices(all, ud.p_tree->get_branch_lengths(), ud.bounds);

    gamma_model_reconstruction* result = new gamma_model_reconstruction(ud.gene_transcripts, ud.p_replicate_model, _gamma);
    vector<gamma_model_reconstruction::gamma_reconstruction *> recs(ud.gene_transcripts.size());

    auto sigmas = _p_all_transcripts_likelihood->get_sigmas();
    for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
    {
        auto id = ud.gene_transcripts[i].id();
        recs[i] = &result->_reconstructions[id];
        auto& cl = result->_reconstructions[id]._category_likelihoods;
        cl.resize(_n_gamma_cats);
        transform(sigmas.begin(), sigmas.end(), cl.begin(), [&](const sigma_squared* s) {
            return _p_all_transcripts_likelihood->get_value(s, i);
        });
        result->_reconstructions[id].category_reconstruction.resize(_n_gamma_cats);
    }

    for (size_t k = 0; k<sigmas.size(); ++k)
    {
        VLOG(1) << "Reconstructing for multiplier " << sigmas[k];

        inference_pruner tr(sigmas[k], ud.p_tree, ud.p_replicate_model, calc, ud.bounds);

        for (size_t i = 0; i < ud.gene_transcripts.size(); ++i)
        {
            recs[i]->category_reconstruction[k] = tr.reconstruct(ud.gene_transcripts[i]);
        }
    }

    for (auto reconstruction : recs)
    {
        // multiply every reconstruction by gamma_cat_prob
        reconstruction->reconstruction = get_weighted_averages(reconstruction->category_reconstruction, _gamma);
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
        auto& cat_l = _reconstructions[gf.id()]._category_likelihoods;
        vector<double> log_likelihoods(cat_l.size());
        transform(cat_l.begin(), cat_l.end(), log_likelihoods.begin(), [](double x) { return log(x); });
        ostream_iterator<double> ct(ost, "\t");
        copy(log_likelihoods.begin(), log_likelihoods.end(), ct);
        ost << endl;
    }
}

void gamma_model_reconstruction::print_additional_data(std::string output_prefix)
{
    std::ofstream cat_likelihoods(filename("category_likelihoods", output_prefix));
    print_category_likelihoods(cat_likelihoods);

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

    discretized_gamma gamma(0.5, 2);
    auto avg = get_weighted_averages({ rc1, rc2 }, gamma);
    CHECK_EQ(15.0, avg[&c1]);
    CHECK_EQ(5.0, avg[&c2]);
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
    gmr._reconstructions["Family5"]._category_likelihoods = { exp(0.01), exp(0.03), exp(0.09), exp(0.07) };
    ostringstream ost;
    gmr.print_category_likelihoods(ost);
    CHECK_STREAM_CONTAINS(ost, "Transcript ID\t0.0442801\t0.454936\t1.91267\t\n");
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
    CHECK_STREAM_CONTAINS(ost, "  0.101531;");
    CHECK_STREAM_CONTAINS(ost, "  1.3233;");
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

TEST_CASE("likelihood_combiner")
{
    sigma_squared s1(0.1);
    likelihood_combiner lc(5);
    lc.add_sigma(&s1);
    const double e1 = exp(0.1);
    lc.set_value(&s1, 0, e1);
    lc.set_value(&s1, 1, e1);
    lc.set_value(&s1, 2, e1);
    lc.set_value(&s1, 3, e1);
    lc.set_value(&s1, 4, e1);
    CHECK_EQ(doctest::Approx(-0.5), lc.final_value());
}

TEST_CASE("likelihood_combiner sums logs of each value with a single sigma")
{
    sigma_squared s1(1.0);
    likelihood_combiner lc(5);
    lc.add_sigma(&s1);
    const double e1 = exp(.1);
    lc.set_value(&s1, 0, e1);
    lc.set_value(&s1, 1, e1);
    lc.set_value(&s1, 2, e1);
    lc.set_value(&s1, 3, e1);
    lc.set_value(&s1, 4, e1);
    CHECK_EQ(doctest::Approx(-0.5), lc.final_value());
}

TEST_CASE("likelihood_combiner get_value")
{
    sigma_squared s1(1.0);
    likelihood_combiner lc(5);
    lc.add_sigma(&s1);
    lc.set_value(&s1, 0, 1);
    lc.set_value(&s1, 1, 2);
    lc.set_value(&s1, 2, 3);
    CHECK_EQ(1, lc.get_value(&s1, 0));
    CHECK_EQ(2, lc.get_value(&s1, 1));
    CHECK_EQ(3, lc.get_value(&s1, 2));
}

TEST_CASE("likelihood_combiner get_sigmas")
{
    sigma_squared s1(1.0);
    likelihood_combiner lc(5);
    lc.add_sigma(&s1);
    lc.set_value(&s1, 0, 1);
    lc.set_value(&s1, 1, 2);
    lc.set_value(&s1, 2, 3);
    CHECK_EQ(1, lc.get_value(&s1, 0));
    CHECK_EQ(2, lc.get_value(&s1, 1));
    CHECK_EQ(3, lc.get_value(&s1, 2));
}

TEST_CASE("With multiple sigmas likelihood_combiner sums across sigmas then takes the logs")
{
    sigma_squared s1(1.0);
    sigma_squared s2(2.0);  // values don't matter, they are used as keys

    const double e1 = exp(.1) * .75; // cleverly choose some values that sum to a nice number
    const double e2 = exp(.1) * .25;

    const int num_transcripts = 5;
    likelihood_combiner lc(num_transcripts);
    lc.add_sigma(&s1);
    lc.set_value(&s1, 0, e1);
    lc.set_value(&s1, 1, e1);
    lc.set_value(&s1, 2, e1);
    lc.set_value(&s1, 3, e1);
    lc.set_value(&s1, 4, e1);

    lc.add_sigma(&s2);
    lc.set_value(&s2, 0, e2);
    lc.set_value(&s2, 1, e2);
    lc.set_value(&s2, 2, e2);
    lc.set_value(&s2, 3, e2);
    lc.set_value(&s2, 4, e2);

    CHECK_EQ(doctest::Approx(-0.5), lc.final_value());
}

TEST_CASE("get_priors")
{
    sigma_squared s1(0.1);
    matrix_cache calc;
    user_data ud;
    ud.p_prior = new prior("gamma", 1.0, 1600);
    ud.bounds = boundaries(0, 20);
    auto priors = get_priors(calc, ud);
    REQUIRE_EQ(200, priors.size());
    CHECK_EQ(doctest::Approx(-7.32816), log(priors[0]));
    CHECK_EQ(doctest::Approx(-1.02372), log(priors[75]));
    CHECK_EQ(doctest::Approx(-182.551), log(priors[125]));
    CHECK_EQ(0, priors[199]);
}
