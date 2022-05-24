#include <numeric>
#include <algorithm>
#include <fstream>
#include <random>
#include <vector>
#include <iomanip>

#include "doctest.h"
#include "easylogging++.h"

#include "simulator.h"
#include "user_data.h"
#include "core.h"
#include "matrix_cache.h"
#include "root_equilibrium_distribution.h"
#include "probability.h"
#include "sigma.h"
#include "DiffMat.h"
#include "arguments.h"
#include "newick_ape_loader.h"
#include "proportional_variance.h"

using namespace std;
using namespace Eigen;
namespace pv = proportional_variance;

extern std::mt19937 randomizer_engine; // seeding random number engine

class binner
{
    double _max_value;
    const int _Npts;
public:

    binner(int Npts, double max_value) : _Npts(Npts), _max_value(max_value) {}
    int bin(double value) const
    {
        return value / _max_value * (_Npts -1);
    }
    double value(int bin) const
    {
        double val = bin * _max_value / double(_Npts - 1);

        // subtract a small amount to deal with floating point inaccuracy
        // (The value was occasionally larger than max_value otherwise)
        return val > _max_value ? _max_value - MATRIX_EPSILON : val;

    }
    double max_value() const
    {
        return _max_value;
    }
};

simulator::simulator(user_data& d, const input_parameters& ui) : action(d, ui)
{
#ifdef SILENT
    quiet = true;
#endif
    _p_rootdist = create_rootdist(ui.rootdist_params_or_default(), d.rootdist);
}

void simulator::execute(std::vector<model *>& models)
{
    simulate(models, _user_input);
}

simulated_family create_simulated_family(const clade *p_tree, const sigma_squared* p_sigsqrd, int upper_bound, double root_value, const matrix_cache& cache)
{
    simulated_family sim;
    sim.sigma = p_sigsqrd->get_value_for_clade(p_tree);

    sim.values[p_tree] = root_value;
    boundaries bounds(pv::to_computational_space(0), upper_bound);
    std::function <void(const clade*)> get_child_value;
    get_child_value = [&](const clade* c) {
        auto m = cache.get_matrix(c->get_branch_length(), p_sigsqrd->get_value_for_clade(c), upper_bound);
        VectorXd v(m.cols());
        VectorPos_bounds(sim.values[c->get_parent()], bounds, v);
        VectorXd probs = m * v;
        std::discrete_distribution<int> distribution(probs.data(), probs.data() + probs.size());
        binner b(m.cols(), upper_bound);
        sim.values[c] = b.value(distribution(randomizer_engine));
        c->apply_to_descendants(get_child_value);
    };

    p_tree->apply_to_descendants(get_child_value);

    return sim;
}

std::vector<simulated_family> simulator::simulate_processes(model *p_model) {
    
    if (data.p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    int num_sims = _user_input.nsims;
    if (num_sims == 0)
    {
        num_sims = accumulate(data.rootdist.begin(), data.rootdist.end(), 0,
            [](int acc, std::pair<int, int> p) { return (acc + p.second); });
    }

    std::vector<simulated_family> results(num_sims);

    LOG(INFO) << "Simulating " << results.size() << " families for model " << p_model->name();

    vector<double> root_sizes(results.size());
    int n = 0;
    generate(root_sizes.begin(), root_sizes.end(), [this, &n]() mutable {
        return _p_rootdist->select_root_value(n++);
        });

    unique_ptr<sigma_squared> sim_sigsqd(p_model->get_simulation_lambda());
    unique_ptr<upper_bound_calculator> bound_calculator(upper_bound_calculator::create(sim_sigsqd.get(), data.p_tree));
    int upper_bound = bound_calculator->get_max_bound(root_sizes);
    LOG(DEBUG) << "Upper bound for discretization vector: " << upper_bound;

    matrix_cache cache;
    cache.precalculate_matrices(sim_sigsqd->get_values(), data.p_tree->get_branch_lengths(), upper_bound);

    transform(root_sizes.begin(), root_sizes.end(), results.begin(), [this, &sim_sigsqd, &cache, upper_bound](double root_size) {
        return create_simulated_family(data.p_tree, sim_sigsqd.get(), upper_bound, root_size, cache);
        });

    return results;
}

void print_header(std::ostream& ost, const input_parameters& p, size_t c, const clade *p_tree, const cladevector& order)
{
    auto tm = std::time(nullptr);

    ost << "# Simulated data set created " << std::put_time(std::localtime(&tm), "%Y-%m-%d %H:%M");
    ost << " (" << c << " transcripts)" << endl;
    ost << "# Root distribution: ";
    ost << p.rootdist_params_or_default() << endl;

    if (p_tree)
    {
        auto text_func = [&order](const clade* c) {
            return clade_index_or_name(c, order);
        };
        ost << "# Tree: ";
        p_tree->write_newick(ost, text_func);
        ost << endl;
    }

    ost << "# Sigma: ";
    if (p.fixed_multiple_lambdas.empty())
        ost << p.fixed_lambda << endl;
    else
        ost << p.fixed_multiple_lambdas << endl;
}

/// Simulate
/// \callgraph
void simulator::simulate(std::vector<model *>& models, const input_parameters &my_input_parameters)
{
    LOG(INFO) << "Simulating with " << models.size() << " model(s)";

    if (data.p_tree == nullptr)
        throw std::runtime_error("No tree specified for simulations");

    auto order = get_ape_order(data.p_tree);

    string dir = my_input_parameters.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);

    for (auto p_model : models) {

        auto results = simulate_processes(p_model);

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst2(fname);
        if (!ofst2) throw std::runtime_error("Failed to open " + fname);

        print_header(ofst2, my_input_parameters, results.size(), data.p_tree, order);
        print_simulations(ofst2, false, results, order);
        LOG(INFO) << "Simulated values written to " << fname;

        string truth_fname = filename("simulation_truth", dir);
        std::ofstream ofst(truth_fname);
        if (!ofst) throw std::runtime_error("Failed to open " + truth_fname);

        print_header(ofst, my_input_parameters, results.size(), data.p_tree, order);
        print_simulations(ofst, true, results, order);
        LOG(INFO) << "Simulated values (including internal nodes) written to " << truth_fname;

    }
}


void simulator::print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<simulated_family>& results, const cladevector& order) 
{
    if (results.empty())
    {
        LOG(ERROR) << "No simulations created" << endl;
        return;

    }
    ost << "DESC\tFID";
    for (size_t i = 0; i < order.size(); ++i)
    {
        if (!order[i]) continue;

        if (order[i]->is_leaf())
            ost << '\t' << order[i]->get_taxon_name();
        else if (include_internal_nodes)
            ost << '\t' << i;

    }
    ost << endl;

    for (size_t j = 0; j < results.size(); ++j) {
        auto& transcript = results[j];
        // Printing gene counts
        ost << "SIG" << transcript.sigma << "\ttranscript" << j;
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (!order[i]) continue;

            if (order[i]->is_leaf() || include_internal_nodes)
            {
                ost << '\t';
                ost << pv::to_user_space(transcript.values.at(order[i]));
            }
        }
        ost << endl;
    }
}

root_equilibrium_distribution* create_rootdist(std::string param, const vector<pair<float, int>>& rootdist)
{
    if (param.empty())
        throw std::runtime_error("No root distribution description provided");

    root_equilibrium_distribution* p_dist = nullptr;
    auto tokens = tokenize_str(param, ':');
    if (tokens[0] == "fixed")
    {
        auto fixed_value = stof(tokens[1]);
        LOG(INFO) << "Root distribution fixed at " << fixed_value;
        p_dist = new root_distribution_fixed(fixed_value);
    }
    else if (tokens[0] == "gamma")
    {
        auto k = stof(tokens[1]), theta = stof(tokens[2]);
        LOG(INFO) << "Using gamma root distribution with k=" << k << ", theta=" << theta << ")";
        p_dist = new root_distribution_gamma(k, theta);
    }
    else if (tokens[0] == "file")
    {
        p_dist = new root_distribution_specific(rootdist);
    }
    else
    {
        throw std::runtime_error("Couldn't parse root distribution");
    }

    return p_dist;
}

TEST_CASE("create_trial")
{
    randomizer_engine.seed(10);

    sigma_squared ss(0.25);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    const int ub = 5;

    matrix_cache cache;
    cache.precalculate_matrices(ss.get_values(), p_tree->get_branch_lengths(), ub);

    simulated_family actual = create_simulated_family(p_tree.get(), &ss, ub, 5.0, cache);

    CHECK_EQ(doctest::Approx(5.0), actual.values.at(p_tree.get()));
    CHECK_EQ(doctest::Approx(4.497487), actual.values.at(p_tree->find_descendant("A")));
    CHECK_EQ(doctest::Approx(4.42211), actual.values.at(p_tree->find_descendant("B")));
}

TEST_CASE("binner")
{
    binner b(200, 100);

    CHECK_EQ(5, b.bin(2.7));
    CHECK_EQ(doctest::Approx(60.3015), b.value(120));
    CHECK_EQ(25, b.bin(13));
    CHECK_EQ(doctest::Approx(20.1005), b.value(40));
}

TEST_CASE("binner unbins small values correctly")
{
    binner b(200, 20);
    CHECK_EQ(doctest::Approx(8.0402), b.value(80));
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("print_process_prints_in_order")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    vector<pair<string, double>> values{ {"A",2}, {"B",4},{"AB",6 } };

    std::ostringstream ost;
    clademap<double> t;
    for (auto v : values)
    {
        t[p_tree->find_descendant(v.first)] = pv::to_computational_space(v.second);
    }

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, true, my_trials, get_ape_order(p_tree.get()));

    CHECK_STREAM_CONTAINS(ost, "DESC\tFID\tA\tB\t3");
    CHECK_STREAM_CONTAINS(ost, "SIG0\ttranscript0\t2\t4\t6");

}

TEST_CASE("print_process_can_print_without_internal_nodes")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    vector<pair<string, double>> values{ {"A",2}, {"B",4},{"AB",6 } };

    std::ostringstream ost;
    clademap<double> t;
    for (auto v : values)
    {
        t[p_tree->find_descendant(v.first)] = pv::to_computational_space(v.second);
    }

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, false, my_trials, get_ape_order(p_tree.get()));
    CHECK_STREAM_CONTAINS(ost, "DESC\tFID\tA\tB\n");
    CHECK_STREAM_CONTAINS(ost, "SIG0\ttranscript0\t2\t4\n");

}

TEST_CASE("Check mean and variance of a simulated family leaf")
{
    randomizer_engine.seed(10);

    const double actual_sigma = 10;
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:1):1"));
    sigma_squared sim_sigsqd(actual_sigma);
    auto a = p_tree->find_descendant("A");
    const int ub = 15;

    matrix_cache cache;
    cache.precalculate_matrices(vector<double>{10}, p_tree->get_branch_lengths(), ub);
    size_t sz = 3;
    vector<double> v(sz);
    generate(v.begin(), v.end(), [&]() {
        auto sim = create_simulated_family(p_tree.get(), &sim_sigsqd, ub, 10, cache);
        return sim.values[a];
        });

    auto mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
    auto variance = std::accumulate(v.begin(), v.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));
     });

    CHECK_EQ(doctest::Approx(actual_sigma).epsilon(1), mean);
    CHECK_EQ(doctest::Approx(1.28027), variance);
}

TEST_CASE("print_header")
{
    std::ostringstream ost;
    input_parameters p;
    p.fixed_lambda = 2.5;
    p.rootdist_params = "Fixed:6.0";
    print_header(ost, p, 100, nullptr, cladevector());
    CHECK_STREAM_CONTAINS(ost, "# Simulated data set created ");
    CHECK_STREAM_CONTAINS(ost, "(100 transcripts)");
    CHECK_STREAM_CONTAINS(ost, "# Root distribution: Fixed:6.0");
    CHECK_STREAM_CONTAINS(ost, "# Sigma: 2.5");
}

TEST_CASE("print_header displays tree")
{
    std::ostringstream ost;
    input_parameters p;
    unique_ptr<clade> p_tree(parse_newick("((A:1,B:3):7,(C:11,D:17):23);"));

    print_header(ost, p, 100, p_tree.get(), get_ape_order(p_tree.get()));
    CHECK_STREAM_CONTAINS(ost, "# Tree: ((A<1>,B<2>)<6>,(C<3>,D<4>)<7>)<5>");
}

TEST_CASE("print_header default rootdist")
{
    std::ostringstream ost;
    input_parameters p;
    p.fixed_lambda = 2.5;
    print_header(ost, p, 100, nullptr, cladevector());
    CHECK_STREAM_CONTAINS(ost, "# Root distribution: gamma:0.375:1600.0");
}

TEST_CASE("print_header multiple sigmas")
{
    std::ostringstream ost;
    input_parameters p;
    p.fixed_multiple_lambdas = "1,2,3";
    print_header(ost, p, 100, nullptr, cladevector());
    CHECK_STREAM_CONTAINS(ost, "# Sigma: 1,2,3");
}

TEST_CASE("create_rootdist creates__specifed_distribution_if_given")
{
    randomizer_engine.seed(10);
    input_parameters params;
    unique_ptr<root_equilibrium_distribution> red(create_rootdist("file:rd.txt", { {2,11} }));
    REQUIRE(dynamic_cast<root_distribution_specific*>(red.get()));
    CHECK_EQ(doctest::Approx(pv::to_computational_space(2.0f)), red->select_root_value(0));
}

TEST_CASE("create_rootdist creates gamma distribution if given distribution")
{
    randomizer_engine.seed(10);
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("gamma:1:3", {}));
    auto gamma = dynamic_cast<root_distribution_gamma*>(rd.get());
    REQUIRE(gamma != nullptr);
    CHECK_EQ(doctest::Approx(1.87704f), gamma->get_raw_root_value(0));
}

TEST_CASE("create_rootdist throws error if nothing set")
{
    CHECK_THROWS_WITH(create_rootdist("", {}), "No root distribution description provided");
}

TEST_CASE("create_rootdist creates fixed root if requested")
{
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("fixed:6", {}));
    REQUIRE(dynamic_cast<root_distribution_fixed*>(rd.get()));
    CHECK_EQ(doctest::Approx(pv::to_computational_space(6.0)), rd->select_root_value(0));
}

class mock_model : public model {
    // Inherited via model
    virtual std::string name() const override { return "mockmodel"; }
    virtual void write_family_likelihoods(std::ostream& ost) override {}
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache* p_calc) override { return nullptr; }
    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<string>& sample_groups, const std::gamma_distribution<double>& prior) override { return nullptr; }
    bool _invalid_likelihood = false;
public:
    mock_model(sigma_squared* s) : model(s, NULL, NULL) {}
    virtual double infer_family_likelihoods(const user_data& ud, const sigma_squared* p_lambda, const std::gamma_distribution<double>& prior) override { return 0; }
};

TEST_CASE("Simulation, simulate_processes")
{
    sigma_squared lam(0.05);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    mock_model m(&lam);

    user_data ud;
    ud.p_tree = p_tree.get();
    ud.p_lambda = &lam;
    //    ud.p_prior = new root_distribution_uniform(size_t(100));

    input_parameters ip;
    ip.nsims = 100;
    simulator sim(ud, ip);

    auto sims = sim.simulate_processes(&m);
    CHECK_EQ(100, sims.size());
}

