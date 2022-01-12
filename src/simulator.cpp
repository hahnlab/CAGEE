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

using namespace std;
using namespace Eigen;
using root_distribution = std::map<int, float>;

extern std::mt19937 randomizer_engine; // seeding random number engine


class binner
{
    double _max_value;
public:
    binner(const sigma* p_lambda, const clade* p_tree, double root_size)
    {
        double t = p_tree->distance_from_root_to_tip();
        double sigma = p_lambda->get_value_for_clade(p_tree);
        _max_value = double(root_size) + 4.5 * sigma * sqrt(t);
        VLOG(SIMULATOR) << "Root size: " << root_size << " => max value: " << _max_value << " (Tree length: " << t << ", Sigma: " << sigma << ")";
    }

    int bin(double value) const
    {
        return value / _max_value * (DISCRETIZATION_RANGE-1);
    }
    double value(int bin) const
    {
        return bin * _max_value / double(DISCRETIZATION_RANGE-1);
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
    _p_rootdist = create_rootdist(ui.rootdist_params, d.rootdist);
}

void simulator::execute(std::vector<model *>& models)
{
    simulate(models, _user_input);
}

simulated_family create_simulated_family(const clade *p_tree, const sigma* p_sigma, double root_value)
{
    simulated_family sim;
    sim.lambda = p_sigma->get_value_for_clade(p_tree);

    binner b(p_sigma, p_tree, root_value);

    sim.values[p_tree] = root_value;

    std::function <void(const clade*)> get_child_value;
    get_child_value = [&](const clade* c) {
        double sigma = p_sigma->get_value_for_clade(c);
        MatrixXd m = ConvProp_bounds(c->get_branch_length(), sigma * sigma / 2, DiffMat::instance(), boundaries(0.0, b.max_value()));
        VectorXd v = VectorPos_bounds(sim.values[c->get_parent()], DISCRETIZATION_RANGE, boundaries(0, b.max_value()));
        VectorXd probs = m * v;
        std::discrete_distribution<int> distribution(probs.data(), probs.data() + probs.size());
        sim.values[c] = b.value(distribution(randomizer_engine));
        c->apply_to_descendants(get_child_value);
    };

    p_tree->apply_to_descendants(get_child_value);

    return sim;
}

// At the root, we have a vector of length DISCRETIZATION_RANGE. This has probability 1 at the size of the root
// and 0 everywhere else
// for each child, generate the transition matrix and multiply
simulated_family simulator::create_trial(const sigma*p_sigma, int family_number) {
    double root_size = _p_rootdist->select_root_value(family_number);

    if (data.p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    return create_simulated_family(data.p_tree, p_sigma, root_size);
}

void simulator::simulate_processes(model *p_model, std::vector<simulated_family>& results) {

    if (_user_input.nsims > 0)
    {
        results.resize(_user_input.nsims);
    }
    else
    {
        results.resize(accumulate(data.rootdist.begin(), data.rootdist.end(), 0,
            [](int acc, std::pair<int, int> p) { return (acc + p.second); }));
    }

    LOG(INFO) << "Simulating " << results.size() << " families for model " << p_model->name();

    for (size_t i = 0; i < results.size(); i+= LAMBDA_PERTURBATION_STEP_SIZE)
    {
        unique_ptr<sigma> sim_lambda(p_model->get_simulation_lambda());

        int n = 0;

        auto end_it = i + LAMBDA_PERTURBATION_STEP_SIZE > results.size() ? results.end() : results.begin() + i + LAMBDA_PERTURBATION_STEP_SIZE;
        generate(results.begin()+i, end_it, [this, &sim_lambda, i, &n]() mutable {
            return create_trial(sim_lambda.get(), i+n++);
        });
    }
}

void print_header(std::ostream& ost, const input_parameters& p, size_t c)
{
    auto tm = std::time(nullptr);

    ost << "# Simulated data set created " << std::put_time(std::localtime(&tm), "%Y-%m-%d %H:%M");
    ost << " (" << c << " transcripts)" << endl;
    ost << "# Root distribution: ";
    if (p.rootdist_params.empty())
        ost << "gamma:0.75:0.033" << endl;
    else
        ost << p.rootdist_params << endl;

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

    std::vector<const clade *> order;
    for_each(data.p_tree->reverse_level_begin(), data.p_tree->reverse_level_end(), [&order](const clade* c) { order.push_back(c); });

    string dir = my_input_parameters.output_prefix;
    if (dir.empty()) dir = "results";
    create_directory(dir);

    for (auto p_model : models) {

        std::vector<simulated_family> results;

        simulate_processes(p_model, results);

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst2(fname);
        print_header(ofst2, my_input_parameters, results.size());
        print_simulations(ofst2, false, results);
        LOG(INFO) << "Simulated values written to " << fname << endl;

        string truth_fname = filename("simulation_truth", dir);
        std::ofstream ofst(truth_fname);
        print_header(ofst, my_input_parameters, results.size());
        print_simulations(ofst, true, results);
        LOG(INFO) << "Simulated values (including internal nodes) written to " << truth_fname << endl;

    }
}


void simulator::print_simulations(std::ostream& ost, bool include_internal_nodes, const std::vector<simulated_family>& results) {

    std::vector<const clade *> order;
    auto fn = [&order](const clade *c) { order.push_back(c); };
    for_each(data.p_tree->reverse_level_begin(), data.p_tree->reverse_level_end(), fn);

    if (results.empty())
    {
        LOG(ERROR) << "No simulations created" << endl;
        return;

    }
    ost << "DESC\tFID";
    for (size_t i = 0; i < order.size(); ++i)
    {
        if (order[i]->is_leaf())
            ost << '\t' << order[i]->get_taxon_name();
        else if (include_internal_nodes)
            ost << '\t' << i;

    }
    ost << endl;

    for (size_t j = 0; j < results.size(); ++j) {
        auto& fam = results[j];
        // Printing gene counts
        ost << "L" << fam.lambda << "\ttranscript" << j;
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (order[i]->is_leaf() || include_internal_nodes)
            {
                ost << '\t';
                ost << fam.values.at(order[i]);
            }
        }
        ost << endl;
    }
}

root_equilibrium_distribution* create_rootdist(std::string param, const std::map<int, float>& rootdist)
{
    root_equilibrium_distribution* p_dist = nullptr;
    auto tokens = tokenize_str(param, ':');
    if (tokens.empty())
    {
        LOG(INFO) << "Using default gamma root distribution (0.75, 1.0/30.0)";
        p_dist = new root_distribution_gamma(0.75, 1.0 / 30.0);
    }
    else if (tokens[0] == "fixed")
    {
        auto fixed_value = stof(tokens[1]);
        LOG(INFO) << "Root distribution fixed at " << fixed_value;
        p_dist = new root_distribution_fixed(fixed_value);
    }
    else if (tokens[0] == "gamma")
    {
        auto alpha = stof(tokens[1]), beta = stof(tokens[2]);
        LOG(INFO) << "Using user provided gamma root distribution (" << alpha << ", " << beta << ")";
        p_dist = new root_distribution_gamma(alpha, beta);
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

    sigma lam(0.25);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    user_data data;
    data.p_tree = p_tree.get();
    data.rootdist[1] = 1;
    data.rootdist[2] = 1;
    data.rootdist[5] = 1;
    input_parameters params;
    params.rootdist_params = "file:rd.txt";
    simulator sim(data, params);

    simulated_family actual = sim.create_trial(&lam, 2);

    CHECK_EQ(doctest::Approx(5.0), actual.values.at(p_tree.get()));
    CHECK_EQ(doctest::Approx(4.85931), actual.values.at(p_tree->find_descendant("A")));
    CHECK_EQ(doctest::Approx(4.988327), actual.values.at(p_tree->find_descendant("B")));
}

TEST_CASE("binner")
{
    sigma lam(0.25);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    binner b(&lam, p_tree.get(), 5);

    CHECK_EQ(62, b.bin(2.7));
    CHECK_EQ(doctest::Approx(5.16033), b.value(120));
    CHECK_EQ(30, b.bin(1.3));
    CHECK_EQ(doctest::Approx(1.720113), b.value(40));
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("print_process_prints_in_order")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    std::ostringstream ost;
    clademap<double> t;
    t[p_tree->find_descendant("B")] = 4;
    t[p_tree->find_descendant("A")] = 2;
    t[p_tree->find_descendant("AB")] = 6;

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, true, my_trials);

    CHECK_STREAM_CONTAINS(ost, "DESC\tFID\tB\tA\t2");
    CHECK_STREAM_CONTAINS(ost, "L0\ttranscript0\t4\t2\t6");

}

TEST_CASE("print_process_can_print_without_internal_nodes")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    std::ostringstream ost;
    clademap<double> t;
    t[p_tree->find_descendant("B")] = 4;
    t[p_tree->find_descendant("A")] = 2;
    t[p_tree->find_descendant("AB")] = 6;

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, false, my_trials);
    CHECK_STREAM_CONTAINS(ost, "DESC\tFID\tB\tA\n");
    CHECK_STREAM_CONTAINS(ost, "L0\ttranscript0\t4\t2\n");

}

TEST_CASE("Check mean and variance of a simulated family leaf")
{
    randomizer_engine.seed(10);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:1):1"));
    sigma sigma(10);
    auto a = p_tree->find_descendant("A");

    matrix_cache cache;
    size_t sz = 3;  // make this larger when simulations are faster
    vector<double> v(sz);
    generate(v.begin(), v.end(), [&]() {
        auto sim = create_simulated_family(p_tree.get(), &sigma, 10);
        return sim.values[a];
        });

    auto mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
    auto variance = std::accumulate(v.begin(), v.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));
     });

    CHECK_EQ(doctest::Approx(9.37455), mean);
    CHECK_EQ(doctest::Approx(9.90501), variance);

}

TEST_CASE("print_header")
{
    std::ostringstream ost;
    input_parameters p;
    p.fixed_lambda = 2.5;
    p.rootdist_params = "Fixed:6.0";
    print_header(ost, p, 100);
    CHECK_STREAM_CONTAINS(ost, "# Simulated data set created ");
    CHECK_STREAM_CONTAINS(ost, "(100 transcripts)");
    CHECK_STREAM_CONTAINS(ost, "# Root distribution: Fixed:6.0");
    CHECK_STREAM_CONTAINS(ost, "# Sigma: 2.5");
}

TEST_CASE("print_header default rootdist")
{
    std::ostringstream ost;
    input_parameters p;
    p.fixed_lambda = 2.5;
    print_header(ost, p, 100);
    CHECK_STREAM_CONTAINS(ost, "# Root distribution: gamma:0.75:0.033");
}

TEST_CASE("print_header multiple sigmas")
{
    std::ostringstream ost;
    input_parameters p;
    p.fixed_multiple_lambdas = "1,2,3";
    print_header(ost, p, 100);
    CHECK_STREAM_CONTAINS(ost, "# Sigma: 1,2,3");
}

TEST_CASE("create_rootdist creates__specifed_distribution_if_given")
{
    randomizer_engine.seed(10);
    input_parameters params;
    root_distribution rd;
    rd[2] = 11;
    unique_ptr<root_equilibrium_distribution> red(create_rootdist("file:rd.txt", rd));
    REQUIRE(dynamic_cast<root_distribution_specific*>(red.get()));
    CHECK_EQ(2.0f, red->select_root_value(0));
}

TEST_CASE("create_rootdist creates gamma distribution if given distribution")
{
    randomizer_engine.seed(10);
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("gamma:1:3", root_distribution()));
    REQUIRE(dynamic_cast<root_distribution_gamma*>(rd.get()));
    CHECK_EQ(doctest::Approx(1.87704f), rd->select_root_value(0));
}

TEST_CASE("create_rootdist returns gamma distribution if nothing set")
{
    randomizer_engine.seed(10);
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("", root_distribution()));
    REQUIRE(dynamic_cast<root_distribution_gamma*>(rd.get()));
    CHECK_EQ(doctest::Approx(0.03538f), rd->select_root_value(0));
}

TEST_CASE("create_rootdist creates fixed root if requested")
{
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("fixed:6", root_distribution()));
    REQUIRE(dynamic_cast<root_distribution_fixed*>(rd.get()));
    CHECK_EQ(6.0, rd->select_root_value(0));

}
