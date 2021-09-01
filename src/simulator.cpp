#include <numeric>
#include <algorithm>
#include <fstream>
#include <random>
#include <vector>

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
    double root_size = data.p_prior->select_root_value(family_number);

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
        print_simulations(ofst2, false, results);
        LOG(INFO) << "Simulated values written to " << fname << endl;

        string truth_fname = filename("simulation_truth", dir);
        std::ofstream ofst(truth_fname);
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
        ost << "L" << fam.lambda << "\tsimfam" << j;
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
    data.max_family_size = 10;
    data.max_root_family_size = 10;

    data.p_prior = new root_distribution_specific(data.rootdist);
    input_parameters params;
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

#define STRCMP_CONTAINS(x, y) CHECK(strstr(y,x) != nullptr)

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

    STRCMP_CONTAINS("DESC\tFID\tB\tA\t2", ost.str().c_str());
    STRCMP_CONTAINS("L0\tsimfam0\t4\t2\t6", ost.str().c_str());

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
    STRCMP_CONTAINS("DESC\tFID\tB\tA\n", ost.str().c_str());
    STRCMP_CONTAINS("L0\tsimfam0\t4\t2\n", ost.str().c_str());

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
