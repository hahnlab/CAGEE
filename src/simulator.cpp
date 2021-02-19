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

using namespace std;

extern std::mt19937 randomizer_engine; // seeding random number engine

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

// At the root, we have a vector of length DISCRETIZATION_RANGE. This has probability 1 at the size of the root
// and 0 everywhere else
// for each child, generate the transition matrix and multiply
// The tip values would then be taken to be the highest probability entries in the vector for that node#endif

simulated_family simulator::create_trial(const lambda *p_lambda, int family_number, const matrix_cache& cache) {

    if (data.p_tree == NULL)
        throw runtime_error("No tree specified for simulation");

    simulated_family result;
    result.lambda = get_lambda_values(p_lambda)[0];

    clademap<vector<double>> probs;
    data.p_tree->apply_prefix_order([&probs](const clade*c) { probs[c].resize(DISCRETIZATION_RANGE); });

    probs.at(data.p_tree)[data.prior.select_root_size(family_number)] = 1;

    std::function <void(const clade*)> get_child_probability_vector;
    get_child_probability_vector = [&](const clade* c) {
        p_lambda->calculate_child_factor(cache, c, probs.at(c->get_parent()), 0, DISCRETIZATION_RANGE - 1, 0, DISCRETIZATION_RANGE - 1, probs.at(c).data());
        c->apply_to_descendants(get_child_probability_vector);
    };

    data.p_tree->apply_to_descendants(get_child_probability_vector);

    for (auto it = data.p_tree->reverse_level_begin(); it != data.p_tree->reverse_level_end(); ++it)
    {
        auto max = max_element(probs[*it].begin(), probs[*it].end());
        result.values[*it] = max - probs[*it].begin();
    }
    return result;
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
        unique_ptr<lambda> sim_lambda(p_model->get_simulation_lambda());
        
        matrix_cache cache;
        //cache.precalculate_matrices(get_lambda_values(sim_lambda.get()), this->data.p_tree->get_branch_lengths());
        p_model->prepare_matrices_for_simulation(cache);

        int n = 0;

        auto end_it = i + LAMBDA_PERTURBATION_STEP_SIZE > results.size() ? results.end() : results.begin() + i + LAMBDA_PERTURBATION_STEP_SIZE;
        generate(results.begin()+i, end_it, [this, &sim_lambda, i, &cache, &n]() mutable {
            return create_trial(sim_lambda.get(), i+n++, cache);
        });
    }
}

extern void write_average_multiplier(std::ostream& ost);

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

    single_lambda lam(0.25);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));

    user_data data;
    data.p_tree = p_tree.get();
    data.rootdist[1] = 1;
    data.rootdist[2] = 1;
    data.rootdist[5] = 1;
    data.max_family_size = 10;
    data.max_root_family_size = 10;

    data.prior = root_equilibrium_distribution(data.rootdist);
    input_parameters params;
    simulator sim(data, params);

    matrix_cache cache;
    cache.precalculate_matrices(get_lambda_values(&lam), { 1,3,7 });

    simulated_family actual = sim.create_trial(&lam, 2, cache);

    CHECK_EQ(5, actual.values.at(p_tree.get()));
    CHECK_EQ(5, actual.values.at(p_tree->find_descendant("A")));
    CHECK_EQ(5, actual.values.at(p_tree->find_descendant("B")));
}

