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
#ifdef MODEL_GENE_EXPRESSION_LOGS
        // root_size is already in log space here, so we un-log it in order to calculate the log of the whole statement
        _max_value = log(exp(root_size) + 4.5 * sigma * sqrt(t) + LOG_OFFSET);
#else
        _max_value = root_size + 4.5 * sigma * sqrt(t);
#endif
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
    boundaries bounds(log(LOG_OFFSET), b.max_value());

    sim.values[p_tree] = root_value;

    std::function <void(const clade*)> get_child_value;
    get_child_value = [&](const clade* c) {
        double sigma = p_sigma->get_value_for_clade(c);
        MatrixXd m = ConvProp_bounds(c->get_branch_length(), sigma * sigma / 2, DiffMat::instance(), bounds);
        VectorXd v = VectorPos_bounds(sim.values[c->get_parent()], DISCRETIZATION_RANGE, bounds);
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

void print_header(std::ostream& ost, const input_parameters& p, size_t c, const clade *p_tree, const cladevector& order)
{
    auto tm = std::time(nullptr);

    ost << "# Simulated data set created " << std::put_time(std::localtime(&tm), "%Y-%m-%d %H:%M");
    ost << " (" << c << " transcripts)" << endl;
    ost << "# Root distribution: ";
    if (p.rootdist_params.empty())
        ost << "gamma:0.75:0.033" << endl;
    else
        ost << p.rootdist_params << endl;

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

        std::vector<simulated_family> results;

        simulate_processes(p_model, results);

        string fname = filename("simulation", my_input_parameters.output_prefix);
        std::ofstream ofst2(fname);
        print_header(ofst2, my_input_parameters, results.size(), data.p_tree, order);
        print_simulations(ofst2, false, results, order);
        LOG(INFO) << "Simulated values written to " << fname << endl;

        string truth_fname = filename("simulation_truth", dir);
        std::ofstream ofst(truth_fname);
        print_header(ofst, my_input_parameters, results.size(), data.p_tree, order);
        print_simulations(ofst, true, results, order);
        LOG(INFO) << "Simulated values (including internal nodes) written to " << truth_fname << endl;

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
        ost << "L" << transcript.lambda << "\ttranscript" << j;
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (!order[i]) continue;

            if (order[i]->is_leaf() || include_internal_nodes)
            {
                ost << '\t';
#ifdef MODEL_GENE_EXPRESSION_LOGS
                ost << exp(transcript.values.at(order[i]));
#else
                ost << transcript.values.at(order[i]);
#endif
            }
        }
        ost << endl;
    }
}

root_equilibrium_distribution* create_rootdist(std::string param, const vector<pair<float, int>>& rootdist)
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
    data.rootdist = vector<pair<float, int>>({ {1, 1}, { 2,1 }, { 5,1 } });

    input_parameters params;
    params.rootdist_params = "file:rd.txt";
    simulator sim(data, params);

    simulated_family actual = sim.create_trial(&lam, 2);

#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(doctest::Approx(1.7918).epsilon(0.0001), actual.values.at(p_tree.get()));
    CHECK_EQ(doctest::Approx(1.658).epsilon(0.0001), actual.values.at(p_tree->find_descendant("A")));
    CHECK_EQ(doctest::Approx(1.7765).epsilon(0.0001), actual.values.at(p_tree->find_descendant("B")));
#else
    CHECK_EQ(doctest::Approx(5.0), actual.values.at(p_tree.get()));
    CHECK_EQ(doctest::Approx(4.85931), actual.values.at(p_tree->find_descendant("A")));
    CHECK_EQ(doctest::Approx(4.988327), actual.values.at(p_tree->find_descendant("B")));
#endif
}

TEST_CASE("binner")
{
    sigma lam(0.25);
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    binner b(&lam, p_tree.get(), 5);

#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(106, b.bin(2.7));
    CHECK_EQ(doctest::Approx(3.03331), b.value(120));
    CHECK_EQ(51, b.bin(1.3));
    CHECK_EQ(doctest::Approx(1.0111), b.value(40));
#else
    CHECK_EQ(62, b.bin(2.7));
    CHECK_EQ(doctest::Approx(5.16033), b.value(120));
    CHECK_EQ(30, b.bin(1.3));
    CHECK_EQ(doctest::Approx(1.720113), b.value(40)); 
#endif
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("print_process_prints_in_order")
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    const bool use_logs = true;
#else
    const bool use_logs = false;
#endif

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    vector<pair<string, double>> values{ {"A",2}, {"B",4},{"AB",6 } };

    std::ostringstream ost;
    clademap<double> t;
    for (auto v : values)
    {
        t[p_tree->find_descendant(v.first)] = use_logs ? log(v.second) : v.second;
    }

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, true, my_trials, get_ape_order(p_tree.get()));

    CHECK_STREAM_CONTAINS(ost, "DESC\tFID\tA\tB\t3");
    CHECK_STREAM_CONTAINS(ost, "L0\ttranscript0\t2\t4\t6");

}

TEST_CASE("print_process_can_print_without_internal_nodes")
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    const bool use_logs = true;
#else
    const bool use_logs = false;
#endif

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    vector<pair<string, double>> values{ {"A",2}, {"B",4},{"AB",6 } };

    std::ostringstream ost;
    clademap<double> t;
    for (auto v : values)
    {
        t[p_tree->find_descendant(v.first)] = use_logs ? log(v.second) : v.second;
    }

    vector<simulated_family> my_trials(1);
    my_trials[0].values = t;

    user_data data;
    data.p_tree = p_tree.get();
    input_parameters params;
    simulator sim(data, params);
    sim.print_simulations(ost, false, my_trials, get_ape_order(p_tree.get()));
    CHECK_STREAM_CONTAINS(ost, "DESC\tFID\tA\tB\n");
    CHECK_STREAM_CONTAINS(ost, "L0\ttranscript0\t2\t4\n");

}

TEST_CASE("Check mean and variance of a simulated family leaf")
{
    randomizer_engine.seed(10);

    unique_ptr<clade> p_tree(parse_newick("(A:1,B:1):1"));
    sigma sigma(10);
    auto a = p_tree->find_descendant("A");

    matrix_cache cache;
    size_t sz = 3;
    vector<double> v(sz);
    generate(v.begin(), v.end(), [&]() {
        auto sim = create_simulated_family(p_tree.get(), &sigma, 10);
        return sim.values[a];
        });

    auto mean = std::accumulate(v.begin(), v.end(), 0.0) / sz;
    auto variance = std::accumulate(v.begin(), v.end(), 0.0, [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));
     });

#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(doctest::Approx(4.4569), mean);
    CHECK_EQ(doctest::Approx(2.0525), variance);
#else
    CHECK_EQ(doctest::Approx(9.37455), mean);
    CHECK_EQ(doctest::Approx(9.90501), variance);
#endif
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
    CHECK_STREAM_CONTAINS(ost, "# Root distribution: gamma:0.75:0.033");
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
#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(doctest::Approx(log(3.0f)), red->select_root_value(0));
#else
    CHECK_EQ(doctest::Approx(2.0f), red->select_root_value(0));
#endif
}

TEST_CASE("create_rootdist creates gamma distribution if given distribution")
{
    randomizer_engine.seed(10);
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("gamma:1:3", {}));
    REQUIRE(dynamic_cast<root_distribution_gamma*>(rd.get()));
#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(doctest::Approx(1.05676f), rd->select_root_value(0));
#else
    CHECK_EQ(doctest::Approx(1.87704f), rd->select_root_value(0));
#endif
}

TEST_CASE("create_rootdist returns gamma distribution if nothing set")
{
    randomizer_engine.seed(10);
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("", {}));
    REQUIRE(dynamic_cast<root_distribution_gamma*>(rd.get()));
#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(doctest::Approx(0.03477f), rd->select_root_value(0));
#else
    CHECK_EQ(doctest::Approx(0.03538f), rd->select_root_value(0));
#endif
}

TEST_CASE("create_rootdist creates fixed root if requested")
{
    user_data ud;
    unique_ptr<root_equilibrium_distribution> rd(create_rootdist("fixed:6", {}));
    REQUIRE(dynamic_cast<root_distribution_fixed*>(rd.get()));
#ifdef MODEL_GENE_EXPRESSION_LOGS
    CHECK_EQ(doctest::Approx(1.9459f), rd->select_root_value(0));
#else
    CHECK_EQ(6.0, rd->select_root_value(0));
#endif
}
