#include <string>
#include <fstream>
#include <iterator>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "doctest.h"
#include "easylogging++.h"

#include "arguments.h"
#include "io.h"

using namespace std;

template<typename T>
void maybe_set(const po::variables_map& vm, string key, T& value)
{
    if (vm.find(key) != vm.end())
        value = vm[key].as<T>();
}

input_parameters read_arguments(int argc, char* const argv[])
{
    input_parameters my_input_parameters;
    if (argc == 1)
    {
        my_input_parameters.help = true;
        return my_input_parameters;
    }

    string config_file;
    po::options_description generic("Generic options");
    generic.add_options()
        ("help", "produce help message")
        ("version,v", "print version string")
        ("config,c", po::value<string>(&config_file),
            "Configuration file containing additional options");

    po::options_description config("Configuration options (May be specified on command line or in file)");
    config.add_options()
        ("infile,i", po::value<string>(), "input file")
        ("tree,t", po::value<string>(), "Tree file in Newick format")
        ("output_prefix,o", po::value<string>(), "Output prefix")
        ("fixed_multiple_sigmas,m", po::value<string>(), "Output prefix")
        ("pvalue,P", po::value<double>(), "PValue")
        ("cores", po::value<int>(), "Number of processing cores to use, requires an integer argument. Default=All available cores.")
        ("optimizer_expansion,E", po::value<double>())
        ("optimizer_reflection,R", po::value<double>())
        ("optimizer_iterations,I", po::value<int>())
        ("error,e", po::value<string>()->default_value("false")->implicit_value("true"))
        ("rootdist", po::value<string>())
        ("prior", po::value<string>())
        ("verbose", po::value<int>())
        ("sample_group", po::value<vector<string>>())
        ("zero_root,z", po::value<bool>()->implicit_value(true))
        ("sigma_tree,y", po::value<string>(), "Path to sigma tree, for use with multiple sigmas")
        ("simulate,s", po::value<string>()->default_value("false"), "Simulate families. Optionally provide the number of simulations to generate")
        ("fixed_sigma,l", po::value<double>(), "Value for a single user provided sigma value, otherwise sigma is estimated.")
        ;
    
    po::options_description config_file_options;
    config_file_options.add(config);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
        options(cmdline_options).run(), vm);
    notify(vm);

    if (vm.count("help")) {
        my_input_parameters.help = true;
    }

    if (vm.count("version")) {
        my_input_parameters.help = true;
    }

    if (!config_file.empty())
    {
        ifstream ifs(config_file.c_str());
        if (!ifs)
        {
            throw std::runtime_error("Config file not found: " + config_file);
        }
        else
        {
            try
            {
                store(po::parse_config_file(ifs, config_file_options), vm);
                notify(vm);
            }
            catch (boost::program_options::unknown_option& e)
            {
                throw std::runtime_error("Unknown option: " + e.get_option_name());
            }
        }
    }

    //maybe_set(vm, "fixed_root_value", my_input_parameters.fixed_root_value);
    maybe_set(vm, "fixed_sigma", my_input_parameters.fixed_lambda);
    maybe_set(vm, "pvalue", my_input_parameters.pvalue);
    maybe_set(vm, "infile", my_input_parameters.input_file_path);
    maybe_set(vm, "output_prefix", my_input_parameters.output_prefix);
    maybe_set(vm, "tree", my_input_parameters.tree_file_path);
    maybe_set(vm, "zero_root", my_input_parameters.exclude_zero_root_families);
    maybe_set(vm, "cores", my_input_parameters.cores);
    maybe_set(vm, "optimizer_expansion", my_input_parameters.optimizer_params.neldermead_expansion);
    maybe_set(vm, "optimizer_reflection", my_input_parameters.optimizer_params.neldermead_reflection);
    maybe_set(vm, "optimizer_iterations", my_input_parameters.optimizer_params.neldermead_iterations);
    maybe_set(vm, "fixed_multiple_sigmas", my_input_parameters.fixed_multiple_lambdas);
    maybe_set(vm, "sigma_tree", my_input_parameters.lambda_tree_file_path);
    maybe_set(vm, "verbose", my_input_parameters.verbose_logging_level);
    maybe_set(vm, "rootdist", my_input_parameters.rootdist_params);
    maybe_set(vm, "sample_group", my_input_parameters.sample_groups);

    if (vm.find("prior") != vm.end())
    {
        auto tokens = tokenize_str(vm["prior"].as<string>(), ':');
        auto alpha = stof(tokens[1]), beta = stof(tokens[2]);
        my_input_parameters.prior = gamma_distribution<double>(alpha, beta);
    }
    string simulate_string = vm["simulate"].as<string>();
    my_input_parameters.is_simulating = simulate_string != "false";
    if (my_input_parameters.is_simulating)
    {
        if (simulate_string.find_first_not_of("0123456789") == std::string::npos)
            my_input_parameters.nsims = stoi(simulate_string);
        else
            my_input_parameters.nsims = 0;
    }

    string error = vm["error"].as<string>();
    my_input_parameters.use_error_model = error != "false";
    if (my_input_parameters.use_error_model && error != "true")
    {
        my_input_parameters.error_model_file_path = error;
    }

    my_input_parameters.check_input(); // seeing if options are not mutually exclusive              

    return my_input_parameters;
}

void input_parameters::check_input() {
    vector<string> errors;

    //! Options -l and -m cannot both specified.
    if (fixed_lambda > 0.0 && !fixed_multiple_lambdas.empty()) {
        errors.push_back("Options -l and -m are mutually exclusive.");
    }

    //! Option -m requires a lambda tree (-y)
    if (!fixed_multiple_lambdas.empty() && lambda_tree_file_path.empty()) {
        errors.push_back("Multiple lambda values (-m) specified with no lambda tree (-y)");
    }

    //! Options -l and -i have to be both specified (if estimating and not simulating).
    if (fixed_lambda > 0.0 && input_file_path.empty() && !is_simulating) {
        errors.push_back("Options -l and -i must both be provided an argument.");
    }

    if (is_simulating)
    {
        // Must specify a lambda
        if (fixed_lambda <= 0.0 && fixed_multiple_lambdas.empty()) {
            errors.push_back("Cannot simulate without initial sigma values");
        }

        if (fixed_alpha <= 0.0 && this->n_gamma_cats > 1) {
            errors.push_back("Cannot simulate gamma clusters without an alpha value");
        }

        //! Options -i and -f cannot be both specified. Either one or the other is used to specify the root eq freq distr'n.
        if (!input_file_path.empty()) {
            errors.push_back("A families file was provided while simulating");
        }
    }
    else
    {
        if (fixed_alpha >= 0.0 && n_gamma_cats == 1) {
            errors.push_back("Alpha specified with 1 gamma category.");
        }

        if (!rootdist_params.empty()) {
            errors.push_back("A root distribution was provided while estimating");
        }

        if (lambda_per_family)
        {
            if (input_file_path.empty())
                errors.push_back("No family file provided");
            if (tree_file_path.empty())
                errors.push_back("No tree file provided");
        }

        if (n_gamma_cats > 1 && use_error_model && error_model_file_path.empty())
        {
            errors.push_back("Estimating an error model with a gamma distribution is not supported at this time");
        }
    }

    if (errors.size() == 1)
    {
        throw std::runtime_error(errors[0]);
    }
    else if (!errors.empty())
    {
        ostringstream ost;
        ostream_iterator<string> lm(ost, "\n");
        copy(errors.begin(), errors.end(), lm);
        throw std::runtime_error(ost.str());
    }
}

input_parameters read_arguments(int argc, char* const argv[]);

struct option_test
{
    char* values[100];
    size_t argc;

    option_test(vector<string> arguments)
    {
        argc = arguments.size();
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            values[i] = strdup(arguments[i].c_str());
        }
    }

    ~option_test()
    {
        for (size_t i = 0; i < argc; ++i)
        {
            free(values[i]);
        }
    }
};

TEST_CASE("read_arguments translates short values ") {
    option_test c({ "cafe5", "-ifile" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}

TEST_CASE("read_arguments translates long values ") {
    option_test c({ "cafe5", "--infile", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}

TEST_CASE("Options, input_short_space_separated")
{
    option_test c({ "cafe5", "-i", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}
TEST_CASE("Options, simulate_long")
{
    option_test c({ "cafe5", "--simulate=1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, simulate_short")
{
    option_test c({ "cafe5", "-s1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, pvalue_long")
{
    option_test c({ "cafe5", "--pvalue=0.01" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.01, actual.pvalue);
}

TEST_CASE("Options, pvalue_short")
{
    option_test c({ "cafe5", "-P0.01" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.01, actual.pvalue);
}

TEST_CASE("Options, optimizer_long")
{
    option_test c({ "cafe5", "--optimizer_expansion=0.05", "--optimizer_reflection=3.2", "--optimizer_iterations=5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.neldermead_iterations);
}

TEST_CASE("Options, optimizer_short")
{
    option_test c({ "cafe5", "-E", "0.05", "-R", "3.2", "-I", "5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.neldermead_iterations);
}

TEST_CASE("Options, cores_long")
{
    option_test c({ "cafe5", "--cores=6" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(6, actual.cores);
}

TEST_CASE("Options, multiple_sigmas_long")
{
    option_test c({ "cafe5", "--fixed_multiple_sigmas=5,10,15", "--sigma_tree=foo.txt" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("5,10,15", actual.fixed_multiple_lambdas);
}

TEST_CASE("Options: errormodel_accepts_argument")
{
    option_test c({ "cafe5", "-eerror.txt" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path == "error.txt");
}

TEST_CASE("Options: fixed_root_value")
{
    option_test c({ "cafe5", "--rootdist=fixed:12.7", "--simulate=100", "--fixed_sigma=0.01"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("fixed:12.7", actual.rootdist_params);
}

TEST_CASE("Options: sample_group")
{
    option_test c({ "cafe5", "--sample_group=heart,lungs", "--sample_group=brain" });

    auto actual = read_arguments(c.argc, c.values);

    REQUIRE_EQ(2, actual.sample_groups.size());
    CHECK_EQ("heart,lungs", actual.sample_groups[0]);
    CHECK_EQ("brain", actual.sample_groups[1]);
}

TEST_CASE("Options, errormodel_accepts_no_argument")
{
    option_test c({ "cafe5", "-e" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path.empty());
}

TEST_CASE("Options: zero_root_familes")
{
    input_parameters by_default;
    CHECK_FALSE(by_default.exclude_zero_root_families);

    option_test c({ "cafe5", "-z" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.exclude_zero_root_families);
}

TEST_CASE("Options: must_specify_sigma_for_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate without initial sigma values");
}

TEST_CASE("Options: must_specify_lambda_and_input_file_for_estimator")
{
    input_parameters params;
    params.fixed_lambda = 0.05;
    CHECK_THROWS_WITH(params.check_input(), "Options -l and -i must both be provided an argument.");
}

TEST_CASE("Options: must_specify_alpha_for_gamma_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_lambda = 0.05;
    params.n_gamma_cats = 3;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate gamma clusters without an alpha value");
}

TEST_CASE("Options: must_specify_alpha_and_k_for_gamma_inference")
{
    input_parameters params;
    params.fixed_alpha = 0.7;
    CHECK_THROWS_WITH(params.check_input(), "Alpha specified with 1 gamma category.");
}

TEST_CASE("Options: can_specify_alpha_without_k_for_gamma_simulation")
{
    input_parameters params;
    params.fixed_alpha = 0.7;
    params.fixed_lambda = 0.01;
    params.is_simulating = true;
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: check_input_does_not_throw_when_simulating_with_multiple_sigmas")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_multiple_lambdas = "0.01,0.05";
    params.lambda_tree_file_path = "./tree";
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: per_family_must_provide_families")
{
    input_parameters params;
    params.lambda_per_family = true;
    CHECK_THROWS_WITH_AS(params.check_input(), "No family file provided\nNo tree file provided\n", runtime_error);
}

TEST_CASE("Options: per_family_must_provide_tree")
{
    input_parameters params;
    params.lambda_per_family = true;
    params.input_file_path = "/tmp/test";
    CHECK_THROWS_WITH_AS(params.check_input(), "No tree file provided", runtime_error);
}

TEST_CASE("Options: cannot_estimate_error_and_gamma_together")
{
    input_parameters params;
    params.n_gamma_cats = 3;
    params.use_error_model = true;
    params.error_model_file_path = "model.txt";
    params.check_input();
    CHECK(true);

    params.error_model_file_path.clear();
    CHECK_THROWS_WITH_AS(params.check_input(), "Estimating an error model with a gamma distribution is not supported at this time", runtime_error);

}

TEST_CASE("Specifying a --rootdist argument without the --simulate argument will result in an error")
{
    input_parameters params;
    params.fixed_lambda = 10;
    params.rootdist_params = "file:rd.txt";
    params.is_simulating = true;
    params.check_input();
    CHECK(true);

    params.is_simulating = false;
    params.input_file_path = "transcripts.txt";
    CHECK_THROWS_WITH_AS(params.check_input(), "A root distribution was provided while estimating", runtime_error);

}

TEST_CASE("Cannot specify an input file when simulating")
{
    input_parameters params;
    params.fixed_lambda = 10;
    params.is_simulating = true;
    params.check_input();
    CHECK(true);

    params.input_file_path = "transcripts.txt";
    CHECK_THROWS_WITH_AS(params.check_input(), "A families file was provided while simulating", runtime_error);

}
