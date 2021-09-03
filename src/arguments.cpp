#include <string>
#include <fstream>

#ifdef Boost_FOUND
#include <boost/program_options.hpp>

namespace po = boost::program_options;
#elif defined(HAVE_GETOPT_H)
#include <getopt.h>
#else
#error Either Boost or GetOpt must be available
#endif

#include "doctest.h"
#include "easylogging++.h"

#include "arguments.h"

using namespace std;


#ifndef Boost_FOUND
#include <getopt.h>

struct option longopts[] = {
  { "infile", required_argument, NULL, 'i' },
  { "error_model", optional_argument, NULL, 'e' },
  { "output_prefix", required_argument, NULL, 'o'},
  { "tree", required_argument, NULL, 't' },
  { "fixed_sigma", required_argument, NULL, 'l' },
  { "fixed_multiple_sigmas", required_argument, NULL, 'm' },
  { "sigma_tree", required_argument, NULL, 'y' },
  { "n_gamma_cats", required_argument, NULL, 'k' },
  { "fixed_alpha", required_argument, NULL, 'a' },
  { "simulate", optional_argument, NULL, 's' },
  { "pvalue", required_argument, NULL, 'P' },
  { "zero_root", no_argument, NULL, 'z' },
  { "cores", required_argument, NULL, 'c' },
  { "lambda_per_family", no_argument, NULL, 'b' },
  { "log_config", required_argument, NULL, 'L' },
  { "optimizer_expansion", optional_argument, NULL, 'E' },
  { "optimizer_reflection", optional_argument, NULL, 'R' },
  { "optimizer_iterations", optional_argument, NULL, 'I' },
  { "help", no_argument, NULL, 'h'},
  { 0, 0, 0, 0 }
};
#else

template<typename T>
void maybe_set(const po::variables_map& vm, string key, T& value)
{
    if (vm.find(key) != vm.end())
        value = vm[key].as<T>();
}
#endif

#ifndef Boost_FOUND
input_parameters read_arguments(int argc, char* const argv[])
{
    input_parameters my_input_parameters;
    if (argc == 1)
    {
        my_input_parameters.help = true;
        return my_input_parameters;
    }

    int args; // getopt_long returns int or char
    int prev_arg;

    while (prev_arg = optind, (args = getopt_long(argc, argv, "c:v:i:e::o:t:y:n:f:E:F:R:L:P:I:l:m:k:a:g:s::p::zbh", longopts, NULL)) != -1) {
        if (optind == prev_arg + 2 && optarg && *optarg == '-') {
            LOG(ERROR) << "You specified option " << argv[prev_arg] << " but it requires an argument. Exiting..." << endl;
            exit(EXIT_FAILURE);
        }

        switch (args) {
        case 'a':
            my_input_parameters.fixed_alpha = atof(optarg);
            break;
        case 'b':
            my_input_parameters.lambda_per_family = true;
            break;
        case 'c':
            my_input_parameters.cores = atoi(optarg);
            break;
        case 'e':
            my_input_parameters.use_error_model = true;
            if (optarg)
                my_input_parameters.error_model_file_path = optarg;
            break;
        case 'f':
            my_input_parameters.rootdist = optarg;
            break;
        case 'h':
            my_input_parameters.help = true;
            break;
        case 'i':
            my_input_parameters.input_file_path = optarg;
            break;
        case 'k':
            if (optarg != NULL) { my_input_parameters.n_gamma_cats = atoi(optarg); }
            break;
        case 'l':
            my_input_parameters.fixed_lambda = atof(optarg);
            break;
        case 'm':
            my_input_parameters.fixed_multiple_lambdas = optarg;
            break;
        case 'o':
            my_input_parameters.output_prefix = optarg;
            break;
        case 'p':
            my_input_parameters.use_poisson_dist_for_prior = true; // If the user types '-p', the root eq freq dist will not be a uniform
                                                             // If the user provides an argument to -p, then we do not estimate it
            if (optarg != NULL) { my_input_parameters.poisson_lambda = atof(optarg); }
            break;
        case 's':
            // Number of fams simulated defaults to 0 if -f is not provided
            my_input_parameters.is_simulating = true;
            if (optarg != NULL) { my_input_parameters.nsims = atoi(optarg); }
            break;
        case 't':
            my_input_parameters.tree_file_path = optarg;
            break;
        case 'v':
            my_input_parameters.verbose_logging_level = atoi(optarg);
            break;
        case 'y':
            my_input_parameters.lambda_tree_file_path = optarg;
            break;
        case 'z':
            my_input_parameters.exclude_zero_root_families = true;
            break;
        case 'E':
            my_input_parameters.optimizer_params.neldermead_expansion = atof(optarg);
            break;
        case 'F':
            my_input_parameters.fixed_root_value = atof(optarg);
            break;
        case 'I':
            my_input_parameters.optimizer_params.neldermead_iterations = atoi(optarg);
            break;
        case 'L':
            my_input_parameters.log_config_file = optarg;
            break;
        case 'P':
            my_input_parameters.pvalue = atof(optarg);
            break;
        case 'R':
            my_input_parameters.optimizer_params.neldermead_reflection = atof(optarg);
            break;
        case ':':   // missing argument
            fprintf(stderr, "%s: option `-%c' requires an argument",
                argv[0], optopt);
            break;
        default: // '?' is parsed (
            throw std::runtime_error(string("Unrecognized parameter: '") + (char)args + "'");

        }
    }

    if (optind < argc)
    {
        throw std::runtime_error(string("Unrecognized parameter: '") + argv[optind] + "'");
    }

    my_input_parameters.check_input(); // seeing if options are not mutually exclusive              

    return my_input_parameters;
}
#else   // Boost version 
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
            store(po::parse_config_file(ifs, config_file_options), vm);
            notify(vm);
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

    if (vm.find("rootdist") != vm.end())
        my_input_parameters.rootdist_params = rootdist_options(vm["rootdist"].as<string>());

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
#endif

void input_parameters::check_input() {
    //! Options -l and -m cannot both specified.
    if (fixed_lambda > 0.0 && !fixed_multiple_lambdas.empty()) {
        throw runtime_error("Options -l and -m are mutually exclusive.");
    }

    //! Option -m requires a lambda tree (-y)
    if (!fixed_multiple_lambdas.empty() && lambda_tree_file_path.empty()) {
        throw runtime_error("Multiple lambda values (-m) specified with no lambda tree (-y)");
    }

    //! Options -l and -i have to be both specified (if estimating and not simulating).
    if (fixed_lambda > 0.0 && input_file_path.empty() && !is_simulating) {
        throw runtime_error("Options -l and -i must both be provided an argument.");
    }

    if (is_simulating)
    {
        // Must specify a lambda
        if (fixed_lambda <= 0.0 && fixed_multiple_lambdas.empty()) {
            throw runtime_error("Cannot simulate without initial sigma values");
        }

        if (fixed_alpha <= 0.0 && this->n_gamma_cats > 1) {
            throw runtime_error("Cannot simulate gamma clusters without an alpha value");
        }

        //! Options -i and -f cannot be both specified. Either one or the other is used to specify the root eq freq distr'n.
        if (!input_file_path.empty() && rootdist_params.type == rootdist_options::file) {
            throw runtime_error("Options -i and -f are mutually exclusive.");
        }
    }
    else
    {
        if (fixed_alpha >= 0.0 && n_gamma_cats == 1) {
            throw runtime_error("Alpha specified with 1 gamma category.");
        }


        if (lambda_per_family)
        {
            if (input_file_path.empty())
                throw runtime_error("No family file provided");
            if (tree_file_path.empty())
                throw runtime_error("No tree file provided");
        }

        if (n_gamma_cats > 1 && use_error_model && error_model_file_path.empty())
        {
            throw runtime_error("Estimating an error model with a gamma distribution is not supported at this time");
        }

    }
}

input_parameters read_arguments(int argc, char* const argv[]);

struct option_test
{
    char* values[100];
    size_t argc;

    option_test(vector<string> arguments)
    {
#ifdef HAVE_GETOPT_H
        optind = 0;
#endif
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
    option_test c({ "cafe5", "--rootdist=fixed:12.7"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(doctest::Approx(12.7), actual.rootdist_params.fixed_value);
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

TEST_CASE("Options: check_input_does_not_throw_when_simulating_with_multiple_lambdas")
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
    CHECK_THROWS_WITH_AS(params.check_input(), "No family file provided", runtime_error);
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

TEST_CASE("Options: Cannot specify rootdist for simulations with rootdist file and transcript file")
{
    input_parameters params;
    params.fixed_lambda = 10;
    params.input_file_path = "transcripts.txt";
    params.rootdist_params.type = rootdist_options::file;
    params.check_input();
    CHECK(true);

    params.is_simulating = true;
    CHECK_THROWS_WITH_AS(params.check_input(), "Options -i and -f are mutually exclusive.", runtime_error);

}

