#include <string>
#include <fstream>
#include <iterator>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "doctest.h"
#include "easylogging++.h"

#include "arguments.h"
#include "io.h"

#ifdef INCLUDE_GIT_IN_VERSION_INFO
#include "../git_version.h"
#endif

using namespace std;

template<typename T>
void maybe_set(const po::variables_map& vm, string key, T& value)
{
    if (vm.find(key) != vm.end())
        value = vm[key].as<T>();
}

void show_version()
{
    string desc = "\nUsage: cagee [options]\n\n"
        "CAGEE is a software that provides a statistical foundation for evolutionary inferences about changes in gene expression.\n "
        "The program employs a Brownian motion process to model gene expression across a user-specified phylogenetic tree,\n "
        "thus accounting for the species phylogenetic history. The distribution of gene expression generated under this model can\n "
        "provide a basis for assessing the significance of the observed expression differences among taxa.\n\n";

    cout << PROJECT_NAME " " PROJECT_VER << endl;
    cout << desc;

#ifdef INCLUDE_GIT_IN_VERSION_INFO
    if (GitMetadata::Populated()) {
        cout << "Last git commit: \n";
        if (GitMetadata::AnyUncommittedChanges()) {
            std::cerr << "WARN: there were uncommitted changes at build-time." << std::endl;
        }
        std::cout << "commit " << GitMetadata::CommitSHA1() << " (HEAD)\n"
            << "describe " << GitMetadata::Describe() << "\n"
            << "Author: " << GitMetadata::AuthorName() << " <" << GitMetadata::AuthorEmail() << ">\n"
            << "Date: " << GitMetadata::CommitDate() << "\n\n"
            << GitMetadata::CommitSubject() << "\n" << GitMetadata::CommitBody() << std::endl;
    }
#endif
}
void show_help(const po::options_description& gen, const po::options_description& required, const po::options_description& common, const po::options_description& rare)
{
    show_version();
    std::cout << gen << endl;
    std::cout << required << endl;
    std::cout << common << endl;
    std::cout << rare << endl;
}


input_parameters read_arguments(int argc, char* const argv[])
{
    input_parameters my_input_parameters;
    my_input_parameters.command_line = std::accumulate(argv, argv + argc, std::string(), [](std::string x, std::string y) { return x + y + " "; });

    if (argc == 1)
    {
        show_version();
        my_input_parameters.help = true;
        return my_input_parameters;
    }

    string config_file;
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce help message")
        ("version,v", "print version string")
        ("config,c", po::value<string>(&config_file),
            "Configuration file containing additional options");

    po::options_description required("Required options (May be specified on command line or in file)");
    required.add_options()
        ("infile,i", po::value<string>(), "Path to tab delimited gene families file to be analyzed - Required for estimation")
        ("tree,t", po::value<string>(), "Path to file containing newick formatted tree - Required for estimation.");
        
    po::options_description common("Configuration options (May be specified on command line or in file)");
    common.add_options()
        ("cores,@", po::value<int>(), "Number of processing cores to use, requires an integer argument. Default=All available cores.")
        ("n_gamma_cats,@", po::value<int>(), "Number of gamma categories to use, requires an integer argument. Default=1 (No gamma modelling)")
        ("error,e", po::value<string>()->default_value("false")->implicit_value("true"), "Run with no file name to estimate the global error model file. This file can be provided"
            "in subsequent runs by providing the path to the Error model file with no spaces(e.g. - eBase_error_model.txt).")
        ("output_prefix,o", po::value<string>(), " Output directory - Name of directory automatically created for output. Default=results.")
        ("fixed_multiple_sigmas,m", po::value<string>(), "Multiple sigma values, comma separated.")
        ("rootdist", po::value<string>(), "Distribution of the root in Simulation Mode (mutually exclusive with --prior in Inference Mode). Can be gamma:[k]:[theta], fixed:[count], or a path/to/tab_sep_file.txt with two columns: trascripts names and their counts. Default=gamma:0.375:1600.0")
        ("prior", po::value<string>(), "Expected distribution of the root in Inference Mode (mutually exclusive with --rootdist in Simulation Mode). Must be gamma:[k]:[theta].  Default=gamma:0.375:1600.0")
        ("verbose", po::value<int>())
        ("replicate_map", po::value<string>(), "Filename of a file containing a list of specie replicates to be combined into a single species")
        ("sample_group", po::value<vector<string>>(), "Specifies sample groups (if any) for which to infer sigma^2.  Each sample and sigma^2 estimate requires a --sample_group [your_sample_A] arg, or combine them with comma: --sample_group [your_sample_A,your_sample_B,...].  Optional, no default.")
        ("sigma_tree,y", po::value<string>(), "Path to sigma tree, for use with multiple sigmas")
        ("simulate,s", po::value<string>()->default_value("false"), "Simulate families. Optionally provide the number of simulations to generate")
        ("fixed_sigma,l", po::value<double>(), "Value for a single user provided sigma value, otherwise sigma is estimated.")
        ("fixed_alpha", po::value<double>(), "Value for a single user provided alpha value, otherwise alpha is estimated.")
        ("count_all_changes", po::value<bool>()->implicit_value(true), "Reconstruction will count all changes rather than only credible changes")
        ("ratio", po::value<bool>()->implicit_value(true), "The input file contains ratios of gene expression values rather than absolute values")
        ;
    
    po::options_description rare("Less Common Options");
    rare.add_options()
        ("discretization_size,D", po::value<int>()->default_value(200), "Size (length) of the discretization vector, Default=200. Can increase resolution at the cost of computation time.")
        ("zero_root,z", po::value<bool>()->implicit_value(true), "Exclude gene families that don't exist at the root, not recommended.")
        ("free_rate", po::value<string>()->default_value("")->implicit_value("global"), "Calculate values using Free Rate Model")
        ("optimizer_expansion,E", po::value<double>(), "Expansion parameter for Nelder-Mead optimizer, Default=2.")
        ("optimizer_reflection,R", po::value<double>(), "Reflection parameter for Nelder-Mead optimizer, Default=1.")
        ("optimizer_iterations,I", po::value<int>(), "Maximum number of iterations that will be performed in "
            "sigma search.Default = 300 (increase this number if likelihood is still improving when limit is hit).");

    po::options_description config_file_options;
    config_file_options.add(required).add(common).add(rare);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(common).add(rare);

    po::variables_map vm;
    try
    {
        store(po::command_line_parser(argc, argv).
            options(cmdline_options).run(), vm);
        notify(vm);

        if (vm.count("help")) {
            show_help(generic, required, common, rare);
            my_input_parameters.help = true;
        }

        if (vm.count("version")) {
            show_version();
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
    }
    catch (boost::program_options::unknown_option& e)
    {
        throw std::runtime_error(e.what());
    }

    //maybe_set(vm, "fixed_root_value", my_input_parameters.fixed_root_value);
    maybe_set(vm, "fixed_sigma", my_input_parameters.fixed_sigma);
    maybe_set(vm, "infile", my_input_parameters.input_file_path);
    maybe_set(vm, "output_prefix", my_input_parameters.output_prefix);
    maybe_set(vm, "tree", my_input_parameters.tree_file_path);
    maybe_set(vm, "zero_root", my_input_parameters.exclude_zero_root_transcripts);
    maybe_set(vm, "cores", my_input_parameters.cores);
    maybe_set(vm, "optimizer_expansion", my_input_parameters.optimizer_params.neldermead_expansion);
    maybe_set(vm, "optimizer_reflection", my_input_parameters.optimizer_params.neldermead_reflection);
    maybe_set(vm, "optimizer_iterations", my_input_parameters.optimizer_params.max_iterations);
    maybe_set(vm, "fixed_multiple_sigmas", my_input_parameters.fixed_multiple_sigmas);
    maybe_set(vm, "sigma_tree", my_input_parameters.sigma_tree_file_path);
    maybe_set(vm, "verbose", my_input_parameters.verbose_logging_level);
    maybe_set(vm, "rootdist", my_input_parameters.rootdist_params);
    maybe_set(vm, "sample_group", my_input_parameters.sample_groups);
    maybe_set(vm, "prior", my_input_parameters.prior);
    maybe_set(vm, "discretization_size", my_input_parameters.discretization_size);
    maybe_set(vm, "count_all_changes", my_input_parameters.count_all_changes);
    maybe_set(vm, "ratio", my_input_parameters.input_file_has_ratios);
    maybe_set(vm, "replicate_map", my_input_parameters.replicate_model_file_path);
    maybe_set(vm, "n_gamma_cats", my_input_parameters.n_gamma_cats);
    maybe_set(vm, "fixed_alpha", my_input_parameters.fixed_alpha);
    maybe_set(vm, "free_rate", my_input_parameters.free_rate);

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

    if (!my_input_parameters.free_rate.empty() && vm.find("optimizer_iterations") == vm.end()) {
        my_input_parameters.optimizer_params.max_iterations = 2;
    }

    my_input_parameters.check_input(); // seeing if options are not mutually exclusive              

    return my_input_parameters;
}

void input_parameters::check_input() {
    vector<string> errors;

    //! Options -l and -m cannot both specified.
    if (fixed_sigma > 0.0 && !fixed_multiple_sigmas.empty()) {
        errors.push_back("Options -l and -m are mutually exclusive.");
    }

    //! Option -m requires a sigma tree (-y)
    if (!fixed_multiple_sigmas.empty() && sigma_tree_file_path.empty()) {
        errors.push_back("Multiple sigma values (-m) specified with no sigma tree (-y)");
    }

    //! Options -l and -i have to be both specified (if estimating and not simulating).
    if (fixed_sigma > 0.0 && input_file_path.empty() && !is_simulating) {
        errors.push_back("Options -l and -i must both be provided an argument.");
    }

    if (is_simulating)
    {
        // Must specify a sigma
        if (fixed_sigma <= 0.0 && fixed_multiple_sigmas.empty()) {
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

TEST_CASE("read_arguments stores the command line") {
    option_test c({ "cagee", "--infile", "foo.txt", "-R", "7"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("cagee --infile foo.txt -R 7 ", actual.command_line);
}

TEST_CASE("read_arguments translates short values ") {
    option_test c({ "cagee", "-ifile" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}

TEST_CASE("read_arguments translates long values ") {
    option_test c({ "cagee", "--infile", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}

TEST_CASE("Options, input_short_space_separated")
{
    option_test c({ "cagee", "-i", "file" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.input_file_path.compare("file") == 0);
}
TEST_CASE("Options, simulate_long")
{
    option_test c({ "cagee", "--simulate=1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, simulate_short")
{
    option_test c({ "cagee", "-s1000", "-l", "0.05" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(1000, actual.nsims);
}

TEST_CASE("Options, optimizer_long")
{
    option_test c({ "cagee", "--optimizer_expansion=0.05", "--optimizer_reflection=3.2", "--optimizer_iterations=5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.max_iterations);
}

TEST_CASE("Options, optimizer_short")
{
    option_test c({ "cagee", "-E", "0.05", "-R", "3.2", "-I", "5" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(0.05, actual.optimizer_params.neldermead_expansion);
    CHECK_EQ(3.2, actual.optimizer_params.neldermead_reflection);
    CHECK_EQ(5, actual.optimizer_params.max_iterations);
}

TEST_CASE("Options, cores_long")
{
    option_test c({ "cagee", "--cores=6" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(6, actual.cores);
}

TEST_CASE("Options, multiple_sigmas_long")
{
    option_test c({ "cagee", "--fixed_multiple_sigmas=5,10,15", "--sigma_tree=foo.txt" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("5,10,15", actual.fixed_multiple_sigmas);
}

TEST_CASE("Options: errormodel_accepts_argument")
{
    option_test c({ "cagee", "-eerror.txt" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path == "error.txt");
}

TEST_CASE("Options: fixed_root_value")
{
    option_test c({ "cagee", "--rootdist=fixed:12.7", "--simulate=100", "--fixed_sigma=0.01"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("fixed:12.7", actual.rootdist_params);
}

TEST_CASE("Options: sample_group")
{
    option_test c({ "cagee", "--sample_group=heart,lungs", "--sample_group=brain" });

    auto actual = read_arguments(c.argc, c.values);

    REQUIRE_EQ(2, actual.sample_groups.size());
    CHECK_EQ("heart,lungs", actual.sample_groups[0]);
    CHECK_EQ("brain", actual.sample_groups[1]);
}

TEST_CASE("Options, errormodel_accepts_no_argument")
{
    option_test c({ "cagee", "-e" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.use_error_model);
    CHECK(actual.error_model_file_path.empty());
}

TEST_CASE("Options, replicate_map needs argument")
{
    option_test c({ "cagee", "--replicate_map", "map.tsv"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK(!actual.replicate_model_file_path.empty());
}

TEST_CASE("Options: zero_root_familes")
{
    input_parameters by_default;
    CHECK_FALSE(by_default.exclude_zero_root_transcripts);

    option_test c({ "cagee", "-z" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.exclude_zero_root_transcripts);
}

TEST_CASE("Options: count_all_changes")
{
    input_parameters by_default;
    CHECK_FALSE(by_default.exclude_zero_root_transcripts);

    option_test c({ "cagee", "--count_all_changes" });

    auto actual = read_arguments(c.argc, c.values);
    CHECK(actual.count_all_changes);
}

TEST_CASE("Options: must_specify_sigma_for_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    CHECK_THROWS_WITH(params.check_input(), "Cannot simulate without initial sigma values");
}

TEST_CASE("Options: must_specify sigma and_input_file_for_estimator")
{
    input_parameters params;
    params.fixed_sigma = 0.05;
    CHECK_THROWS_WITH(params.check_input(), "Options -l and -i must both be provided an argument.");
}

TEST_CASE("Options: must_specify_alpha_for_gamma_simulation")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_sigma = 0.05;
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
    params.fixed_sigma = 0.01;
    params.is_simulating = true;
    params.check_input();
    CHECK(true);
}

TEST_CASE("Options: check_input_does_not_throw_when_simulating_with_multiple_sigmas")
{
    input_parameters params;
    params.is_simulating = true;
    params.fixed_multiple_sigmas = "0.01,0.05";
    params.sigma_tree_file_path = "./tree";
    params.check_input();
    CHECK(true);
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
    params.fixed_sigma = 10;
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
    params.fixed_sigma = 10;
    params.is_simulating = true;
    params.check_input();
    CHECK(true);

    params.input_file_path = "transcripts.txt";
    CHECK_THROWS_WITH_AS(params.check_input(), "A families file was provided while simulating", runtime_error);

}

TEST_CASE("Prior params defaults to gamma if not unbounded")
{
    input_parameters params;
    CHECK_EQ("gamma:0.375:1600.0", params.prior_params_or_default());
}

TEST_CASE("Prior params defaults to fisher if unbounded")
{
    input_parameters params;
    params.input_file_has_ratios = true;
    CHECK_EQ("fisher:0.75:0.75", params.prior_params_or_default());
}

TEST_CASE("Options: n_gamma_cats")
{
    input_parameters by_default;
    CHECK_EQ(1, by_default.n_gamma_cats);

    option_test c({ "cagee", "--n_gamma_cats", "5"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(5, actual.n_gamma_cats);
}

TEST_CASE("Options: fixed_alpha")
{
    input_parameters by_default;
    CHECK_EQ(-1, by_default.fixed_alpha);

    option_test c({ "cagee", "--n_gamma_cats", "2", "--fixed_alpha", "5"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(5, actual.fixed_alpha);
}

TEST_CASE("Options: freerate_model")
{
    input_parameters by_default;
    CHECK_EQ("", by_default.free_rate);

    option_test c({ "cagee", "--free_rate"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("global", actual.free_rate);
}

TEST_CASE("Options: freerate paired model")
{
    input_parameters by_default;
    CHECK_EQ("", by_default.free_rate);

    option_test c({ "cagee", "--free_rate", "paired"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ("paired", actual.free_rate);
}

TEST_CASE("Options: If Freerate is set, max_iterations defaults to 2")
{
    input_parameters by_default;
    CHECK_EQ(300, by_default.optimizer_params.max_iterations);

    option_test c({ "cagee", "--free_rate", "paired"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(2, actual.optimizer_params.max_iterations);
}

TEST_CASE("Options: If Freerate is set, user can still set max iterations")
{
    input_parameters by_default;
    CHECK_EQ(300, by_default.optimizer_params.max_iterations);

    option_test c({ "cagee", "--free_rate", "paired", "--optimizer_iterations", "10"});

    auto actual = read_arguments(c.argc, c.values);
    CHECK_EQ(10, actual.optimizer_params.max_iterations);
}
