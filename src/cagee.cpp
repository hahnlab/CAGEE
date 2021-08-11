#include <map>
#include <random>
#include <numeric>
#include <fstream>
#include <omp.h>


#include "easylogging++.h"

#include "execute.h"
#include "simulator.h"

#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "core.h"
#include "arguments.h"


using namespace std;

void init_lgamma_cache();

action* get_executor(input_parameters& user_input, user_data& data)
{
    if (user_input.is_simulating) {
        return new simulator(data, user_input);
    }
    else
    {
        return new estimator(data, user_input);
    }

    return NULL;
}

void show_help()
{
    const char *text = ""
        "\n\nUsage: cagee [options]\n\n"
        "CAGEE is a software that provides a statistical foundation for evolutionary inferences about changes in gene family size.\n "
        "The program employs a birth and death process to model gene gain and loss across a user-specified phylogenetic tree,\n "
        "thus accounting for the species phylogenetic history. The distribution of family sizes generated under this model can\n "
        "provide a basis for assessing the significance of the observed family size differences among taxa.\n\n"
        "OPTIONS:\n"
        "  Required Options:\n"
        "   --infile, -i\t\t\tPath to tab delimited gene families file to be analyzed - Required for estimation.\n"
        "   --tree, -t\t\t\tPath to file containing newick formatted tree - Required for estimation.\n\n"

        "  Common Options:\n"
        "   --help, -h\t\t\tThis help menu.\n"
        "   --cores, -c\t\t\tNumber of processing cores to use, requires an integer argument. Default=All available cores.\n"
        "   --error_model, -e\t\tRun with no file name to estimate the global error model file. This file can be provided\n \t\t\t\t    in subsequent runs by providing the path to the Error model file with no spaces (e.g. -eBase_error_model.txt).\n"
//        "   --logfile, -g\t\t\tFilename to which run log will be written, Default=cafe.log (requires -L option to specify log configuration file).\n"
        "   --log_config, -L\t\tTurn on logging, provide name of the configuration file for logging (see example log.config file).\n"
        "   --n_gamma_cats, -k\t\tNumber of gamma rate categories to use. If specified, the Gamma model will be used to run\n \t\t\t\t    calculations, otherwise the Base model will be used.\n"
        "   --output_prefix, -o\t\tOutput directory - Name of directory automatically created for output. Default=results.\n"
        "   --lambda_tree, -y\t\tPath to lambda tree, for use with multiple lambdas.\n"
        "   --fixed_alpha, -a\t\tAlpha value of the discrete gamma distribution to use in category calculations. If not\n \t\t\t\t    specified, the alpha parameter will be estimated by maximum likelihood.\n"
        "   --fixed_lambda, -l\t\tValue (between 0 and 1) for a single user provided lambda value, otherwise lambda is estimated.\n"
        "   --fixed_multiple_lambdas, -m\tMultiple lambda values, comma separated.\n"
        "   --poisson, -p\t\tUse a Poisson distribution for the root frequency distribution. Without specifying this, a\n \t\t\t\t    uniform distribution will be used. A value can be specified -p10 (no space) or --poisson = 10,\n \t\t\t\t    otherwise the distribution will be estimated from the gene families.\n"
        "   --simulate, -s\t\tSimulate families. Optionally provide the number of simulations to generate (-s100 no space, or --simulate = 100)\n"
        "   --rootdist, -f\t\tPath to root distribution file for simulating datasets.\n\n"

        "  Less Common Options:\n"
        "   --lambda_per_family, -b\tEstimate lambda by family (for testing purposes only).\n"
//        "   --chisquare_compare, -r\tChi square compare.\n"
        "   --pvalue, -P\t\t\tP-value to use for determining significance of family size change, Default=0.05.\n"
        "   --zero_root, -z\t\tExclude gene families that don't exist at the root, not recommended.\n"
        "   --Expansion, -E\t\tExpansion parameter for Nelder-Mead optimizer, Default=2.\n"
		"   --Iterations, -I\t\tMaximum number of iterations that will be performed in \n \t\t\t\t    lambda search. Default=300 (increase this number if likelihood is still improving when limit is hit).\n"
        "   --Reflection, -R\t\tReflection parameter for Nelder-Mead optimizer, Default=1.\n\n\n";

        std::cout << text;
}

/// The main function. Evaluates arguments, calls processes
/// \callgraph
int cagee(int argc, char *const argv[]) {
    init_lgamma_cache();

    try {
        input_parameters user_input = read_arguments(argc, argv);

        if (user_input.help)
        {
            show_help();
            return 0;
        }
        if (user_input.cores > 0)
        {
            omp_set_num_threads(user_input.cores);
        }
        user_data data;
        data.read_datafiles(user_input);

        auto cmd = std::accumulate(argv, argv + argc, std::string(), [](std::string x, std::string y) { return x + y + " "; });
        LOG(INFO) << "\nCommand line: " << cmd;

        if (user_input.exclude_zero_root_families)
        {
            auto rem = std::remove_if(data.gene_families.begin(), data.gene_families.end(), [&data](const gene_transcript& fam) {
                return !fam.exists_at_root(data.p_tree);
            });

            int fmsize = data.gene_families.size();
            data.gene_families.erase(rem, data.gene_families.end());
            LOG(INFO) << "\nFiltering families not present at the root from: " << fmsize << " to " << data.gene_families.size();

        }

        data.create_prior(user_input);

        // When computing or simulating, only base or gamma model is used. When estimating, base and gamma model are used (to do: compare base and gamma w/ LRT)
        // Build model takes care of -f
        vector<model *> models = build_models(user_input, data);

        unique_ptr<action> act(get_executor(user_input, data));
        if (act)
        {
            act->execute(models);
        }

        return 0;
    }
    catch (runtime_error& err) {
        LOG(ERROR) << err.what() << endl;
        return EXIT_FAILURE;
    }
} 
