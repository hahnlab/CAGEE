#include <map>
#include <random>
#include <numeric>
#include <fstream>
#include <omp.h>


#include "easylogging++.h"

#include "execute.h"
#include "simulator.h"
#include "matrix_cache.h"
#include "user_data.h"
#include "root_equilibrium_distribution.h"
#include "core.h"
#include "arguments.h"

using namespace std;

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

/// The main function. Evaluates arguments, calls processes
/// \callgraph
int cagee(int argc, char *const argv[]) {
    try {
        input_parameters user_input = read_arguments(argc, argv);

        if (user_input.help)
        {
            return 0;
        }
        if (user_input.cores > 0)
        {
            omp_set_num_threads(user_input.cores);
        }
        user_data data;
        data.read_datafiles(user_input);

        LOG(INFO) << "Command line: " << user_input.command_line;

        if (user_input.exclude_zero_root_transcripts)
        {
            auto rem = std::remove_if(data.gene_transcripts.begin(), data.gene_transcripts.end(), [&data](const gene_transcript& fam) {
                return !fam.exists_at_root(data.p_tree);
            });

            int fmsize = data.gene_transcripts.size();
            data.gene_transcripts.erase(rem, data.gene_transcripts.end());
            LOG(INFO) << "\nFiltering families not present at the root from: " << fmsize << " to " << data.gene_transcripts.size();

        }

        gene_transcript::remove_ungrouped_transcripts(user_input.sample_groups, data.gene_transcripts);

        matrix_cache::initialize(user_input.discretization_size);
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
