#include <random>

#include "optimizer.h"

struct input_parameters {
public:
    std::string input_file_path;
    // std::string error_model_file_path;
    std::string replicate_model_file_path;
    std::string output_prefix;
    std::string tree_file_path;
    std::string sigma_tree_file_path;
    std::string fixed_multiple_sigmas;
    std::string log_config_file;
    std::string command_line;
    double fixed_sigma = 0.0;
    double fixed_alpha = -1.0;
    bool is_simulating = false;
    int nsims = 0;
    int n_gamma_cats = 1;
    bool exclude_zero_root_transcripts = false;
    bool use_error_model = false;
    bool count_all_changes = false;
    int verbose_logging_level = 0;
    int cores = 0;
    optimizer_parameters optimizer_params;
    std::string rootdist_params;
    std::vector<std::string> sample_groups;
    std::string prior;
    bool help = false;
    int discretization_size = 200;
    bool input_file_has_ratios = false;

    //! Check calls
    void check_input();

    std::string rootdist_params_or_default() const
    {
        return rootdist_params.empty() ? "gamma:0.375:1600.0" : rootdist_params;
    }

    std::string prior_params_or_default() const
    {
        if (!prior.empty())
            return prior;

        return input_file_has_ratios ? "fisher:0.75:0.75" : "gamma:0.375:1600.0";
    }
};

input_parameters read_arguments(int argc, char* const argv[]);
