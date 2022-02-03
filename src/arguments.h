#include <random>

#include "optimizer.h"
#include "root_equilibrium_distribution.h"

struct input_parameters {
public:
    std::string input_file_path;
    std::string error_model_file_path;
    std::string output_prefix;
    std::string tree_file_path;
    std::string lambda_tree_file_path;
    std::string fixed_multiple_lambdas;
    std::string log_config_file;
    double fixed_lambda = 0.0;
    double fixed_alpha = -1.0;
    double pvalue = 0.05;
    bool is_simulating = false;
    int nsims = 0;
    int n_gamma_cats = 1;
    bool exclude_zero_root_families = false;
    bool lambda_per_family = false;
    bool use_error_model = false;
    int verbose_logging_level = 0;
    int cores = 0;
    optimizer_parameters optimizer_params;
    std::string rootdist_params;
    std::vector<std::string> sample_groups;
    std::gamma_distribution<double> prior = std::gamma_distribution<double>(0.75, 30.0);
    bool help = false;

    //! Check calls
    void check_input();

    std::string rootdist_params_or_default() const
    {
        return rootdist_params.empty() ? "gamma:0.75:30.0" : rootdist_params;
    }
};

input_parameters read_arguments(int argc, char* const argv[]);
