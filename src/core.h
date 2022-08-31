#ifndef CORE_H
#define CORE_H

#include <set>
#include <random>

#include "easylogging++.h"

#include "clade.h"

class simulation_data;
class inference_process;
class gene_transcript_reconstructor;
class reconstruction;
class user_data;
class gene_transcript;
class matrix_cache;
class root_distribution_gamma;
class sigma_optimizer_scorer;
class sigma_squared;
class error_model;

struct input_parameters;

struct family_info_stash {
    family_info_stash() : lambda_multiplier(0.0), category_likelihood(0.0), family_likelihood(0.0), 
        posterior_probability(0.0), significant(false) {}
    family_info_stash(std::string fam, double lam, double cat_lh, double fam_lh, double pp, bool signif) : 
        family_id(fam), lambda_multiplier(lam), category_likelihood(cat_lh), family_likelihood(fam_lh),
        posterior_probability(pp), significant(signif) {}
    std::string family_id;
    double lambda_multiplier;
    double category_likelihood;
    double family_likelihood;
    double posterior_probability;
    bool significant;
};

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r);

class event_monitor : public el::Loggable
{
    std::map<std::string, int> failure_count;
    int attempts = 0;
    int rejects = 0;
public:
    virtual void log(el::base::type::ostream_t& os) const;

    void Event_InferenceAttempt_Started();
    void Event_InferenceAttempt_InvalidValues() { rejects++; }
    void Event_InferenceAttempt_Saturation(std::string family) { failure_count[family]++; }
};

/*! @brief Describes the actions that are taken when estimating or simulating data

    A Model represents a way to calculate or simulate values in the data.
*/
class model {
protected:
    std::ostream & _ost; 
    sigma_squared*_p_sigma;
    error_model* _p_error_model;
    std::vector<std::vector<int> > _rootdist_bins; // holds the distribution for each lambda bin

    /// Used to track gene families with identical species counts
    std::vector<size_t> references;

    std::vector<family_info_stash> results;

    event_monitor _monitor;

public:
    model(sigma_squared* p_lambda,
        const std::vector<gene_transcript>* p_gene_families,
        error_model *p_error_model);
    
    virtual ~model() {}
    
    sigma_squared* get_sigma() const {
        return _p_sigma;
    }

    //! Returns a lambda suitable for creating a simulated family. Default case is simply to return the lambda provided by the user.
    virtual sigma_squared* get_simulation_lambda();

    virtual double infer_family_likelihoods(const user_data& ud, const sigma_squared*p_lambda, const std::gamma_distribution<double>& prior) = 0;  // return vector of likelihoods
    
    virtual std::string name() const = 0;
    virtual void write_vital_statistics(std::ostream& ost, const clade *p_tree, double final_likelihood);
    void write_error_model(int max_family_size, std::ostream& ost) const;

    //! Based on the model parameters, attempts to reconstruct the most likely counts of each family at each node
    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc) = 0;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups, const std::gamma_distribution<double>& prior) = 0;

    const event_monitor& get_monitor() { return _monitor;  }
};

//! @brief Creates a list of families that are identical in all values
//!
//! With this information we can reduce the number of calculations required
//! and speed up the overall performance
std::vector<size_t> build_reference_list(const std::vector<gene_transcript>& families);

std::vector<model *> build_models(const input_parameters& my_input_parameters, user_data& user_data);

inline std::string filename(std::string base, std::string suffix, std::string extension)
{
    return (suffix.empty() ? std::string("results") : suffix) + "/" + base + "." + extension;
}

inline std::string filename(std::string base, std::string suffix)
{
    return filename(base, suffix, "txt");
}

int upper_bound_from_transcript_values(const std::vector<gene_transcript>& transcripts);

//! Create a sigma based on the sigma tree model the user passed.
/// Called when the user has provided no sigma value and one must
/// be estimated. If the sigma tree is NULL, uses a single
/// sigma; otherwise uses the number of unique sigmas in the provided
/// tree
sigma_squared* initialize_search_sigma(clade* p_sigma_tree, const std::vector<std::string>& sample_groups);

#endif /* CORE_H */

