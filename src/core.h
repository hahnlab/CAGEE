#ifndef CORE_H
#define CORE_H

#include <set>

#include "easylogging++.h"

class reconstruction;
class user_data;
class gene_transcript;
class matrix_cache;
class sigma_optimizer_scorer;
class clade;
class sigma_squared;
class error_model;
class prior;

struct input_parameters;

using boundaries = std::pair<double, double>;

class event_monitor : public el::Loggable
{
    std::map<std::string, int> failure_count;
    int attempts = 0;
    int rejects = 0;
public:
    virtual void log(el::base::type::ostream_t& os) const;

    void Event_InferenceAttempt_Started();
    void Event_InferenceAttempt_InvalidValues() { rejects++; }
    void Event_InferenceAttempt_Saturation(std::string transcript) { failure_count[transcript]++; }
};

/*! @brief Describes the actions that are taken when estimating or simulating data

    A Model represents a way to calculate or simulate values in the data.
*/
class model {
protected:
    std::ostream & _ost; 
    sigma_squared*_p_sigma;
    error_model* _p_error_model;
    std::vector<std::vector<int> > _rootdist_bins;

    /// Used to track gene transcripts with identical sizes
    std::vector<size_t> references;

    event_monitor _monitor;

public:
    model(sigma_squared* p_sigma,
        const std::vector<gene_transcript>* p_gene_transcripts,
        error_model *p_error_model);
    
    virtual ~model() {}
    
    sigma_squared* get_sigma() const {
        return _p_sigma;
    }

    virtual sigma_squared* get_simulation_sigma();

    virtual double infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma) = 0;  // return vector of likelihoods
    
    void write_vital_statistics(std::ostream& ost, const clade *p_tree, double final_likelihood, const input_parameters& p);
    void write_error_model(int max_transcript_size, std::ostream& ost) const;

    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc) = 0;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups) = 0;

    virtual void write_extra_vital_statistics(std::ostream& ost) {}

    const event_monitor& get_monitor() { return _monitor;  }
};

std::vector<size_t> build_reference_list(const std::vector<gene_transcript>& transcripts);
bool verify_results(const std::vector<size_t>& references, const std::vector<std::vector<double>>& partial_likelihoods, int& error);

std::vector<model *> build_models(const input_parameters& my_input_parameters, user_data& user_data);

inline std::string filename(std::string base, std::string suffix, std::string extension)
{
    return (suffix.empty() ? std::string("results") : suffix) + "/" + base + "." + extension;
}

inline std::string filename(std::string base, std::string suffix)
{
    return filename(base, suffix, "txt");
}


#endif /* CORE_H */

