#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "error_model.h"

using namespace std;

class error_model {

private:

    string _model_type = "normal";
    int _vector_length = 200;
    long _upper_bound = 1e6;
    double _elem_width = _upper_bound / _vector_length;
    vector<double> _elem_vals = vector<double>(200);
    double _likelihood_cutoff = 0.000001
    
    // vals for normal mean vs. variance fit
    double _a = 0.314021;
    double _r = -0.553507;
    double _inv_sqrt_2pi = 1 / sqrt(2 * M_PI);

    void set_elem_vals();
    double calc_variance(double log_counts);
    double calc_density(double element_val, double log_counts);

public:
    error_model();
    error_model(int vector_length, long upper_bound);
    error_model(string model_type, int vector_length, long upper_bound, vector<double> model_params);
    Eigen::VectorXd get_error_vector(double log_counts);
    void print_info();
};

// ORIGINAL (CAFE) ERROR MODEL HEADER:
/*
#ifndef ERROR_MODEL_H
#define ERROR_MODEL_H

#include <vector>
#include <string>
#include <map>

class error_model {
private:
    size_t _max_family_size;

    std::vector<int> _deviations; //!< Deviations from the true gene family (e.g., -1 0 1)

    std::vector<std::vector<double> > _error_dists; //!< Each vector element will be a gene family size; the vector of doubles inside (e.g., 0.1 0.8 0.1) will be the probs of deviating from the true value

public:
    error_model();

    //! Set max family size for which deviations apply
    void set_max_family_size(size_t max_cnt);

    //! Set deviations
    void set_deviations(std::vector<std::string> deviations);

    //! Set deviation probability vector for a certain family size
    void set_probabilities(size_t fam_size, std::vector<double>);

    //! Get deviation probability vector for a certain family size
    std::vector<double> get_probs(size_t fam_size) const;

    size_t n_deviations() const {
        return _deviations.size();
    }

    size_t get_max_family_size() const {
        return _error_dists.size();
    }
    std::vector<double> get_epsilons() const;
    void replace_epsilons(std::map<double, double>* new_epsilons);
    void update_single_epsilon(double new_epsilon);

    friend void write_error_model_file(std::ostream& ost, error_model& errormodel);
};

#endif
*/
