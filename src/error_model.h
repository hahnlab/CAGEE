#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "user_data.h"

using namespace std;

class error_model {

private:

    string _model_type = "normal";
    int _vector_length;
    int _upper_bound;
    double _elem_width;
    vector<double> _elem_vals = vector<double>(200);
    double _likelihood_cutoff = 0.000001;
    
    // vals for log counts vs. variance fit
    double _a = 0.314021;
    double _r = -0.553507;
    double _inv_sqrt_2pi = 1 / sqrt(2 * M_PI);

    void set_elem_vals();
    double calc_variance(double log_counts) const;
    double calc_density(double element_val, double log_counts) const;


public:

    error_model(int vector_length, int upper_bound);
    error_model(int vector_length, int upper_bound, model_params model_params);
    Eigen::VectorXd get_error_vector(double log_counts) const;
    void print_info();
};


class fixed_error : public error_model {

private:

    double _variance;

public:

    fixed_error(int vector_length, int upper_bound, double variance);
    double calc_variance(double log_counts) const;
    void print_info();
};