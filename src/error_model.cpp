#include <cassert>
#include <stdexcept>
#include "error_model.h"
#include "user_data.h"
#include "easylogging++.h"

using namespace std;

// set element values of likelihood vector for pdf input
void error_model::set_elem_vals() {
    _elem_vals = vector<double>(_vector_length);
    _elem_width = static_cast<double>(_upper_bound) / static_cast<double>(_vector_length);
    _elem_vals[0] = 0.5 * _elem_width;
    for(size_t i = 1; i < _elem_vals.size(); i++){
        _elem_vals[i] = _elem_vals[i - 1] + _elem_width;
    }
}

// current normal log counts variance fit = a * exp(r * log_counts)
double error_model::calc_variance(double log_counts) const { 
    double variance = _a * exp(_r * log_counts);
    // VLOG(9) << "calc_variance(" << log_counts << ") = " << variance;
    return variance;
}

// normal pdf = (1 / s * sqrt(2 * pi)) * exp(-0.5 * ((x - mu) / s)^2)
double error_model::calc_density(double element_val, double log_counts) const {
    double s = sqrt(calc_variance(log_counts));
    double A = _inv_sqrt_2pi * s;
    double B = (element_val - log_counts) / s;
    return A * exp(-0.5 * B * B);
}

// for default normal with user-specified likelihood vector
error_model::error_model(int vector_length, int upper_bound)
    : _vector_length{ vector_length }, _upper_bound{ upper_bound }
{
    LOG(INFO) << "initializing default parametric error model" << endl;
    set_elem_vals();
    print_info();
}

// for default normal with custom model parameters, or future non-normal model
error_model::error_model(int vector_length, int upper_bound, model_params model_params)
    : _vector_length{ vector_length }, _upper_bound{ upper_bound }
{
    if(model_params.model_type == "normal") {
        LOG(INFO) << "initializing parametric error model with custom parameters..." << endl;
        _model_type = model_params.model_type;
        _a = model_params.a;
        _r = model_params.r;
        set_elem_vals();
        print_info();
    } else {
        throw invalid_argument("only 'NORMAL' parametric error model currently supported");
    }
}

Eigen::VectorXd error_model::get_error_vector(double log_counts) const {
    Eigen::VectorXd likelihood = Eigen::VectorXd::Zero(_vector_length);
    // for loop over vector elements and their rep values, add PDF density
    assert(_elem_vals.size() == (size_t) likelihood.size());
    for(size_t i = 0; i < (size_t) likelihood.size(); i++) {
        likelihood(i) = calc_density(_elem_vals[i], log_counts);
    }
    // renormalize total likelihood to 1
    likelihood = likelihood / likelihood.sum();
    // cutoff off very small values to zero, renormalize to 1
    return (_likelihood_cutoff < likelihood.array()).select(likelihood, 0.0f);
}

void error_model::print_info() {
    LOG(INFO) << "initialized " << _model_type << " error model with parameters:";
    LOG(INFO) << "  fitting variance = a * exp(r * log_counts)";
    LOG(INFO) << "      a = " << _a;
    LOG(INFO) << "      r = " << _r;
    VLOG(7) << "  likelihood vector length = " << _vector_length;
    VLOG(7) << "  upper_bound = " << _upper_bound;
    VLOG(7) << "  likelihood vector element width = " << _elem_width;
    VLOG(7) << "  likelihood vector element values (for pdf):";
    cout << "    [";
    cout.precision(4);
    for(auto val : _elem_vals) cout << val << ", ";
    cout << "]" << endl << endl;
}
