#include <cassert>
#include <stdexcept>
//#include <utility>
#include "error_model.h"

using namespace std;

// set element values of likelihood vector for pdf input
void error_model::set_elem_vals() {
    _elem_vals = vector<double> (_vector_length);
    _elem_width = _upper_bound / _vector_length;
    _elem_vals[0] = 0.5 * _elem_width;
    for(size_t i = 1; i < _elem_vals.size(); i++){
        _elem_vals[i] = _elem_vals[i - 1] + _elem_width;
    }
}

// current normal log counts variance fit = a * exp(r * log_counts)
double error_model::calc_variance(double log_counts) const { 
    return _a * exp(_r * log_counts);
}

// normal pdf = (1 / s * sqrt(2 * pi)) * exp(-0.5 * ((x - mu) / s)^2)
double error_model::calc_density(double element_val, double log_counts) const {
    double s = sqrt(calc_variance(log_counts));
    double A = _inv_sqrt_2pi * s;
    double B = (element_val - log_counts) / s;
    return A * exp(-0.5 * B * B);
}

// for default normal with user-specified likelihood vector
error_model::error_model(int vector_length, double upper_bound)
    : _vector_length{vector_length}, _upper_bound{upper_bound}
{
    cout << "initialized default parametric error model" << endl;
    set_elem_vals();
    print_info();
}

// for default normal with custom model parameters, or future non-normal model
error_model::error_model(int vector_length, long upper_bound, model_params model_params)
    : _vector_length{vector_length}, _upper_bound{upper_bound}
{
    if(model_params.model_type == "normal") {
        _model_type = model_params.model_type;
        _a = model_params.a;
        _r = model_params.r;
        set_elem_vals();
        cout << "initialized parametric error model with custom parameters" << endl;
        print_info();
    } else {
        throw invalid_argument("only NORMAL parametric error model currently supported");
    }
}

Eigen::VectorXd error_model::get_error_vector(double log_counts) const {
    Eigen::VectorXd likelihood = Eigen::VectorXd::Zero(_vector_length);
    // for loop over vector elements and their rep values, add PDF density
    assert(_elem_vals.size() == likelihood.size());
    for(size_t i = 0; i < likelihood.size(); i++) {
        likelihood(i) = calc_density(_elem_vals[i], log_counts);
    }
    // renormalize total likelihood to 1
    likelihood = likelihood / likelihood.sum();
    // cutoff off very small values to zero, renormalize to 1
    return (_likelihood_cutoff < likelihood.array()).select(likelihood, 0.0f);
}

void error_model::print_info() {
    cout << "initialized " << _model_type << " error model with parameters:" << endl;
    cout << "  fitting variance = a * exp(r * log_counts)" << endl;
    cout << "      a = " << _a << endl;
    cout << "      r = " << _r << endl;
    cout << "  likelihood vector length = " << _vector_length << endl;
    cout << "  upper_bound = " << _upper_bound << endl;
    cout << "  likelihood vector element width = " << _elem_width << endl;
    cout << "  likelihood vector element values (for pdf):" << endl << "    ";
    cout.precision(4);
    for(auto val : _elem_vals) cout << val << " | ";
    cout << endl << endl;
}

// ORIGINAL (CAFE) ERROR MODEL DEFINITIONS
/*
#include <numeric>
#include <set>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "error_model.h"

using namespace std;

error_model::error_model()
{
    _deviations = { -1, 0, 1 };
}

void error_model::set_max_family_size(size_t max_cnt) {
    _max_family_size = max_cnt;
}

void error_model::set_deviations(std::vector<std::string> deviations) {
    _deviations.resize(deviations.size());
    transform(deviations.begin(), deviations.end(), _deviations.begin(), [](const string& s) {return std::stoi(s); });
}

inline bool is_nearly_equal(double x, double y)
{
    const double epsilon = 0.01;
    return std::abs(x - y) <= epsilon * std::abs(x);
}

void error_model::set_probabilities(size_t fam_size, std::vector<double> probs_deviation) {
    if ((fam_size == 0 || _error_dists.empty()) && !is_nearly_equal(probs_deviation[0], 0.0))
    {
        throw std::runtime_error("Cannot have a non-zero probability for family size 0 for negative deviation");
    }

    if (!is_nearly_equal(accumulate(probs_deviation.begin(), probs_deviation.end(), 0.0), 1.0))
    {
        throw std::runtime_error("Sum of probabilities must be equal to one");
    }

    if (_error_dists.empty())
        _error_dists.push_back(probs_deviation);

    if (_error_dists.size() <= fam_size)
    {
        _error_dists.resize(fam_size + 1, _error_dists.back());
    }
    _error_dists[fam_size] = probs_deviation; // fam_size starts at 0 at tips, so fam_size = index of vector
}

std::vector<double> error_model::get_probs(size_t fam_size) const {
    if (fam_size >= _error_dists.size() && fam_size <= _max_family_size)
        return _error_dists.back();

    return _error_dists[fam_size];
}

std::vector<double> error_model::get_epsilons() const {
    set<double> unique_values;
    for (auto& vec : _error_dists)
        unique_values.insert(vec.back());

    vector<double> result(unique_values.size());
    copy(unique_values.begin(), unique_values.end(), result.begin());
    return result;
}

// simple case where we have a single epsilon value in the tree
void error_model::update_single_epsilon(double new_epsilon)
{
    auto epsilons = get_epsilons();
    assert(epsilons.size() == 1);
    map<double, double> replacements;
    replacements[epsilons[0]] = new_epsilon;
    replace_epsilons(&replacements);
}

void error_model::replace_epsilons(std::map<double, double>* new_epsilons)
{
    vector<double> vec = _error_dists[0];
    assert(vec.size() == 3);
    for (auto kv : *new_epsilons)
    {
        if (is_nearly_equal(kv.first, vec.back()))
        {
            vec.back() = kv.second;
            vec[1] = 1 - kv.second;
            set_probabilities(0, vec);
        }
    }

    for (size_t i = 1; i < _error_dists.size(); ++i)
    {
        vector<double> vec = _error_dists[i];
        assert(vec.size() == 3);

        for (auto kv : *new_epsilons)
        {
            if (is_nearly_equal(kv.first, vec.back()))
            {
                vec.back() = kv.second;
                vec.front() = kv.second;
                vec[1] = 1 - (kv.second * 2);
                set_probabilities(i, vec);
            }
        }
    }
}
*/
