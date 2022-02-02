#ifndef root_equilibrium_distribution_h
#define root_equilibrium_distribution_h

#include <vector>
#include <map>
#include <random>
#include <iosfwd>

class clade;
struct input_parameters;
class gene_transcript;

class root_equilibrium_distribution
{
    virtual float get_raw_root_value(int family_number) = 0;
public:
    /// return the prior probability of root value being n based on the given root distribution
    virtual void resize(size_t new_size) = 0;
    float select_root_value(int family_number);
};

class root_distribution_fixed : public root_equilibrium_distribution
{
    double _fixed_root_value = -1;
public:
    root_distribution_fixed(double fixed_root_value) : _fixed_root_value(fixed_root_value) {}

    void resize(size_t new_size)  override;
    float get_raw_root_value(int family_number) override;
};

class root_distribution_uniform : public root_equilibrium_distribution
{
    int _max_value;
public:
    root_distribution_uniform(int max_value) : _max_value(max_value) {}

    void resize(size_t new_size) override;
    float get_raw_root_value(int family_number) override;
};

class root_distribution_poisson : public root_equilibrium_distribution
{
    std::vector<float> _vectorized_distribution;
    std::vector<double> _frequency_percentage;
public:
    root_distribution_poisson(double poisson_lambda, size_t num_values);
    root_distribution_poisson(const std::vector<gene_transcript>& gene_families, size_t num_values);

    void resize(size_t new_size) override;
    float get_raw_root_value(int family_number) override;
}; 

class root_distribution_specific : public root_equilibrium_distribution
{
    std::vector<float> _vectorized_distribution;
    std::vector<double> _frequency_percentage;
public:
    root_distribution_specific(std::vector<std::pair<float, int>> distribution);
    void resize(size_t new_size) override;
    float get_raw_root_value(int family_number) override;
};

class root_distribution_gamma : public root_equilibrium_distribution
{
    std::gamma_distribution<double> _dist;
public:
    root_distribution_gamma(double alpha, double beta);

    void resize(size_t new_size) override;
    float get_raw_root_value(int family_number) override;
};

#endif
