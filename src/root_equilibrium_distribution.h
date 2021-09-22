#ifndef root_equilibrium_distribution_h
#define root_equilibrium_distribution_h

#include <vector>
#include <map>

class clade;
struct input_parameters;
class root_distribution;
class gene_transcript;

struct rootdist_options
{
    rootdist_options() : type(gamma)
    {

    }
    rootdist_options(std::string cfg);
    enum rootdist_type { gamma, fixed, file };
    double gamma_alpha = 0.25, gamma_beta = 10.0;
    double fixed_value = 0.0;
    std::string filename;
    rootdist_type type;
};

class root_equilibrium_distribution
{
public:
    /// return the prior probability of root size being n based on the given root distribution
    virtual float compute(const gene_transcript& t, size_t n) const = 0;
    virtual void resize(size_t new_size) = 0;
    virtual float select_root_value(int family_number) const = 0;
};

class root_distribution_fixed : public root_equilibrium_distribution
{
    double _fixed_root_value = -1;
public:
    root_distribution_fixed(double fixed_root_value) : _fixed_root_value(fixed_root_value) {}

    float compute(const gene_transcript& t, size_t n) const;
    void resize(size_t new_size);
    float select_root_value(int family_number) const;
};

class root_distribution_uniform : public root_equilibrium_distribution
{
    int _max_value;
public:
    root_distribution_uniform(int max_value) : _max_value(max_value) {}

    float compute(const gene_transcript& t, size_t n) const;
    void resize(size_t new_size);
    float select_root_value(int family_number) const;
};

class root_distribution_poisson : public root_equilibrium_distribution
{
    std::vector<float> _vectorized_distribution;
    std::vector<double> _frequency_percentage;
public:
    root_distribution_poisson(double poisson_lambda, size_t num_values);
    root_distribution_poisson(const std::vector<gene_transcript>& gene_families, size_t num_values);

    float compute(const gene_transcript& t, size_t n) const;
    void resize(size_t new_size);
    float select_root_value(int family_number) const;
}; 

class root_distribution_specific : public root_equilibrium_distribution
{
    std::vector<float> _vectorized_distribution;
    std::vector<double> _frequency_percentage;
public:
    root_distribution_specific(std::map<int, float> distribution);
    float compute(const gene_transcript& t, size_t n) const;
    void resize(size_t new_size);
    float select_root_value(int family_number) const;
};

class root_distribution_gamma : public root_equilibrium_distribution
{
    std::vector<float> _vectorized_distribution;
    std::vector<double> _frequency_percentage;
public:
    root_distribution_gamma(double alpha, double beta, size_t num_values);
    root_distribution_gamma(const std::vector<gene_transcript>& gene_families, size_t num_values);

    float compute(const gene_transcript& t, size_t n) const;
    void resize(size_t new_size);
    float select_root_value(int family_number) const;
};

#endif
