#ifndef root_equilibrium_distribution_h
#define root_equilibrium_distribution_h

#include <vector>
#include <map>

class clade;
struct input_parameters;
class root_distribution;
class gene_transcript;

class root_equilibrium_distribution
{
    std::vector<int> _vectorized_distribution;
    std::vector<double> _frequency_percentage;
    double _fixed_root_value = -1;
    void build_percentages();
    void create_from_poisson(double poisson_lambda, size_t num_values);
public:
    /// Create a distribution matching that in the map
    root_equilibrium_distribution(const std::map<int, int>& root_distribution);

    /// Create a uniform distribution up to the given size
    root_equilibrium_distribution(size_t max_size);

    root_equilibrium_distribution(double fixed_root_value);

    /// Create a Poisson distribution with the given lambda
    root_equilibrium_distribution(double poisson_lambda, size_t num_values);

    /// Estimate a Poisson distribution from the given families
    root_equilibrium_distribution(const std::vector<gene_transcript>& gene_families, size_t num_values);

    /// Move constructor
    root_equilibrium_distribution(root_equilibrium_distribution&& other) noexcept
    {
        *this = std::move(other);
    }

    /// return the prior probability of root size being n based on the given root distribution
    float compute(const gene_transcript& t, size_t n) const;

    int select_root_size(int family_number) const;

    // this is the move assignment operator
    root_equilibrium_distribution& operator=(root_equilibrium_distribution&& other) noexcept
    {
        _vectorized_distribution = std::move(other._vectorized_distribution);
        _frequency_percentage = std::move(other._frequency_percentage);
        _fixed_root_value = other._fixed_root_value;
        return *this;
    }

    void resize(size_t new_size);

};

#endif
