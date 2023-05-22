#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>
#include <map>

class clade;
class gene_transcript;

class ss_resolver;

class sigma_squared {
private:
    ss_resolver* _p_resolver;
    std::vector<double> _values;

public:
    sigma_squared(double lam);
    sigma_squared(ss_resolver *p_resolver, size_t sz);
    sigma_squared(ss_resolver* p_resolver, const std::vector<double>& v);

    ~sigma_squared();

    sigma_squared* multiply(double factor) const;
    void update(const double* values) ;
    int count() const;

    std::string to_string() const;
    double get_named_value(const clade* c, const gene_transcript& t) const;

    bool is_valid() const;

    std::vector<double> get_values() const;

    sigma_squared* clone() const;

    static sigma_squared* create(clade* p_sigma_tree, const std::vector<std::string>& sample_groups);
    static sigma_squared* create(clade* p_sigma_tree, const std::vector<double>& values);

};

inline std::ostream& operator<<(std::ostream& ost, const sigma_squared& sigma)
{
    ost << sigma.to_string();
    return ost;
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */


#endif
