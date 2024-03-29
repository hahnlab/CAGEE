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
    sigma_squared(const sigma_squared& other);
    sigma_squared(ss_resolver *p_resolver, size_t sz);
    sigma_squared(ss_resolver* p_resolver, const std::vector<double>& v);

    // multiplier constructor
    sigma_squared(const sigma_squared& other, double mutiplier);

    ~sigma_squared();

    void update(const double* values) ;
    int count() const;

    double get_named_value(const clade* c, const gene_transcript& t) const;

    bool is_valid() const;

    std::vector<double> get_values() const;

//    sigma_squared* clone() const;

    static sigma_squared* create(clade* p_sigma_tree, const std::vector<std::string>& sample_groups);
    static sigma_squared* create(clade* p_sigma_tree, const std::vector<double>& values);

    friend std::ostream& operator<<(std::ostream& ost, const sigma_squared& sigma);
};

std::ostream& operator<<(std::ostream& ost, const sigma_squared& sigma);

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */


#endif
