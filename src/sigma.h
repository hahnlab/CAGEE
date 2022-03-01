#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>
#include <map>

class clade;
class gene_transcript;

enum class sigma_type { uniform, lineage_specific, sample_specific };

class sigma_squared {
private:
    std::map<std::string, int> _node_name_to_sigma_index;
    std::vector<double> _values;
    sigma_type _type = sigma_type::uniform;
public:
    sigma_squared(double lam)
    {
        _values.push_back(lam);
    }

    sigma_squared(std::map<std::string, int> nodename_index_map, std::vector<double> value_vector, sigma_type t) :
        _node_name_to_sigma_index(nodename_index_map), _values(value_vector), _type(t) { } //!< Constructor

    sigma_squared* multiply(double factor) const
    {
        auto npi = _values;

        for (auto& i : npi)
            i *= factor;

        return new sigma_squared(_node_name_to_sigma_index, npi, _type);
    }
    void update(const double* values) ;
    int count() const  {
        return _values.size();
    }
    std::string to_string() const;
    double get_value_for_clade(const clade* c) const ;
    double get_named_value(const clade* c, const gene_transcript& t) const;

    bool is_valid() const ;

    std::vector<double> get_values() const {
        return _values;
    }
    virtual sigma_squared* clone() const  {
        return new sigma_squared(_node_name_to_sigma_index, _values, _type);
    }

};

inline std::ostream& operator<<(std::ostream& ost, const sigma_squared& sigma)
{
    ost << sigma.to_string();
    return ost;
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */


#endif
