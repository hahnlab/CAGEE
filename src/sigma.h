#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>
#include <map>

class clade;

class sigma {
private:
    std::map<std::string, int> _node_name_to_lambda_index;
    std::vector<double> _lambdas;

public:
    sigma(double lam)
    {
        _lambdas.push_back(lam);
    }

    sigma(std::map<std::string, int> nodename_index_map, std::vector<double> lambda_vector) :
        _node_name_to_lambda_index(nodename_index_map), _lambdas(lambda_vector) { } //!< Constructor

    sigma* multiply(double factor) const
    {
        auto npi = _lambdas;

        for (auto& i : npi)
            i *= factor;

        return new sigma(_node_name_to_lambda_index, npi);
    }
    void update(const double* values) ;
    int count() const  {
        return _lambdas.size();
    }
    std::string to_string() const;
    double get_value_for_clade(const clade* c) const ;

    bool is_valid() const ;

    std::vector<double> get_lambdas() const {
        return _lambdas;
    }
    virtual sigma* clone() const  {
        return new sigma(_node_name_to_lambda_index, _lambdas);
    }

};

inline std::ostream& operator<<(std::ostream& ost, const sigma& sigma)
{
    ost << sigma.to_string();
    return ost;
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */


#endif
