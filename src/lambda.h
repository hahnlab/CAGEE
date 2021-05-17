#ifndef LAMBDA_H
#define LAMBDA_H

#include <vector>
#include <map>

class clade;
class matrix_cache;
class root_equilibrium_distribution;
class model;

struct FMinSearch;

/* START: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

//! Abstract class to hold lambda value.
/*!
  The main feature of this abstract class is its pure virtual member calculate_child_factor. This method is required because one or more lambda may be specified, and this determines how the pruner should compute the likelihoods. So this abstract class is what decides how the likelihood is computed for each descendant.
 */
class lambda {
public:
    virtual lambda* multiply(double factor) const = 0;
    virtual void update(const double* values) = 0;
    virtual int count() const = 0;
    virtual std::string to_string() const = 0;
    virtual double get_value_for_clade(const clade* c) const = 0;
    virtual bool is_valid() const = 0;
    virtual lambda* clone() const = 0;

    virtual ~lambda() {}
};

//! (lambda) Derived class 1: one lambda for whole tree
class single_lambda : public lambda {
private:
    double _lambda;

public:
    single_lambda(double lam) : _lambda(lam) { } //!< Constructor 
    double get_single_lambda() const { return _lambda; }

    virtual lambda* multiply(double factor) const override
    {
        return new single_lambda(_lambda * factor);
    }
    virtual void update(const double* values) override { _lambda = *values; }

    virtual int count() const override {
        return 1;
    }
    virtual std::string to_string() const override;
    virtual double get_value_for_clade(const clade* c) const override {
        return _lambda;
    }
    virtual bool is_valid() const override {
        return _lambda > 0;
    }
    virtual lambda* clone() const override {
        return new single_lambda(_lambda);
    }
};

//! (lambda) Derived class 2: multiple lambdas
class multiple_lambda : public lambda {
private:
    std::map<std::string, int> _node_name_to_lambda_index;
    std::vector<double> _lambdas;

public:
    multiple_lambda(std::map<std::string, int> nodename_index_map, std::vector<double> lambda_vector) :
        _node_name_to_lambda_index(nodename_index_map), _lambdas(lambda_vector) { } //!< Constructor

    virtual lambda* multiply(double factor) const override
    {
        auto npi = _lambdas;

        for (auto& i : npi)
            i *= factor;

        return new multiple_lambda(_node_name_to_lambda_index, npi);
    }
    virtual void update(const double* values) override;
    virtual int count() const override {
        return _lambdas.size();
    }
    virtual std::string to_string() const override;
    virtual double get_value_for_clade(const clade* c) const override;

    virtual bool is_valid() const override;

    std::vector<double> get_lambdas() const {
        return _lambdas;
    }
    virtual lambda* clone() const override {
        return new multiple_lambda(_node_name_to_lambda_index, _lambdas);
    }

};

inline std::ostream& operator<<(std::ostream& ost, const lambda& lambda)
{
    ost << lambda.to_string();
    return ost;
}

std::vector<double> get_lambda_values(const lambda* p_lambda);


/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */


#endif
