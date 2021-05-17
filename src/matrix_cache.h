#ifndef MATRIX_CACHE_H
#define MATRIX_CACHE_H

#include <map>
#include <vector>
#include <set>

#include <Eigen/Dense>

class matrix_cache_key {
    size_t _size;
    long _lambda;
    long _branch_length;
public:
    matrix_cache_key(int size, double some_lambda, double some_branch_length) :
        _size(size),
        _lambda(long(some_lambda * 1000000000)),    // keep 9 significant digits
        _branch_length(long(some_branch_length * 1000)) {} // keep 3 significant digits

    bool operator<(const matrix_cache_key &o) const {
        return std::tie(_size, _branch_length, _lambda) < std::tie(o._size, o._branch_length, o._lambda);
    }
    double lambda() const {
        return double(_lambda) / 1000000000.0;
    }
    double branch_length() const {
        return double(_branch_length) / 1000.0;
    }
};

//! Computation of the probabilities of moving from a family size (parent) to another (child)
/*!
Contains a map (_cache) that serves as a hash table to store precalculated values.
If the given parameters have already been calculated, will return the cached value rather than calculating the value again.
*/
class matrix_cache {
private:
    std::map<matrix_cache_key, Eigen::MatrixXd*> _matrix_cache; //!< nested map that stores transition probabilities for a given lambda and branch_length (outer), then for a given parent and child size (inner)
    int _matrix_size;
public:
    void precalculate_matrices(const std::vector<double>& lambdas, const std::set<double>& branch_lengths);
    const Eigen::MatrixXd* get_matrix(double branch_length, double lambda) const;

    int get_cache_size() const {
        return _matrix_cache.size();
    }

    int get_matrix_size() const {
        return _matrix_size;
    }

    matrix_cache();
    ~matrix_cache();

    friend std::ostream& operator<<(std::ostream& ost, matrix_cache& c);

};

std::ostream& operator<<(std::ostream& ost, matrix_cache& c);
#endif
