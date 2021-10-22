#ifndef MATRIX_CACHE_H
#define MATRIX_CACHE_H

#include <map>
#include <vector>
#include <set>

#include "DiffMat.h"

class sigma;

class matrix_cache_key {
    std::pair<long, long> _bounds;
    long _branch_length;
    long _sigma;
public:
    matrix_cache_key(boundaries bounds, double sigma, double branch_length) :
        _sigma(long(sigma * 1000000000)),    // keep 9 significant digits
        _bounds(long(bounds.first * 1000000000), long(bounds.second * 1000000000)),    // keep 9 significant digits
        _branch_length(long(branch_length * 1000)) {} // keep 3 significant digits

    bool operator<(const matrix_cache_key& o) const {
        return std::tie(_bounds.first, _bounds.second, _branch_length, _sigma) < std::tie(o._bounds.first, o._bounds.second, o._branch_length, o._sigma);
    }
    boundaries bounds() const {
        return boundaries(double(_bounds.first / 1000000000.0), double(_bounds.second / 1000000000.0));
    }
    double branch_length() const {
        return double(_branch_length) / 1000.0;
    }

    double sigma() const {
        return double(_sigma) / 1000000000.0;
    }
};

//! Computation of the probabilities of moving from a family size (parent) to another (child)
/*!
Contains a map (_cache) that serves as a hash table to store precalculated values.
If the given parameters have already been calculated, will return the cached value rather than calculating the value again.
*/
class matrix_cache {
private:
    std::map<matrix_cache_key, Eigen::MatrixXd> _matrix_cache; //!< nested map that stores transition probabilities for a given lambda and branch_length (outer), then for a given parent and child size (inner)
public:
    void precalculate_matrices(const std::vector<double>& sigmas, const std::set<boundaries>& boundses, const std::set<double>& branch_lengths);
    const Eigen::MatrixXd& get_matrix(double branch_length, double sigma, boundaries bounds) const;
    void set_matrix(double branch_length, double sigma, boundaries bounds, const Eigen::MatrixXd& m);

    int get_cache_size() const {
        return _matrix_cache.size();
    }

    friend std::ostream& operator<<(std::ostream& ost, matrix_cache& c);

};

std::ostream& operator<<(std::ostream& ost, matrix_cache& c);
#endif
