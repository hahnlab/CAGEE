#ifndef MATRIX_CACHE_H
#define MATRIX_CACHE_H

#include <map>
#include <vector>
#include <set>

#include <Eigen/Dense>

class sigma;
class DiffMat;

class matrix_cache_key {
    long _branch_length;
    long _sigma;
    const double SIGNIFICANT_DIGITS = double(1e9);

public:
    matrix_cache_key(double sigma, double branch_length)
    {
        _branch_length = long(branch_length * SIGNIFICANT_DIGITS);
        _sigma = long(sigma * SIGNIFICANT_DIGITS);
    }

    bool operator<(const matrix_cache_key& o) const {
        return std::tie(_branch_length, _sigma) < std::tie(o._branch_length, o._sigma);
    }
    double branch_length() const {
        return double(_branch_length) / SIGNIFICANT_DIGITS;
    }

    double sigma() const {
        return double(_sigma) / SIGNIFICANT_DIGITS;
    }
};

using boundaries = std::pair<double, double>;

//! Computation of the probabilities of moving from a family size (parent) to another (child)
/*!
Contains a map (_cache) that serves as a hash table to store precalculated values.
If the given parameters have already been calculated, will return the cached value rather than calculating the value again.
*/
class matrix_cache {
private:
    std::map<matrix_cache_key, Eigen::MatrixXd> _matrix_cache; //!< nested map that stores transition probabilities for a given lambda and branch_length (outer), then for a given parent and child size (inner)
public:
    void precalculate_matrices(const std::vector<double>& sigmas, const std::set<double>& branch_lengths, boundaries bounds);
    const Eigen::MatrixXd& get_matrix(double branch_length, double sigma) const;
    void set_matrix(double branch_length, double sigma, const Eigen::MatrixXd& m);

    Eigen::VectorXd create_vector() const;

    int get_cache_size() const {
        return _matrix_cache.size();
    }

    static DiffMat* diffusion_matrix;
    static void initialize(int Npts);
    friend std::ostream& operator<<(std::ostream& ost, matrix_cache& c);

};

std::ostream& operator<<(std::ostream& ost, matrix_cache& c);
#endif
