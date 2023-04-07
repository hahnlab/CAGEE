#include <algorithm>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <random>

#include <Eigen/Dense>

#include "matrix_cache.h"
#include "doctest.h"
#include "easylogging++.h"
#include "DiffMat.h"
#include "sigma.h"

extern std::mt19937 randomizer_engine;
using namespace Eigen;
using namespace std;

DiffMat* matrix_cache::diffusion_matrix = nullptr;

std::ostream& operator<<(std::ostream& ost, const matrix_cache_key& k)
{
    ost << "(" << k.branch_length() << ", " << k.branch_length() << ");";
    return ost;
}

void matrix_cache::initialize(int Npts)
{
    LOG(INFO) << "Discretization range set to " << Npts;
    diffusion_matrix = new DiffMat(Npts);
    diffusion_matrix->create_or_read_eigenvectors();
}

const Eigen::MatrixXd& matrix_cache::get_matrix(double branch_length, double sigma) const
{
    // cout << "Matrix request " << size << "," << branch_length << "," << lambda << endl;

    matrix_cache_key key(sigma, branch_length);
    if (_matrix_cache.find(key) != _matrix_cache.end())
    {
       return _matrix_cache.at(key);
    }
    else
    {
        ostringstream ost;
        ost << "Failed to find matrix for " << key;
        throw std::runtime_error(ost.str());
    }
}

void matrix_cache::precalculate_matrices(const std::vector<double>& squared_sigmas, const std::set<double>& branch_lengths, boundaries bounds)
{
    // build a list of required matrices
    vector<matrix_cache_key> keys;
    for (double branch_length : branch_lengths)
    {
        for (double sigsqd : squared_sigmas)
        {
            matrix_cache_key key(sigsqd, branch_length);
            if (_matrix_cache.find(key) == _matrix_cache.end())
            {
                keys.push_back(key);
            }
        }
    }

    // calculate matrices in parallel
    size_t i = 0;
    size_t num_keys = keys.size();
    vector<double> vBranches(keys.size());
    vector<double> vSigSqd(keys.size());
    transform(keys.begin(), keys.end(), vBranches.begin(), [](matrix_cache_key k) { return k.branch_length(); });
    transform(keys.begin(), keys.end(), vSigSqd.begin(), [](matrix_cache_key k) { return k.sigma() / 2; });
    auto matrices = ConvProp_bounds_batched(vBranches, vSigSqd, *diffusion_matrix, bounds);
    for (i = 0; i < num_keys; ++i)
    {
        _matrix_cache[keys[i]] = matrices[i];
    }

}

void matrix_cache::set_matrix(double branch_length, double sigma, const Eigen::MatrixXd& m)
{
    matrix_cache_key key(sigma, branch_length);
    _matrix_cache[key] = m;
}

Eigen::VectorXd matrix_cache::create_vector() const 
{ 
    return Eigen::VectorXd::Zero(diffusion_matrix->Diff.rows()); 
}

std::ostream& operator<<(std::ostream& ost, matrix_cache& c)
{
    ost << c.get_cache_size() << " matrices. Keys: ";
    for (auto& kv : c._matrix_cache)
    {
        ost << kv.first;
    }
    return ost;
}


TEST_CASE(" matrix_cache_key handles floating point imprecision")
{
    set<matrix_cache_key> keys;
    double t = 0.3;
    for (int i = 0; i < 31; i++)
    {
        t += 0.1;
        matrix_cache_key key(0.01, t);
        keys.insert(key);
    }
    CHECK_EQ(31, keys.size());

    matrix_cache_key key(0.01, 1.4);
    CHECK_EQ(1, keys.count(key));
}

TEST_CASE("get_matrix returns correct matrix based on key")
{
    matrix_cache cache;
    Matrix3d m1;
    cache.set_matrix(3.2642504711034, 85, m1);

    Matrix2d m2;
    cache.set_matrix(3.1010379475482, 85, m2);

    CHECK_EQ(4, cache.get_matrix(3.1010379475482, 85).size());
    CHECK_EQ(9, cache.get_matrix(3.2642504711034, 85).size());
}

TEST_CASE("matrix_cache_key checks sigma for less than")
{
    matrix_cache_key key(3.2642504711034, 44);
    matrix_cache_key key2(3.1010379475482, 44);

    CHECK(key2 < key);
}

TEST_CASE("Probability:matrices_take_fractional_branch_lengths_into_account")
{
    // This test is not right. It should check that the matrix contains different values for the two different branch lengths
    sigma_squared sigsq(9.1);
    matrix_cache calc;
    std::set<double> branch_lengths{ 68, 68.7105 };
    calc.precalculate_matrices(sigsq.get_values(), branch_lengths, boundaries(0,3));
    CHECK_EQ(doctest::Approx(0.005), calc.get_matrix(68.7105, sigsq.get_values()[0])(100, 100)); // a value 
    CHECK_EQ(doctest::Approx(0.005), calc.get_matrix(68,      sigsq.get_values()[0])(100, 100));
}

TEST_CASE("Probability: probability_of_matrix")
{
    sigma_squared sigsq(0.05);
    matrix_cache calc;
    std::set<double> branch_lengths{ 5 };
    calc.precalculate_matrices(sigsq.get_values(), branch_lengths, boundaries(0,3));
    MatrixXd actual = calc.get_matrix(5, sigsq.get_values()[0]);
    CHECK_EQ(200, actual.rows());
    CHECK_EQ(200, actual.cols());

    CHECK_EQ(doctest::Approx(0.02187), actual(10, 10));
#if 0
    MatrixXd expected(5, 5);
    double values[5][5] = {
    {1,0,0,0,0},
    { 0.2,0.64,0.128,0.0256,0.00512 },
    { 0.04,0.256,0.4608,0.17408,0.0512 },
    { 0.008,0.0768,0.26112,0.36352,0.187392 },
    { 0.0016,0.02048,0.1024,0.249856,0.305562 } };
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            expected(i, j) = values[i][j];
    CHECK(actual == expected);

    // a second call should get the same results as the first
    actual = calc.get_matrix(5, sigsq.get_values()[0], 0);
    CHECK(actual == expected);
#endif
}

