#include <algorithm>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <random>

#include <Eigen/Dense>

#include "matrix_cache.h"
#include "probability.h"
#include "doctest.h"
#include "easylogging++.h"
#include "DiffMat.h"
#include "sigma.h"

#ifdef BLAS_FOUND
#ifdef HAVE_OPENBLAS
#include "cblas.h"
#else
#include "mkl.h"
#endif
#endif

extern std::mt19937 randomizer_engine;
using namespace Eigen;
using namespace std;

std::ostream& operator<<(std::ostream& ost, const matrix_cache_key& k)
{
    ost << "(" << k.branch_length() << ", " << k.branch_length() << ", (" << k.bounds().first << "," << k.bounds().second << "));";
    return ost;
}


const Eigen::MatrixXd& matrix_cache::get_matrix(double branch_length, double sigma, boundaries bounds) const
{
    // cout << "Matrix request " << size << "," << branch_length << "," << lambda << endl;

    matrix_cache_key key(bounds, sigma, branch_length);
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

void matrix_cache::precalculate_matrices(const std::vector<double>& sigmas, const set<boundaries>& boundses, const std::set<double>& branch_lengths)
{
    // build a list of required matrices
    vector<matrix_cache_key> keys;
    for (auto bounds : boundses)
    {
        for (double branch_length : branch_lengths)
        {
            for (double sigma : sigmas)
            {
                matrix_cache_key key(bounds, sigma, branch_length);
                if (_matrix_cache.find(key) == _matrix_cache.end())
                {
                    keys.push_back(key);
                }
            }
        }
    }

    // calculate matrices in parallel
    size_t i = 0;
    size_t num_keys = keys.size();
    vector<boundaries> vBounds(keys.size());
    vector<double> vBranches(keys.size());
    vector<double> vSigmas(keys.size());
    transform(keys.begin(), keys.end(), vBounds.begin(), [](matrix_cache_key k) { return k.bounds(); });
    transform(keys.begin(), keys.end(), vBranches.begin(), [](matrix_cache_key k) { return k.branch_length(); });
    transform(keys.begin(), keys.end(), vSigmas.begin(), [](matrix_cache_key k) { return k.sigma() * k.sigma() / 2; });
    auto matrices = ConvProp_bounds_batched(vBranches, vSigmas, DiffMat::instance(), vBounds);
    for (i = 0; i < num_keys; ++i)
    {
        _matrix_cache[keys[i]] = matrices[i];
    }

}

void matrix_cache::set_matrix(double branch_length, double sigma, boundaries bounds, const Eigen::MatrixXd& m)
{
    matrix_cache_key key(bounds, sigma, branch_length);
    _matrix_cache[key] = m;
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
    double t = 0.0;
    for (int i = 0; i < 31; i++)
    {
        t += 0.1;
        matrix_cache_key key(boundaries(0, t), 0.01, 0.3);
        keys.insert(key);
    }
    CHECK_EQ(31, keys.size());

    matrix_cache_key key(boundaries(0, 3.0), 0.01, 0.3);
    CHECK_EQ(1, keys.count(key));
}

TEST_CASE("get_matrix returns correct matrix based on key")
{
    matrix_cache cache;
    Matrix3d m1;
    cache.set_matrix(44, 3.2642504711034, boundaries(0, 84.9), m1);

    Matrix2d m2;
    cache.set_matrix(44, 3.1010379475482, boundaries(0, 84.9), m2);

    CHECK_EQ(4, cache.get_matrix(44, 3.1010379475482, boundaries(0, 84.9)).size());
    CHECK_EQ(9, cache.get_matrix(44, 3.2642504711034, boundaries(0, 84.9)).size());
}

TEST_CASE("matrix_cache_key checks sigma for less than")
{
    matrix_cache_key key(boundaries(0, 84.9), 3.2642504711034, 44);
    matrix_cache_key key2(boundaries(0, 84.9), 3.1010379475482, 44);

    CHECK(key2 < key);
}
