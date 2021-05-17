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
#include "lambda.h"

#ifdef HAVE_BLAS
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
    ost << "(" << k.branch_length() << ", (" << k.bounds().first << "," << k.bounds().second << "));";
    return ost;
}

matrix_cache::matrix_cache(const lambda *p_sigma)
{
    _sigma_squared = dynamic_cast<const single_lambda*>(p_sigma)->get_single_lambda();
}

matrix_cache::~matrix_cache()
{

}

const Eigen::MatrixXd& matrix_cache::get_matrix(double branch_length, boundaries bounds) const
{
    // cout << "Matrix request " << size << "," << branch_length << "," << lambda << endl;

    matrix_cache_key key(bounds, branch_length);
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

void matrix_cache::precalculate_matrices(const set<boundaries>& boundses, const std::set<double>& branch_lengths)
{
    // build a list of required matrices
    vector<matrix_cache_key> keys;
    for (auto bounds : boundses)
    {
        for (double branch_length : branch_lengths)
        {
            matrix_cache_key key(bounds, branch_length);
            if (_matrix_cache.find(key) == _matrix_cache.end())
            {
                keys.push_back(key);
            }
        }
    }

    // calculate matrices in parallel
    size_t i = 0;
    size_t num_keys = keys.size();

    for (i = 0; i < num_keys; ++i)
    {
        double branch_length = keys[i].branch_length();
        _matrix_cache[keys[i]] = ConvProp_bounds(branch_length, _sigma_squared * _sigma_squared /2, DiffMat::instance(), keys[i].bounds());
    }

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

