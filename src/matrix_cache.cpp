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

matrix_cache::matrix_cache() : _matrix_size(DISCRETIZATION_RANGE)
{
}

matrix_cache::~matrix_cache()
{
    for (auto m : _matrix_cache)
    {
        delete m.second;
    }
}

const Eigen::MatrixXd* matrix_cache::get_matrix(double branch_length, double lambda) const {
    // cout << "Matrix request " << size << "," << branch_length << "," << lambda << endl;

    Eigen::MatrixXd*result = NULL;
    matrix_cache_key key(_matrix_size, lambda, branch_length);
    if (_matrix_cache.find(key) != _matrix_cache.end())
    {
        result = _matrix_cache.at(key);
    }

    if (result == NULL)
    {
        ostringstream ost;
        ost << "Failed to find matrix for " << _matrix_size << "," << branch_length << "," << lambda;
        throw std::runtime_error(ost.str());
    }
    return result;
}

void matrix_cache::precalculate_matrices(const std::vector<double>& lambdas, const std::set<double>& branch_lengths)
{
    // build a list of required matrices
    vector<matrix_cache_key> keys;
    for (double lambda : lambdas)
    {
        for (double branch_length : branch_lengths)
        {
            matrix_cache_key key(_matrix_size, lambda, branch_length);
            if (_matrix_cache.find(key) == _matrix_cache.end())
            {
                keys.push_back(key);
            }
        }
    }

    // calculate matrices in parallel
    vector<Eigen::MatrixXd*> matrices(keys.size());
    generate(matrices.begin(), matrices.end(), [this] { return new Eigen::MatrixXd(this->_matrix_size, this->_matrix_size); });

    size_t i = 0;
    size_t num_keys = keys.size();

    for (i = 0; i < num_keys; ++i)
    {
        double sigma = keys[i].lambda();
        double branch_length = keys[i].branch_length();
        MatrixXd mxd = ConvProp_bounds(branch_length, sigma*sigma/2, DiffMat::instance(), pair<double, double>(0.0, _matrix_size));

        Eigen::MatrixXd* m = matrices[i];
        for (int j = 0; j < DISCRETIZATION_RANGE; ++j)
            for (int k = 0; k < DISCRETIZATION_RANGE; ++k)
                (*m)(j, k) = mxd(j, k);
    }

    // copy matrices to our internal map
    for (size_t i = 0; i < keys.size(); ++i)
    {
        _matrix_cache[keys[i]] = matrices[i];
    }
}

std::ostream& operator<<(std::ostream& ost, matrix_cache& c)
{
    ost << c.get_cache_size() << " matrices. Keys: ";
    for (auto& kv : c._matrix_cache)
    {
        ost << "(" << kv.first.branch_length() << "," << kv.first.lambda() << "),";
    }
    return ost;
}

