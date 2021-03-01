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

#ifdef HAVE_BLAS
#ifdef HAVE_OPENBLAS
#include "cblas.h"
#else
#include "mkl.h"
#endif
#endif

extern std::mt19937 randomizer_engine;
using namespace Eigen;

class DiffMat {
public:
    MatrixXd Diff;
    MatrixXcd passage;
    VectorXcd eig;

    DiffMat(int Npts);
};

MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, pair<double, double> bounds) {
    auto tP = dMat.passage.transpose();

    int Npts = dMat.Diff.cols();
    double tau = pow((bounds.second - bounds.first) / (Npts - 1), 2);
    MatrixXd expD(Npts, Npts);
    expD.setZero();
    for (int i = 0; i < Npts; ++i)
    {
        expD(i, i) = exp(cCoeff * (t / tau) * dMat.eig[i].real());
    }
    MatrixXcd a = dMat.passage * expD * tP;
    return a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
}

DiffMat::DiffMat(int Npts) {
    // Npts is the number of points in which the interval is discretized
    MatrixXd A(Npts, Npts);
    A.setZero();
    for (int i = 0; i < Npts - 1; ++i) {
        A(i, i) = -2;
        A(i, i + 1) = 1;
        A(i + 1, i) = 1;
    }
    A(0, 0) = -1;
    A(Npts - 1, Npts - 1) = -1;

    Diff = A;
    EigenSolver<MatrixXd> es(Diff);
    passage = es.eigenvectors();
    eig = Diff.eigenvalues();
    VLOG(MATRIX) << "Eigenvalues for Diff matrix";
    VLOG(MATRIX) << eig;
    VLOG(MATRIX) << "Eigenvalues end";

}

bool matrix::is_zero() const
{
    return *max_element(values.begin(), values.end()) == 0;
}

//! Take in a matrix and a vector, compute product, return it
/*!
This function returns a likelihood vector by multiplying an initial likelihood vector and a transition probability matrix.
A minimum and maximum on the parent's and child's family sizes is provided. Because the root is forced to be >=1, for example, s_min_family_size for the root could be set to 1.
*/
void matrix::multiply(const vector<double>& v, int s_min_family_size, int s_max_family_size, int c_min_family_size, int c_max_family_size, double* result) const
{
    assert(c_min_family_size < c_max_family_size);
    if (v.size() <= size_t(c_max_family_size - c_min_family_size))
    {
        std::ostringstream ost;
        ost << "Matrix error: size " << v.size() << " less than max family size (" << c_max_family_size << ") - min family size (" << c_min_family_size << ")";
        throw std::runtime_error(ost.str());
    }
#ifdef HAVE_BLAS
    double alpha = 1.0, beta = 0.;
    int m = s_max_family_size - s_min_family_size + 1;
    int k = c_max_family_size - c_min_family_size + 1;
    int n = 1;
    const double *sub = &values[0] + s_min_family_size*_size + c_min_family_size;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, sub, _size, &v[0], n, beta, result, n);
#else
    //cout << "Matrix multiply " << matrix.size() << "x" << v.size() << " (submatrix " << s_min_family_size << ":" << s_max_family_size;
    //cout << " " << c_min_family_size << ":" << c_max_family_size << ")" << endl;

    //assert(s_max_family_size - s_min_family_size == c_max_family_size - c_min_family_size);

    for (int s = s_min_family_size; s <= s_max_family_size; s++) {
        result[s - s_min_family_size] = 0;

        for (int c = c_min_family_size; c <= c_max_family_size; c++) {
            result[s - s_min_family_size] += get(s, c) * v[c - c_min_family_size];
        }
    }
#endif
}

int matrix::select_random_y(int x, int max) const
{
    assert(x < _size);
    assert(max < _size);
    el::Logger* l = el::Loggers::getLogger("default");
    bool enabled = l->typedConfigurations()->enabled(el::Level::Trace);
    if (enabled)
    {
        ostringstream ost;
        ost << "Selecting random value from: ";
        for (int i = 0; i < max; ++i)
        {
            ost << *(values.begin() + x * _size + i) << " ";
        }
        LOG(TRACE) << ost.str();
    }

    std::discrete_distribution<int> distribution(values.begin() + x*_size, values.begin() + x*_size + max);
    return distribution(randomizer_engine);
}

DiffMat* matrix_cache::_p_diffmat = nullptr;

matrix_cache::matrix_cache() : _matrix_size(DISCRETIZATION_RANGE)
{
    if (!_p_diffmat)
        _p_diffmat = new DiffMat(DISCRETIZATION_RANGE);
}

matrix_cache::~matrix_cache()
{
    for (auto m : _matrix_cache)
    {
        delete m.second;
    }
}

const matrix* matrix_cache::get_matrix(double branch_length, double lambda) const {
    // cout << "Matrix request " << size << "," << branch_length << "," << lambda << endl;

    matrix *result = NULL;
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

bool matrix_cache::is_saturated(double branch_length, double lambda)
{
    double alpha = lambda*branch_length / (1 + lambda*branch_length);
    return (1 - 2 * alpha) < 0;
}

MatrixXd matrix_cache::get_matrix(double branch_length, double sigma, double max_value) const
{
    MatrixXd mxd = ConvProp_bounds(branch_length, sigma * sigma / 2, *_p_diffmat, pair<double, double>(0.0, max_value));
    VLOG(MATRIX) << "Matrix for sigma: " << sigma << ", Branch length: " << branch_length << ", Max value: " << max_value;
    VLOG(MATRIX) << mxd;
    VLOG(MATRIX) << "Matrix end";
    return mxd;
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
    vector<matrix*> matrices(keys.size());
    generate(matrices.begin(), matrices.end(), [this] { return new matrix(this->_matrix_size); });

    size_t i = 0;
    size_t num_keys = keys.size();

    for (i = 0; i < num_keys; ++i)
    {
        double sigma = keys[i].lambda();
        double branch_length = keys[i].branch_length();
        MatrixXd mxd = ConvProp_bounds(branch_length, sigma*sigma/2, *_p_diffmat, pair<double, double>(0.0, _matrix_size));

        matrix* m = matrices[i];
        for (int j = 0; j < DISCRETIZATION_RANGE; ++j)
            for (int k = 0; k < DISCRETIZATION_RANGE; ++k)
                m->set(j, k, mxd(j, k));
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

TEST_CASE("DiffMat creates the expected matrix")
{
    DiffMat dMat(3);

    Matrix3d expected;
    expected <<
        -1, 1, 0,
        1, -2, 1,
        0, 1, -1;

    CHECK(dMat.Diff == expected);
}

TEST_CASE("ConvProp_bounds")
{
    DiffMat dMat(3);
    MatrixXd actual = ConvProp_bounds(2.0, 3.0, dMat, pair<double, double>(0.0, 3.0));

    Matrix3d expected;
    expected <<
        0.368131, 0.333222, 0.298648,
        0.333222, 0.333557, 0.333222,
        0.298648, 0.333222, 0.368131;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK(actual(i, j) == doctest::Approx(expected(i, j)));
}
