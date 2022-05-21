#include "doctest.h"
#include "easylogging++.h"

#include "DiffMat.h"

#ifdef HAVE_CUDA
#include "gpu_multiplier.h"
#endif

#ifdef BLAS_FOUND
#include "mkl_multiplier.h"
#endif

using namespace Eigen;
using namespace std;

class eigen_multiplier
{
public:
    vector<MatrixXd> doit(const vector<MatrixXcd>& matrices, const MatrixXcd& transpose);
};

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

#ifdef BLAS_FOUND
    multiplier = new mkl_multiplier();
#elif defined HAVE_CUDA    
    multiplier = new gpu_multiplier(Npts);
#else
    LOG(WARNING) << "No math library available, calculations may be slow";
    multiplier = new eigen_multiplier();
#endif
}

DiffMat* p_diffmat = nullptr;

vector<MatrixXd> ConvProp_bounds_batched(vector<double> vt, vector<double> cCoeff, const DiffMat& dMat, vector<boundaries> vbounds) {
    // Calculate the transition density (dMat.diff to the power of cCoeff * t * (n-1)^2 / (b-a)^2
    // using eigenvectors to speed up the calculation
    size_t count = vt.size();
    vector<MatrixXcd> vTemp(count);
    int Npts = dMat.Diff.cols();

    for (size_t k = 0; k < count; ++k)
    {
        double tau = pow((vbounds[k].second - vbounds[k].first) / (Npts - 1), 2);
        vector<double> expD(Npts);
        for (int i = 0; i < Npts; ++i)
        {
            expD[i] = exp(cCoeff[k] * (vt[k] / tau) * dMat.eig[i].real());
        }
        MatrixXcd temp(Npts, Npts);
#pragma omp parallel for
        for (int i = 0; i < Npts; ++i)
        {
            for (int j = 0; j < Npts; ++j)
                temp(i, j) = dMat.passage(i, j) * expD[j];
        }
        vTemp[k] = temp;
    }
    MatrixXcd transpose = dMat.passage.transpose();

    return dMat.multiplier->doit(vTemp, transpose);
}

MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, boundaries bounds) {
    return ConvProp_bounds_batched(vector<double>({ t }), vector<double>({ cCoeff }), dMat, vector<boundaries>({ bounds }))[0];
}

// This function is used in tight parallel loops so make an effort to avoid having it allocate memory
void VectorPos_bounds(double x, boundaries bounds, VectorXd& result) {
    
    int Npts = result.size();
    if (x < bounds.first || x > bounds.second)
    {
        std::ostringstream ost;
        ost << "VectorPos_bounds error: " << x << " not between " << bounds.first << " and " << bounds.second;
        throw std::runtime_error(ost.str());
    }
    result.setConstant(0);
    if (abs(x - bounds.second) < MATRIX_EPSILON)
    {
        LOG(WARNING) << "Calculating probability for " << x << " with upper bound of " << bounds.second;
        result[Npts - 1] = 1;
    }
    else
    {
        double nx = (Npts - 1) * (x - bounds.first) / double(bounds.second - bounds.first);
        int ix = floor(nx);
        double ux = nx - ix;
        result[ix + 1] = ux;
        result[ix] = 1 - ux;
    }
    transform(result.begin(), result.end(), result.begin(), [Npts, bounds](double x) {return x * (Npts - 1) / double(bounds.second - bounds.first); });
}

vector<MatrixXd> eigen_multiplier::doit(const vector<MatrixXcd>& matrices, const MatrixXcd& transpose)
{
    int Npts = transpose.rows();
    MatrixXcd a(Npts, Npts);

    size_t count = matrices.size();
    vector<MatrixXd> vResult(count);
    for (size_t k = 0; k < count; ++k)
    {
        MatrixXcd a = matrices[k] * transpose;
#ifdef _WIN32
        vResult[k] = MatrixXd(Npts, Npts);
        for (int i = 0; i < Npts; ++i)
            for (int j = 0; j < Npts; ++j)
                vResult[k](i, j) = max(a(i, j).real(), 0.0);
#else
        vResult[k] = a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
#endif
    }
    return vResult;
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
    MatrixXd actual = ConvProp_bounds(2.0, 3.0, dMat, boundaries(0.0, 3.0));

    Matrix3d expected;
    expected <<
        0.368131, 0.333222, 0.298648,
        0.333222, 0.333557, 0.333222,
        0.298648, 0.333222, 0.368131;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK(actual(i, j) == doctest::Approx(expected(i, j)));
}

TEST_CASE("ConvProp_bounds returns different values for different values of sigma")
{
    DiffMat dMat(200);
    MatrixXd actual = ConvProp_bounds(44, 3.2642504711034, dMat, boundaries(0.0, 85.0028));
    MatrixXd a2 = ConvProp_bounds(44, 3.1010379475482, dMat, boundaries(0.0, 85.0028));
    CHECK(actual != a2);
}

TEST_CASE("VectorPos_bounds")
{
    VectorXd actual(20);
    VectorPos_bounds(7.0, pair<double, double>(0, 10), actual);
    vector<double> expected{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.33, 0.57, 0, 0, 0, 0, 0 };
    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(expected[i]), actual[i]);
    }
}

TEST_CASE("VectorPos_bounds at right edge")
{
    VectorXd actual(20);
    VectorPos_bounds(10.0, pair<double, double>(0, 10), actual);
    vector<double> expected{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.9 };
    CHECK_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        CHECK_EQ(doctest::Approx(expected[i]), actual[i]);
    }
}


TEST_CASE("VectorPos_bounds handles negative numbers")
{
    VectorXd v(200);
    VectorPos_bounds(-2.56793, boundaries(-10.0, 20), v);

}

TEST_CASE("VectorPos_bounds throws on x outside boundaries")
{
    VectorXd v(200);
    CHECK_THROWS_WITH(VectorPos_bounds(-2, boundaries(0.0, 20), v), "VectorPos_bounds error: -2 not between 0 and 20");
}
