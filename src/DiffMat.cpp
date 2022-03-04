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

const DiffMat& DiffMat::instance()
{
    if (!p_diffmat)
        p_diffmat = new DiffMat(DISCRETIZATION_RANGE);

    return *p_diffmat;
}

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
        for (int i = 0; i < Npts; ++i)
            for (int j = 0; j < Npts; ++j)
                temp(i, j) = dMat.passage(i, j) * expD[j];
        vTemp[k] = temp;
    }
    MatrixXcd transpose = dMat.passage.transpose();

    return dMat.multiplier->doit(vTemp, transpose);
}

MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, boundaries bounds) {
    return ConvProp_bounds_batched(vector<double>({ t }), vector<double>({ cCoeff }), dMat, vector<boundaries>({ bounds }))[0];
}

VectorXd VectorPos_bounds(double x, int Npts, boundaries bounds) {
    
    if (x < bounds.first || x > bounds.second)
    {
        std::ostringstream ost;
        ost << "VectorPos_bounds error: " << x << " not between " << bounds.first << " and " << bounds.second;
        throw std::runtime_error(ost.str());
    }
    VectorXd X = VectorXd::Zero(Npts);
    if (abs(x - bounds.second) < MATRIX_EPSILON)
    {
        X[Npts - 1] = 1;
    }
    else
    {
        double nx = (Npts - 1) * (x - bounds.first) / double(bounds.second - bounds.first);
        int ix = floor(nx);
        double ux = nx - ix;
        X[ix + 1] = ux;
        X[ix] = 1 - ux;
    }
    return X.unaryExpr([Npts, bounds](double x) {return x * (Npts - 1) / double(bounds.second - bounds.first); });
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


TEST_CASE("VectorPos_bounds handles negative numbers")
{
    VectorPos_bounds(-2.56793, 200, boundaries(-10.0, 20));

}

TEST_CASE("VectorPos_bounds throws on x outside boundaries")
{
    CHECK_THROWS_WITH(VectorPos_bounds(-2, 200, boundaries(0.0, 20)), "VectorPos_bounds error: -2 not between 0 and 20");
}
