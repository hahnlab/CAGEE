#include "DiffMat.h"
#include "easylogging++.h"
#include "doctest.h"

using namespace Eigen;
using namespace std;

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
    MatrixXd result = a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
    VLOG(MATRIX) << "Matrix for cCoeff: " << cCoeff << ", t: " << t << ", Max value: " << bounds.second;
    VLOG(MATRIX) << result;
    VLOG(MATRIX) << "Matrix end";

    return result;
}

VectorXd VectorPos_bounds(int x, int Npts, pair<int, int> bounds) {
    VectorXd X = VectorXd::Zero(Npts);
    if (x == bounds.second)
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
