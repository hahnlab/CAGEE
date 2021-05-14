#include "DiffMat.h"

#include "doctest.h"
#include "easylogging++.h"
#include "mkl.h"
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

using namespace Eigen;
using namespace std;

cublasHandle_t handle;
cuDoubleComplex* d_matrix, * d_matrixT, * d_matrixResult;

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

#ifdef HAVE_CUDA
    cudaMalloc((void**)&d_matrixT, Npts * Npts * sizeof(complex<double>));
    cudaMalloc((void**)&d_matrix, Npts * Npts * sizeof(complex<double>));
    cudaMalloc((void**)&d_matrixResult, Npts * Npts * sizeof(complex<double>));
#endif
}

DiffMat* p_diffmat = nullptr;

const DiffMat& DiffMat::instance()
{
    if (!p_diffmat)
        p_diffmat = new DiffMat(DISCRETIZATION_RANGE);

    return *p_diffmat;
}


MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, boundaries bounds) {
    // Calculate the transition density (dMat.diff to the power of cCoeff * t * (n-1)^2 / (b-a)^2
    // using eigenvectors to speed up the calculation
    int Npts = dMat.Diff.cols();
    double tau = pow((bounds.second - bounds.first) / (Npts - 1), 2);
    vector<double> expD(Npts);
    for (int i = 0; i < Npts; ++i)
    {
        expD[i] = exp(cCoeff * (t / tau) * dMat.eig[i].real());
    }
    MatrixXcd temp(Npts, Npts);
    for (int i = 0; i < Npts; ++i)
        for (int j = 0; j < Npts; ++j)
            temp(i, j) = dMat.passage(i, j) * expD[j];

    MatrixXcd transpose = dMat.passage.transpose();
#ifdef HAVE_CUDA
    const complex<double> alpha = 1.0;
    const complex<double> beta = 0.0;
    cudaMemcpy(d_matrix, temp.data(), Npts * Npts * sizeof(complex<double>), cudaMemcpyHostToDevice);
    cudaMemcpy(d_matrixT, transpose.data(), Npts * Npts * sizeof(complex<double>), cudaMemcpyHostToDevice);
    //auto status = cublasZgemmBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, Npts, Npts, Npts, reinterpret_cast<const cuDoubleComplex*>(&alpha),
    auto status = cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, Npts, Npts, Npts, reinterpret_cast<const cuDoubleComplex*>(&alpha),
        d_matrix, Npts,
        d_matrixT, Npts,
        reinterpret_cast<const cuDoubleComplex*>(&beta), d_matrixResult, Npts);
    MatrixXcd a(Npts, Npts);
    cudaMemcpy(a.data(), d_matrixResult, Npts * Npts * sizeof(complex<double>), cudaMemcpyDeviceToHost);
#elif defined HAVE_MKL
    const complex<double> alpha = 1.0;
    const complex<double> beta = 0.0;
    MatrixXcd a(Npts, Npts);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Npts, Npts, Npts, &alpha,
        temp.data(), Npts, transpose.data(), Npts, &beta, a.data(), Npts);
#else
    MatrixXcd a = temp * dMat.passage.transpose();
#endif

#ifdef _WIN32
    MatrixXd result(Npts, Npts);
    for (int i = 0; i < Npts; ++i)
        for (int j = 0; j < Npts; ++j)
            result(i, j) = max(a(i, j).real(), 0.0);
#else
    MatrixXd result = a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
#endif
    VLOG(MATRIX) << "Matrix for cCoeff: " << cCoeff << ", t: " << t << ", Max value: " << bounds.second;
    VLOG(MATRIX) << result;
    VLOG(MATRIX) << "Matrix end";

    return result;
}

VectorXd VectorPos_bounds(double x, int Npts, boundaries bounds) {
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

