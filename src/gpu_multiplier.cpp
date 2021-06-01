#include "gpu_multiplier.h"

#include <complex>
#include <vector>

using namespace std;
using namespace Eigen;

gpu_multiplier::gpu_multiplier(int Npts) : max_matrix_count(500), d_matrixA(500), d_matrixB(500), d_matrixResult(500)
{
    cublasCreate(&handle);

    for (size_t i = 0; i < max_matrix_count; ++i)
    {
        cudaMalloc((void**)&d_matrixA[i], Npts * Npts * sizeof(complex<double>));
        cudaMalloc((void**)&d_matrixB[i], Npts * Npts * sizeof(complex<double>));
        cudaMalloc((void**)&d_matrixResult[i], Npts * Npts * sizeof(complex<double>));
    }

    cudaMalloc((void**)&d_A, max_matrix_count * sizeof(cuDoubleComplex*));
    cudaMalloc((void**)&d_B, max_matrix_count * sizeof(cuDoubleComplex*));
    cudaMalloc((void**)&d_C, max_matrix_count * sizeof(cuDoubleComplex*));
}

vector<MatrixXd> gpu_multiplier::doit(const vector<MatrixXcd>& matrices, const MatrixXcd& transpose)
{
    size_t count = matrices.size();
    if (count > max_matrix_count)
        throw std::runtime_error("Too many matrices!");

    int Npts = transpose.rows();

    const complex<double> alpha = 1.0;
    const complex<double> beta = 0.0;
    for (size_t i = 0; i < count; ++i)
    {
        cudaMemcpy(d_matrixA[i], matrices[i].data(), Npts * Npts * sizeof(complex<double>), cudaMemcpyHostToDevice);
        cudaMemcpy(d_matrixB[i], transpose.data(), Npts * Npts * sizeof(complex<double>), cudaMemcpyHostToDevice);
    }
    cudaMemcpy(d_A, d_matrixA.data(), count * sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, d_matrixB.data(), count * sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, d_matrixResult.data(), count * sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);
#if 1
    auto status = cublasZgemmBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, Npts, Npts, Npts, reinterpret_cast<const cuDoubleComplex*>(&alpha),
        d_A, Npts,
        d_B, Npts,
        reinterpret_cast<const cuDoubleComplex*>(&beta), d_C, Npts, matrices.size());
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("CUDA failure");
#else
    for (size_t i = 0; i < count; ++i)
    {
        auto status = cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, Npts, Npts, Npts, reinterpret_cast<const cuDoubleComplex*>(&alpha),
            d_matrix[i], Npts,
            d_matrixT[i], Npts,
            reinterpret_cast<const cuDoubleComplex*>(&beta), d_matrixResult[i], Npts);
        if (status != CUBLAS_STATUS_SUCCESS)
            throw std::runtime_error("CUDA failure");
    }

#endif

    vector<MatrixXd> vResult(count);
    for (size_t k = 0; k < count; ++k)
    {
        MatrixXcd a(Npts, Npts);
        vResult[k] = MatrixXd(Npts, Npts);
        cudaMemcpy(a.data(), d_matrixResult[k], Npts * Npts * sizeof(complex<double>), cudaMemcpyDeviceToHost);

        // cout << "Complex array:" << a;

#ifdef _WIN32
        for (int i = 0; i < Npts; ++i)
            for (int j = 0; j < Npts; ++j)
                vResult[k](i, j) = max(a(i, j).real(), 0.0);
#else
        vResult[k] = a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
#endif
        //cout << vResult[k];
    }

    return vResult;
}


