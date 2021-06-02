#pragma once

#ifdef HAVE_CUDA
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <vector>
#include <Eigen/Dense>

class gpu_multiplier
{
    cublasHandle_t handle;
    size_t max_matrix_count;
    std::vector<cuDoubleComplex*> d_matrixA, d_matrixB, d_matrixResult;
    cuDoubleComplex** d_A, ** d_B, ** d_C;
public:
    gpu_multiplier(int Npts);

    std::vector<Eigen::MatrixXd> doit(const std::vector<Eigen::MatrixXcd>& matrices, const Eigen::MatrixXcd& transpose);
};
#endif


