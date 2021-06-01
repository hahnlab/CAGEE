#pragma once

#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <vector>
#include <Eigen/Dense>

class gpu_multiplier
{
    cublasHandle_t handle;
    std::vector<cuDoubleComplex*> d_matrixA, d_matrixB, d_matrixResult;
    cuDoubleComplex** d_A, ** d_B, ** d_C;
    size_t max_matrix_count;
public:
    gpu_multiplier(int Npts);

    std::vector<Eigen::MatrixXd> doit(const std::vector<Eigen::MatrixXcd>& matrices, const Eigen::MatrixXcd& transpose);
};


