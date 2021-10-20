#pragma once

#ifdef BLAS_FOUND
#include "mkl.h"
#include <Eigen/Dense>
#include <vector>

class mkl_multiplier 
{
    size_t max_matrix_count;
public:
    std::vector<Eigen::MatrixXd> doit(const std::vector<Eigen::MatrixXcd>& matrices, const Eigen::MatrixXcd& transpose);
};
#endif

