#ifndef BLAS_MULTIPLIER_H
#define BLAS_MULTIPLIER_H

#ifdef BLAS_FOUND

#ifdef MKL_FOUND
#include "mkl.h"
#else
#include "cblas.h"
#endif
#endif

#include <Eigen/Dense>
#include <vector>

class blas_multiplier 
{
    size_t max_matrix_count;
public:
    std::vector<Eigen::MatrixXd> doit(const std::vector<Eigen::MatrixXcd>& matrices, const Eigen::MatrixXcd& transpose);
};
#endif

