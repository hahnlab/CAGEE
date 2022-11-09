#ifndef BLAS_MULTIPLIER_H
#define BLAS_MULTIPLIER_H


#include <Eigen/Dense>
#include <vector>

class blas_multiplier 
{
    size_t max_matrix_count;
public:
    std::vector<Eigen::MatrixXd> doit(const std::vector<Eigen::MatrixXcd>& matrices, const Eigen::MatrixXcd& transpose);
};
#endif

