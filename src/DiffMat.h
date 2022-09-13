#ifndef DIFFMAT_H
#define DIFFMAT_H

#include <Eigen/Dense>
#include <utility>
#include <vector>

class eigen_multiplier;
class gpu_multiplier;
class blas_multiplier;

using boundaries = std::pair<double, double>;


class DiffMat {
public:
    Eigen::MatrixXd Diff;
    Eigen::MatrixXcd passage;
    Eigen::VectorXcd eig;

    DiffMat(int Npts);
    void create_or_read_eigenvectors();
    void create_eigenvectors();

#ifdef HAVE_CUDA
    gpu_multiplier* multiplier;
#elif defined BLAS_FOUND
    blas_multiplier* multiplier;
#else
    eigen_multiplier* multiplier;
#endif
};

Eigen::MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, boundaries bounds);


void VectorPos_bounds(double x, boundaries bounds, Eigen::VectorXd& result);

std::vector<Eigen::MatrixXd> ConvProp_bounds_batched(std::vector<double> vt, std::vector<double> cCoeff, const DiffMat& dMat, std::vector<boundaries> vbounds);

#endif
