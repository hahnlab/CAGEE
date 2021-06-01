#ifndef DIFFMAT_H
#define DIFFMAT_H

#include <Eigen/Dense>
#include <utility>
#include <vector>

class gpu_multiplier;

using boundaries = std::pair<double, double>;


class DiffMat {
public:
    Eigen::MatrixXd Diff;
    Eigen::MatrixXcd passage;
    Eigen::VectorXcd eig;

    DiffMat(int Npts);

    static const DiffMat& instance();

#ifdef HAVE_CUDA
    gpu_multiplier* multiplier;
#endif
};

Eigen::MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, boundaries bounds);
Eigen::VectorXd VectorPos_bounds(double x, int Npts, boundaries bounds);

std::vector<Eigen::MatrixXd> ConvProp_bounds_batched(std::vector<double> vt, double cCoeff, const DiffMat& dMat, std::vector<boundaries> vbounds);

#endif
