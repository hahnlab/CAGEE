#ifndef DIFFMAT_H
#define DIFFMAT_H

#include <Eigen/Dense>
#include <utility>

using boundaries = std::pair<double, double>;


class DiffMat {
public:
    Eigen::MatrixXd Diff;
    Eigen::MatrixXcd passage;
    Eigen::VectorXcd eig;

    DiffMat(int Npts);

    static const DiffMat& instance();
};

Eigen::MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, boundaries bounds);
Eigen::VectorXd VectorPos_bounds(double x, int Npts, boundaries bounds);

#endif
