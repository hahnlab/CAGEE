#ifndef DIFFMAT_H
#define DIFFMAT_H

#include <Eigen/Dense>
#include <utility>

class DiffMat {
public:
    Eigen::MatrixXd Diff;
    Eigen::MatrixXcd passage;
    Eigen::VectorXcd eig;

    DiffMat(int Npts);

    static const DiffMat& instance();
};

Eigen::MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, std::pair<double, double> bounds);
Eigen::VectorXd VectorPos_bounds(double x, int Npts, std::pair<double, double> bounds);

#endif
