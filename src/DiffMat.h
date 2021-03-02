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
};

Eigen::MatrixXd ConvProp_bounds(double t, double cCoeff, const DiffMat& dMat, std::pair<double, double> bounds);
Eigen::VectorXd VectorPos_bounds(int x, int Npts, std::pair<int, int> bounds);

#endif
