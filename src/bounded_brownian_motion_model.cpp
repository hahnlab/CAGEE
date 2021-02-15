#include <vector>
#include <Eigen/Dense>

#include "doctest.h"

#include "bounded_brownian_motion_model.h"
#include "optimizer_scorer.h"
#include "user_data.h"
#include "matrix_cache.h"

using namespace std;
using namespace Eigen;


vector<size_t> build_reference_list(const vector<gene_family>& families)
{
    const size_t invalid_index = std::numeric_limits<size_t>::max();

    vector<size_t> reff;
    const int num_families = families.size();
    reff.resize(num_families, invalid_index);
    for (int i = 0; i < num_families; ++i) {
        if (reff[i] != invalid_index) continue;

        reff[i] = i;

        for (int j = i + 1; j < num_families; ++j) {
            if (reff[j] == invalid_index)
            {
                if (families[i].species_size_match(families[j]))
                {
                    reff[j] = i;
                }
            }
        }
    }

    return reff;
}


optimizer_scorer* bounded_brownian_motion_model::get_sigma_optimizer(const user_data& data)
{
    return new sigma_optimizer_scorer(data.p_tree, data.gene_families);
}


std::vector<double> sigma_optimizer_scorer::initial_guesses()
{
    return vector<double> { 0.0000001 };
}

class DiffMat {
public:
    MatrixXd Diff;
    MatrixXcd passage;

    void Create(int Npts);
};

MatrixXd ConvProp_bounds(double X, double t, double cCoeff, const DiffMat& dMat, pair<double, double> bounds) {
    VectorXcd eig = dMat.Diff.eigenvalues();

    auto tP = dMat.passage.transpose();

    int Npts = dMat.Diff.cols();
    double tau = pow((bounds.second - bounds.first) / (Npts - 1), 2);
    MatrixXd expD(Npts, Npts);
    for (int i = 0; i < Npts; ++i)
    { 
        expD(i, i) = exp(cCoeff * (t / tau) * eig[i].real());
    }
    MatrixXcd a = dMat.passage * expD * tP * X;
    return a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
}

void DiffMat::Create(int Npts) {
    // Npts is the number of points in which the interval is discretized
    MatrixXd A(Npts, Npts);
    A.setZero();
    for (int i = 0; i < Npts-1; ++i) {
        A(i, i) = -2;
        A(i, i + 1) = 1;
        A(i + 1, i) = 1;
    }
    A(0, 0) = -1;
    A(Npts-1, Npts-1) = -1;
    
    Diff = A;
    EigenSolver<MatrixXd> es(Diff);
    passage = es.eigenvectors();
}

double sigma_optimizer_scorer::calculate_score(const double* values)
{
    DiffMat dMat;
    dMat.Create(200);
    double sigma = values[0];

    LOG(INFO) << "Sigma: " << sigma;

    for (auto& family : _families)
        for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
        {
            if ((*it)->is_leaf())
            {
                double t = (*it)->get_branch_length();
                double X = family.get_species_size((*it)->get_taxon_name());
                ConvProp_bounds(X, t, sigma*sigma/2, dMat, pair<double, double>(0,200));
            }
        }

    return 0.1;
}

TEST_CASE("DiffMat creates the expected matrix")
{
    DiffMat dMat;
    dMat.Create(3);

    Matrix3d expected;
    expected << 
        -1, 1, 0,
        1, -2, 1,
        0, 1, -1;

    CHECK(dMat.Diff == expected);
}

TEST_CASE("ConvProp_bounds")
{
    DiffMat dMat;
    dMat.Create(3);
    MatrixXd actual = ConvProp_bounds(1.0, 2.0, 3.0, dMat, pair<double, double>(0.0, 3.0));
    Matrix3d expected;
    expected <<
        0, 0.502323, 0.298648,
        0.502323, 0.333557,  0.16412,
        0.298648,  0.16412,  1.76198;

    for (int i = 0; i<3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK(actual(i,j) == doctest::Approx(expected(i,j)));
}
