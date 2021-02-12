#include <vector>
#include <Eigen/Dense>

#include "bounded_brownian_motion_model.h"
#include "optimizer_scorer.h"
#include "user_data.h"
#include "matrix_cache.h"

using namespace std;
using namespace Eigen;

//! @brief Scorer that optimizes for sigma
//! \ingroup optimizer
class sigma_optimizer_scorer : public optimizer_scorer
{
    const clade* _p_tree;
    const vector<gene_family>& _families;

public:
    sigma_optimizer_scorer(const clade *p_tree, const vector<gene_family>& families) : _p_tree(p_tree), _families(families)
    {
    }

    virtual std::vector<double> initial_guesses() override;

    virtual double calculate_score(const double* values) override;
};

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
    Matrix<double,Dynamic,Dynamic> Diff;
    DiagonalMatrix<double, Dynamic> diag;
    Matrix<double, Dynamic, Dynamic> passage;

    void Create(int Npts);
};

void ConvProp_bounds(double X, double t, double cCoeff, const DiffMat& dMat, pair<double, double> bounds) {
    auto vDiag = dMat.diag;
    auto P = dMat.passage; 
    auto tP = P.transpose();
    int Npts = dMat.diag.size();
    double tau = pow((bounds.second - bounds.first) / (Npts - 1), 2);
    Matrix<double, Dynamic, Dynamic> expD(Npts, Npts);
    for (int i = 1; i < Npts; ++i)
    { 
        expD(i, i) = exp(cCoeff * (t / tau) * vDiag.diagonal()[i]);
    }
    auto a = P * expD * tP * X;
    a.unaryExpr([](double x) {return max(x, 0.0); });
}

void DiffMat::Create(int Npts) {
    // Npts is the number of points in which the interval is discretized
    auto A = Matrix<double, Dynamic, Dynamic>(Npts, Npts);
    for (int i = 0; i < Npts; ++i) {
        A(i, i) = -2;
        A(i, i + 1) = 1;
        A(i + 1, i) = 1;
    }
    A(0, 0) = -1;
    A(Npts-1, Npts-1) = -1;
    auto eig = A.eigenvalues();
    Diff = A;
    auto diag = eig.asDiagonal();
    diag.diagonal();
    EigenSolver< Matrix<double, Dynamic, Dynamic>> es(A);
    //passage = es.eigenvectors();
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