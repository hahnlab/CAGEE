#include <vector>

#include "bounded_brownian_motion_model.h"
#include "optimizer_scorer.h"
#include "user_data.h"
#include "matrix_cache.h"

using namespace std;

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
    return vector<double> { 0.5};
}

void ConvProp_bounds(double X, double t, double dCoeff, const matrix& dMat, vector<double> bounds) {
#if 0
    vDiag = dMat$diag; P = dMat$passage; tP = t(P);
    Npts = dim(dMat$diag)[1];
    tau = ((bounds[2] - bounds[1]) / (Npts - 1)) ^ 2;
    expD = matrix(0, Npts, Npts);
    for (i in 1 : Npts) { expD[i, i] = exp(exp(dCoeff) * (t / tau) * diag(vDiag)[i]) }
    a = P % *%expD % *%tP % *%X;
    return(apply(a, 1, function(x) max(x, 0)));
#endif
}

void DiffMat(int Npts) {
#if 0
    // Npts is the number of points in which the interval is discretized
    A = matrix(0, Npts, Npts);
    for (i in 1 : (Npts - 1)) {
        A[i, i] = -2;
        A[i, i + 1] = 1;
        A[i + 1, i] = 1;
    }
    A[1, 1] = A[Npts, Npts] = -1;
    double eig = eigen(A)
    return (list(Diff = A, diag = diag(eig$values), passage = eig$vectors))
#endif
}

double sigma_optimizer_scorer::calculate_score(const double* values)
{
    double candidate_sigma = values[0];

    LOG(INFO) << "Sigma: " << candidate_sigma;

    for (auto& family : _families)
        for (auto it = _p_tree->reverse_level_begin(); it != _p_tree->reverse_level_end(); ++it)
        {
            if ((*it)->is_leaf())
            {
                double t = (*it)->get_branch_length();
                double X = family.get_species_size((*it)->get_taxon_name());

                ConvProp_bounds(X, t, 1.0, matrix(10), vector<double>());
            }
        }

    return 0.1;
}