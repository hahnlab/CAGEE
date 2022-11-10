#include "blas_multiplier.h"

#ifdef BLAS_AVAILABLE
#ifdef MKL_FOUND
#include "mkl.h"
#elif defined OpenBLAS_FOUND
#include "cblas.h"
#elif defined GSL_FOUND
#include "gsl/gsl_blas.h"
#endif

using namespace std;
using namespace Eigen;

vector<MatrixXd> blas_multiplier::doit(const vector<MatrixXcd>& matrices, const MatrixXcd& transpose)
{
    int Npts = transpose.rows();
    const complex<double> alpha = 1.0;
    const complex<double> beta = 0.0;
    MatrixXcd a(Npts, Npts);

    size_t count = matrices.size();
    vector<MatrixXd> vResult(count);
#pragma omp parallel for
    for (size_t k = 0; k < count; ++k)
    {
        MatrixXcd a(Npts, Npts);
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Npts, Npts, Npts, &alpha,
            matrices[k].data(), Npts, transpose.data(), Npts, &beta, a.data(), Npts);

#ifdef _WIN32
        vResult[k] = MatrixXd(Npts, Npts);
        for (int i = 0; i < Npts; ++i)
            for (int j = 0; j < Npts; ++j)
                vResult[k](i, j) = max(a(i, j).real(), 0.0);
#else
        vResult[k] = a.unaryExpr([](complex<double> x) {return max(x.real(), 0.0); });
#endif
    }
    return vResult;
}
#endif
