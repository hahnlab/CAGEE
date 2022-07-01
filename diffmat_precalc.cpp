#include <iostream>
#include "zlib.h"

#include "src/DiffMat.h"

#include "src/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

template<class Derived>
void write_binary(const std::string& filename,
    const Eigen::PlainObjectBase<Derived>& matrix)
{
    typedef typename Derived::Index Index;
    typedef typename Derived::Scalar Scalar;

    gzFile out = gzopen(filename.c_str(), "wb");
    Index rows = matrix.rows(), cols = matrix.cols();

    gzwrite(out, (char*)(&rows), sizeof(Index));
    gzwrite(out, (char*)(&cols), sizeof(Index));
    gzwrite(out, (char*)matrix.data(), rows * cols * sizeof(Scalar));
    gzclose(out);
}

int main(int argc, char *const argv[]) {
     if (argc != 2)
     {
        std::cout << "USAGE: diffmat_precalc N\n";
        std::cout << "N is the size of the diffusion matrix you intend to use\n";
        return 0;
     }
     DiffMat diffusion_matrix(std::stoi(argv[1]));
     diffusion_matrix.create_eigenvectors();
     write_binary(std::string("eigenvalues") + argv[1] + ".bin", diffusion_matrix.eig);
     write_binary(std::string("eigenvectors") + argv[1] + ".bin", diffusion_matrix.passage);
     std::cout << "Eigenvectors written\n";
}


