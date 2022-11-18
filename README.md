# CAGEE

<div>
<h3>
Software for <bold>C</bold>omputational <bold>A</bold>nalysis of <bold>G</bold>ene <bold>E</bold>xpression <bold>E</bold>volution
</h3>
</div>

CAGEE analyzes changes in global or sample- or clade-specific gene expression taking into account phylogenetic history, and provides a statistical foundation for evolutionary inferences.
The program (v1.0) uses Brownian motion to model gene expression changes across a user-specified phylogenetic tree.
The reconstructed distribution of counts and their inferred evolutionary rate $\sigma^2$ generated under this model provides a basis for assessing the significance of the observed differences among taxa.

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

## Installation

See the [**release**](https://github.com/hahnlab/CAGEE/releases/) page for the latest official build and precompiled binaries for your OS.

CAGEE has a few [**dependencies**](./docs/manual/cagee_manual.md#dependencies) which must be installed before building or running, including [Boost](https://www.boost.org/users/download/), [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), and a BLAS numeric library such as [MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)

You can also install CAGEE via `conda`: [ command to be added once uploaded to bioconda ]

See the [**manual**](./docs/manual/cagee_manual.md#Installation) for details, including building from source (not recommended)

## Usage

See the [manual](./docs/manual/cagee_manual.md) for [**quick start**](./docs/manual/cagee_manual.md#Quick-Start) and detailed [**usage instructions**](./docs/manual/cagee_manual.md#Starting-with-CAGEE-in-Inference-Mode).
