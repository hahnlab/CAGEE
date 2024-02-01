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

```
conda create -n cagee -c bioconda cagee
```

If you would like to build from source code, download the latest release from the [**release**](https://github.com/hahnlab/CAGEE/releases/) page, and see the [**manual**](./docs/manual/cagee_manual.md#Installation) for details. CAGEE depends on CMake, Boost, and Eigen for building.

There is also a Docker container available.


## Usage

See the [manual](./docs/manual/cagee_manual.md) for [**quick start**](./docs/manual/cagee_manual.md#Quick-Start) and detailed [**usage instructions**](./docs/manual/cagee_manual.md#Starting-with-CAGEE-in-Inference-Mode).
