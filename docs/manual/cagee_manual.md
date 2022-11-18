# CAGEE v1.0 Manual

Beta Release for Review, Nov 2022

## (Brief) Description

CAGEE is software for inferring rates of gene expression evolution (or any other continuous trait) across a phylogenetic tree using a bounded Brownian motion model. In addition to providing estimates for one or more evolutionary rates, $\sigma^2$, it provides ancestral state reconstruction and gene-wise  changes between nodes.

See our publication [ add link to citation here ] for more details.

## Quick Start

1.  Download the appropriate pre-compiled CAGEE binary from the [release page](https://github.com/hahnlab/CAGEE/releases) and save it to your preferred installation directory.  Alternatively, you can install it with dependencies using the [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) environment and package manager (`conda install -n your_env -c bioconda cagee`)
2.  From the command line, make sure the CAGEE installation directory is in your `PATH`  or go to your installation directory and:
3.  Run `cagee -t [path/to/your/newick_tree.nwk] -i [path/to/your/gene_expression_data.tsv]`. The two required inputs are an ultrametric species tree and a file containing the data. Further details on what is required for these are given below.
4.  Review the output printed to the terminal and saved to files in `./results/` for this run with a single $\sigma^2$. Details on how to run models with multiple $\sigma^2$ values are described below.

## Installation

Installation of CAGEE is possible through several avenues, including precompiled binaries, `conda` environments, and building from source.

### Precompiled binaries

Executable files for several common systems are available on the [release page](https://github.com/hahnlab/CAGEE/releases), which will be updated periodically with official releases and patches.  Simply select the appropriate build for your operating system, download, and unpack it to your installation directory.

### `conda` / Anaconda

Using the `conda` environment and package manager included as part of the Anaconda software suite (and available as a standalone manager [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)), you can create your own CAGEE-specific software environment which automagically downloads and installs required dependencies.  To do so, simply run `conda create -n cagee -c bioconda cagee` [ NOTE make sure this actually works when Ben uploads it to bioconda] to install CAGEE from the [bioconda](https://bioconda.github.io/) channel into its own environment, then run `conda activate cagee` before running `cagee` itself.

### Building from source

If you want to build CAGEE from the official release source, or the very latest (and unstable) version from the main GitHub repository, you can do so with `CMake` and the other CAGEE dependencies.

#### Dependencies

CAGEE depends on a number of common libraries which must be installed before attempting to build from source (links current as of October 2022):
  - [CMake](https://cmake.org/download/)
  - [Boost](https://www.boost.org/users/download/)
  - [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
  - BLAS numeric library such as [MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)

#### Compilation

1.  Download your desired source version [ add link to source releases ] or get the latest via `git clone https://github.com/hahnlab/CAGEE.git` (**NOTE** that the released versions contain tested and approved code, while the latest source code may contain experimental and untested features. *It is highly recommended to use a released version instead*)
2. Go to the CAGEE source directory (`cd CAGEE/`), and build CAGEE with the following:
	```
	$ mkdir build
	
	$ cd build/
	
	$ cmake ..
	
	$ make
	```
4. You can then run `make install` to install the binary to the default location (or specify the installation prefix at the CMake step with `-DCMAKE_INSTALL_PREFIX=your/install/prefix/`).

Note that you can speed up the compile step with `make -j [# threads]` if you have multiple cores available.

## Starting with CAGEE in Inference Mode

CAGEE's primary application is to infer one or more evolutionary rates, as represented by $\sigma^2$, across your tree of interest, along with ancestral states for gene expression values at each node of the tree. To do so, it requires at least two valid input files--a tree and gene expression data--as described below.

#### Input File Formats

1. The **tree** file should contain a rooted, ultrametric tree in Newick format, e.g. `((((sp1:2,sp2:2):2,(((sp3:1,sp4:1):1,sp5:2):1,sp6:3):1):1,sp7:5):1,sp8:6)`. We recommend expressing branchlengths in millions of years, if possible.
2. The **gene expression** file contains tab-delimited columns containing a description, unique ID, and expression levels within each species, with one gene on each line (after the header line), e.g.:
	
		DESC    Gene_ID            sp1      sp2      [sp3...]
		SIG1    transcript0    1496.51    14368.5
		SIG2    transcript1    2950.15    3430.33
		SIG3    transcript2    414.778    703.712
		SIG4    transcript3    1193.43    758.884
		SIG5    transcript4    1193.43    111.575
		...

The first two columns must be present in the data file, even if you have no descriptions or unique gene IDs; CAGEE expects the data to start in the third column. The headers of the data columns must match the species names used in the tree file. The data file can contain from one to tens of thousands of genes.

### Estimating a single $\sigma^2$

The default inference mode will determine a single maximum likelihood estimate (MLE) of $\sigma^2$ that best fits your tree.  To run a basic CAGEE inference, you must minimally specify a tree in Newick format (`--tree / -t`) and an input file (`--infile / -i`) of tab-delimited gene expression values in each species as described above.

A minimal command of:

    cagee --tree path/to/your_tree.nwk --infile path/to/your_gene_values.tsv
    
... will estimate a single $\sigma^2$ across your tree as well as ancestral states and changes along each branch. This minimal input will use one thread, the default gamma root distribution paramters ($k = 0.375$, $\theta$ = 1600), and write the inference output files to `results/` in the current working directory.

Optionally, you can also specify the number of cores / threads to use (`--cores / -@`) for increased performance, the parameters of the expected gamma root distribution (`--prior gamma:[k]:[theta]`), and the output directory (`--output_prefix / -o`):

    cagee --cores 4 --tree path/to/your_tree.nwk --infile path/to/your_gene_values.tsv \
        --prior gamma:0.5:1200.0 --output_prefix path/to/output-dir/

You can also provide CAGEE with a static $\sigma^2$ value via `--fixed_sigma / -l`; CAGEE will not estimate $\sigma^2$ in this case, but will instead provide a likelihood value, ancestral states, and changes on each branch using this specific value of $\sigma^2$. (`--fixed_sigma` in this case has a different meaning than its required use in simulation mode [see below])

**NOTE** that `--prior` here in inference mode and `--rootdist` in simulation mode are mutually exclusive: specifying both or the wrong one will return an error.  You also can't mix modes, e.g. run `--simulate` while providing an `--infile`.

### Multiple $\sigma^2$: 1. Multiple branch estimates

CAGEE can estimate multiple $\sigma^2$ values across different branches of the tree.  To do so, you must provide the `--tree` and `--infile` as above, and additionally include `--sigma_tree / -y`. This Newick-formatted tree assigns different $\sigma^2$ parameters to individual branches, in whichever way you want to arrange them. For instance, to specify a model with two $\sigma^2$ parameters, we replace the branch lengths with labels corresponding to your desired rate groups (in this case "1" and "2"), e.g. `((((sp1:1,sp2:1):1,(((sp3:2,sp4:2):2,sp5:2):2,sp6:2):2):1,sp7:1):1,sp8:1`

    cagee --cores [your # threads] --tree path/to/your_tree.nwk --infile path/to/your_gene_values.tsv \
        --sigma_tree path/to/your_sigma_tree.nwk --output_prefix path/to/output-dir/

### Multiple $\sigma^2$: 2. Multiple sample types, tissues, etc.

CAGEE estimates multiple $\sigma^2$ values across various sample types input at each tip; these might be data from different tissues, conditions, timepoints or any other subsample of interest that applies to the given tree.  To use CAGEE in this mode:

1. Include an additional third column labeled `SAMPLETYPE` in your gene expression `--infile`, e.g.

		DESC	Gene_ID		SAMPLETYPE	sp1		sp2      	[sp3...]
		SIG1	transcript0	brain		1496.51    	14368.5
		SIG2	transcript1	brain		2950.15    	3430.33
		SIG3	transcript2	heart		414.778    	703.712
		SIG4	transcript3	heart		1193.43    	758.884
		SIG5	transcript4	lung		1193.43    	111.575
		...

Note that, unlike the `DESC` and `Gene_ID`, the column specifying sample membership must be called `SAMPLETYPE`; it must occur in one of the first three columns. 

2. Specify the sample groups on which to estimate $\sigma^2$ via `--sample_group [SAMPLETYPE_label]` on the command-line.  Each sample group for which you want to estimate a rate must be given a separate argument (`--sample_group brain --sample_group heart`), and groups can be combined for a single rate with a comma-separated list (`--sample_group heart,lung`).

    ```
    cagee --cores [your # threads] --tree path/to/your_tree.nwk --infile path/to/your_gene_values_with_SAMPLETYPE.tsv \
        --sample_group sample_label_1 --sample_group sample_label_2,sample_label_3 --output_prefix path/to/output-dir/
    ```

Given the above input, CAGEE will estimate two values of $\sigma^2$: one for all rows with SAMPLETYPE `sample_label_1` and one for all rows with SAMPLETYPE either `sample_label_2` or `sample_label_3` (since these are both assigned to the same sample group on the command-line). If one instead wanted to estimate a separate value of $\sigma^2$ for sample_label_3, you would simply add a third argument (`--sample_group sample_label_3`).

### Output of a typical run

CAGEE reports its estimates in the terminal and also writes several output files to the specified directory.  These include:
  - `results.txt`, a brief summary of final $\sigma^2$ estimates and likelihoods. This file also contains a tree that labels all ancestral nodes.
  - `ancestral_states.tab`, a tab-separated file of all reconstructed ancestral states.
  - `ancestral_states.tre`, a Nexus-formatted file of all reconstructed ancestral states. 
  - `change.tab`, a tab-separated file of gene expression changes per branch of the tree. "Credible" changes are marked with an `*`.
  - `clade_results.tab`, a summary of the number of genes that have increased or decreased on each branch of the tree. Note: this file only reports credible changes. To see all changes, users should use the flag `--count_all_changes`.
  - `credible_intervals.tab`, a tab-separated file of the inferred 95% credible intervals for ancestral states (and tip values).

 
## Other CAGEE features

CAGEE has a number of additional command line options which you can review with `cagee --help`

### Simulation Mode

To run CAGEE in simulation mode, specify `--simulate` / `-s` with the number of transcripts/genes you want the simulate evolving across the tree. You must also specify your tree (`--tree / -t`) and the overall $\sigma^{2}$ (`--fixed_sigma / -l`):

    cagee --simulate 1000 --tree path/to/your_tree.nwk --fixed_sigma 10

... runs a single simulation across your tree producing 1000 transcripts with data at the tips and stored ancestral states. This call uses a $\sigma^{2}$ value of 10, as well as the default gamma root distribution parameters of $k = 0.375$ and $\theta = 1600.0$, outputting the simulated values to `results/` in the current working directory. The "truth" values are the true values of the ancestral nodes produced during simulation.

Optionally, you can also specify the number of cores / threads to use (`--cores / -@`), the type and/or parameters of the root distribution, (e.g. `--rootdist gamma:k:theta`), and the path of output directory (`--output_prefix / -o`):

    cagee --cores 4 --simulate 1000 --fixed_sigma 10 \
        --tree path/to/your_tree.nwk --rootdist gamma:0.5:1200.0 \
        --output_prefix path/to/output-dir

... runs a simulation of 1000 transcripts on your tree using 4 CPU cores, $\sigma^{2} = 10$, a gamma root distribution with $k = 0.5$ and $\theta = 1200.0$, and is written to your output directory.

Other valid arguments for `--rootdist` include a fixed value at the root for all of your transcripts, or a tab-separated file specifying transcript names (col 1) and root expression values (col 2), useful for starting a simulation from a defined distribution of root values (for instance, those inferred from a different dataset).

NOTE that `--rootdist` in simulation mode and `--prior` in inference mode are mutually exclusive: specifying both or the wrong one will return an error.

To simulate multiple rates across your tree, you must include the `--sigma_tree` option as described above and a comma-separated list of rate values to `--fixed_multiple_sigmas / -m`, e.g. `--fixed_multiple_sigmas 1,3` to simulate rates of $\sigma^{2} = 1$ and $\sigma^{2} = 3$ to branches labeled 1 and 2, respectively.
