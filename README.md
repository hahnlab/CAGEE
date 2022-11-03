# CAGEE

<div>
<h3>
Software for <bold>C</bold>omputational <bold>A</bold>nalysis of <bold>G</bold>ene <bold>E</bold>xpression <bold>E</bold>volution
</h3>
</div>

The purpose of CAGEE is to analyze changes in gene expression in a way that 
accounts for phylogenetic history and provides a statistical foundation for 
evolutionary inferences. The program uses Brownian motion to model gene 
expression levels across a user-specified phylogenetic tree. The distribution of counts (i.e. expression level) 
generated under this model can provide a basis for assessing the significance 
of the observed differences among taxa.

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

Installation
============

### Download

The Github page for CAGEE is https://github.com/hahnlab/CAGEE 

Navigate to a directory that you typically keep source code in and do one of the following:

Download the latest release from the CAGEE release directory https://github.com/hahnlab/CAGEE/releases

If you wish to get the latest version from source, you can run

$ git clone https://github.com/hahnlab/CAGEE.git

Please note that the released versions contain tested and approved code, while the
latest source code may contain experimental and untested features. _It is highly
recommended to use a released version instead._


### Compile

CAGEE requires the Boost and Eigen libraries and a BLAS numeric library such as MKL.

https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html

CAGEE uses the CMake build system. To build, create a "build" directory in the CAGEE
directory and CD to it. Then run "cmake .." followed by "make".

# Running CAGEE 

### Quick Start 

For a typical CAGEE analysis, users are most interested in determining two things:
1) The rate of transcriptome evolution across the tree, among different branches of the tree, or in a specific tissue
2) Which genes have increased or decreased expression, and on which branches of the tree these changes have occurred

This type of analysis requires a minimum of two input files:
1) The **count data file** is a tab-delimited transcript file that contains a column for a description of the gene,
       the unique ID for each transcript, and a column for each taxon that has expression data for each gene.
       (If a functional description is not desired, include this column anyway with a place holder as below (null)).

Example: gene_transcripts.txt
```
Desc	Gene_ID	     human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
ATPase    transcript0     21.5465 1.23451 3.5156  9.09013 14.083  17.4409 19.3906 19.3906
(null)    transcript1     1254.99 123.477 3430.33 375.022 7666.33 185.071 25609.4 686.226
HMG box    transcript2     277.144 1534.61 375.022 339.068 1026.28 91.0759 2537.16 185.071
(null)    transcript3     6.46369 3.99302 167.28  23.9302 0.222631        40.2064 250.549 82.2718
......
....
..
DnaJ    transcript10016     758.884 1876.49 250.549 204.744 1134.89 686.226 686.226 414.778
``` 

2) The **tree file** should contain a binary, rooted, ultrametric, tree in Newick format.  Typically
one obtains this tree using one of several molecular dating methods. If you are unsure if your tree is binary,
rooted, or ultrametric will report this when you try to use it for an analysis. Alternatively, you can use the R package,
Ape with its included functions: is.ultrametric, is.rooted, and is.binary.  

Example: example_tree.txt
```
((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);
```
To get a list of commands just call CAGEE with the -h or --help arguments:

    $ cagee -h

To estimate sigma with no branch rate variation, issue the command:

    $ cagee -i gene_transcripts.txt -t example_tree.txt

To estimate separate sigma values for different lineages in the tree, first identify the branches to which each sigma will apply.
This can be done by making a copy of your tree, and substituting the sigma identifier (1,2,3, etc.) for the branch length values.
For example, to apply a different sigma to the branches leading to human, chimp, and their ancestor, modify the branches as below.

Example sigma tree:
<pre>
((((cat:1,horse:1):1,cow:1):1,(((((<b>chimp:2,human:2):2</b>,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:1,mouse:1):1);
</pre>
For this tree, sigma #2 will be applied to branches leading to human, chimp, and their ancestor while sigma #1 will be applied to all other branches of the tree. 


To run this analysis with both sigmas estimated:


    $ cagee -i gene_transcripts.txt -t example_tree.txt -y example_sigma_tree.txt 

***Caveats***

- **Always** perform multiple runs to ensure convergence, especially if multiple sigmas are used.
- We recommend using the -o flag to assign a unique name to the output directory for each run so that results from previous runs are not overwritten.

----
### Slow Start


CAGEE performs three different operations on either one or two
models. The operations are

-   Estimate Sigma - the traditional function of CAGEE. Takes a tree and
    a file of gene expression values, and performs a maximum likelihood
    calculation to estimate the most likely rate of change across the
    entire tree.

-   Simulate - Given specified values, generate an artificial list of
    gene expressions that matches the values. To generate a simulation,
    pass the --simulate or -s parameter. Either pass a count of transcripts
    to be simulated with the parameter (--simulate=1000) or pass a
    --rootdist (-f) parameter with a file containing the distribution to
    match (see \[rootdist\] for the file format).

All CAGEE options are given at once on the command line. Here is an example:

    cagee -i gene_transcripts.txt -t example_tree.txt

In this example, the -t parameter specifies a file containing the tree
that CAGEE uses; and the -i parameter specifies a list of gene expression values.

Parameters
----------

- 	**--cores, -c**
	
	Number of processing cores to use, requires an integer argument. Default=All available 
	cores.

-   **--Expansion, -E**

    Expansion parameter for Nelder-Mead optimizer, Default=2.

-   **--rootdist, -f**

    Path to root distribution file for simulating datasets.

-   **--help, -h**

    Help menu with a list of all commands.

-   **--infile, -i**

    Path to tab delimited gene transcript file to be analyzed - Required for estimation.

-   **--Iterations, -I**

    Maximum number of iterations that will be performed in sigma search. 
    Default=300 (increase this number if likelihood is still improving when limit is hit).

-   **--fixed\_sigma, -l**

    Value (between 0 and 1) for a single user provided sigma value, otherwise sigma is estimated.

-   **--log\_config, -L**

    Turn on logging, provide name of the configuration file for logging (see example log.config file).

-   **--fixed\_multiple\_sigmas, -m**

    Multiple sigma values, comma separated, must be used in conjunction with sigma tree (-y).

-   **--output\_prefix, -o**

    Output directory - Name of directory automatically created for output. Default=results.

-   **--Reflection, -R**

    Reflection parameter for Nelder-Mead optimizer, Default=1.

-   **--simulate, -s**

    Simulate gene expression values. Either provide an argument of the number of transcripts
	to simulate (-s100, or --simulate=100).

-   **--tree, -t**

    Path to file containing newick formatted tree - Required for estimation.

-   **--sigma\_tree, -y**

    Path to sigma tree, for use with multiple sigmas.
    
-   **--zero\_root, -z**

    Include gene expression values that don't exist at the root, not recommended.


Input files
-----------

- 
- Tree files

    A tree file is specified in Newick format.

        ((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);


    An example may be found in the examples/mammals\_tree.txt file.

-   Transcript files

    Transcript files can be specified in the CAFE input format:

	```
    Desc	Gene_ID	     human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
    ATPase    transcript0     21.5465 1.23451 3.5156  9.09013 14.083  17.4409 19.3906 19.3906
    (null)    transcript1     1254.99 123.477 3430.33 375.022 7666.33 185.071 25609.4 686.226
    HMG box    transcript2     277.144 1534.61 375.022 339.068 1026.28 91.0759 2537.16 185.071
    (null)    transcript3     6.46369 3.99302 167.28  23.9302 0.222631        40.2064 250.549 82.2718
    ......
    ....
    ..
    DnaJ    transcript10016     758.884 1876.49 250.549 204.744 1134.89 686.226 686.226 414.778
	``` 

    The file is tab-separated with a header line giving the order of the
    species. Each line thereafter consists of a description, a gene
    ID, and expression values for each species in the tree.


Examples
========

Sigma Search
-------------

Search for a single sigma value using the example phylogeny and gene expression set in the Examples directory:

    cagee -t example_tree.txt -i gene_transcripts.txt -o singlesigma

Sigma Search with Multiple Sigmas
-----------------------------------

Search for separate sigma values for the chimp/human clade and the rest of the tree separately, using the example phylogeny and gene expression set in the Examples directory:

    cagee -t example_tree.txt -i gene_transcripts.txt -y example_sigma_tree.txt -o doublesigma

Troubleshooting
===============

Logging
-------
More verbose logging can be provided with an EasyLogging log config file.
The file may look like this:

	* GLOBAL:
	   FORMAT               =  "%datetime %msg"
	   FILENAME             =  "cagee.log"
	   ENABLED              =  true
	   TO_FILE              =  true
	   TO_STANDARD_OUTPUT   =  true
	   SUBSECOND_PRECISION  =  6
	   PERFORMANCE_TRACKING =  true
	   MAX_LOG_FILE_SIZE    =  2097152 ## 2MB - Comment starts with two hashes (##)
	   LOG_FLUSH_THRESHOLD  =  100 ## Flush after every 100 logs
	* DEBUG:
	   FORMAT               = "%datetime{%d/%M} %func %msg"

For more information, see https://github.com/amrayn/easyloggingpp#using-configuration-file

Pass the config file to CAGEE with the --log_config flag. For example,

    cagee -c cagee_estimate.cfg --log_config log.config	   

Technical
=========

How does the optimizer work?
----------------------------

The Nelder-Mead optimization algorithm is used. It runs until it can
find a difference of less than 1e-6 in either the calculated score or
the calculated value, or for 10,000 iterations. The parameters that are
used for the optimizer are as follows:

-   rho: 1 (reflection)

-   chi: 2 (expansion)

-   psi: 0.5 (contraction)

-   sigma: 0.5 (shrink)

In some cases, the optimizer suggests values that cannot be calculated
(due to saturation, negative values, or other reasons) In this case, an
infinite score is returned and the optimizer continues.

Acknowledgements
================

Many people have contributed to the CAGEE project, either by code or
ideas. Thanks to:

-   Jason Bertram

-   Ben Fulton

-   Matthew Hahn

-   Mira Han

-	Fabio Mendes

-   Gregg Thomas

-   Jay Tourigny

-	Dan Vanderpool

CAGEE uses the EasyLogging logging framework. https://github.com/amrayn/easyloggingpp

CAGEE uses the DocTest testing framework. https://github.com/onqtam/doctest
