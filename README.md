# CAGEE

<div>
<h3>
Software for <bold>C</bold>omputational <bold>A</bold>nalysis of <bold>G</bold>ene <bold>E</bold>xpression <bold>E</bold>volution
</h3>
</div>

The purpose of CAGEE is to analyze changes in gene expressions in a way that 
accounts for phylogenetic history and provides a statistical foundation for 
evolutionary inferences. The program uses Brownian motion to model gene 
expression across a user-specified phylogenetic tree. The distribution of expressions 
generated under this model can provide a basis for assessing the significance 
of the observed family size differences among taxa.

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

Installation
============

### Download

The Github page for CAGEE is https://github.com/hahnlab/CAGEE 

Navigate to a directory that you typically keep source code in and do one of the following:

Download the latest release from the CAFE release directory https://github.com/hahnlab/CAGEE/releases

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
1) Which gene families are rapidly evolving 
2) The branches of the tree on which these families are rapidly evolving

This type of analysis requires a minimum of two input files:
1) The **count data file** is a tab-delimited family "counts" file that contains a column for a description of the gene family,
       the unique ID for each family, and a column for each taxon that has expression data for each gene.
       This file is acquired by first peforming a clustering analysis, often using software such as 
       OrthoMCL, SwiftOrtho, FastOrtho, OrthAgogue, or OrthoFinder and then parsing the output into a table
       like the one below (Note: if a functional description is not desired, include this column anyway with a place holder as below (null)).

Example: mammal\_gene\_families.txt
```
Desc	Family ID	human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
ATPase	ORTHOMCL1	 52	 55	 54	 57	 54	  56	  56	 53	 52	57	55	 54
(null)	ORTHOMCL2	 76	 51	 41	 39	 45	  36	  37	 67	 79	37	41	 49
HMG box	ORTHOMCL3	 50	 49	 48	 48	 46	  49	  48	 55	 52	51	47	 55
(null)	ORTHOMCL4	 43	 43	 47	 53	 44	  47	  46	 59	 58	51	50	 55
Dynamin	ORTHOMCL5	 43	 40	 43	 44	 31	  46	  33	 79	 70	43	49	 50
......
....
..
DnaJ	ORTHOMCL10016	 45	 46	 50	 46	 46 	  47	  46	 48	 49	45	44	 48
``` 
2) The **tree file** should contain a binary, rooted, ultrametric, tree in Newick format.  Typically
one obtains this tree using one of several molecular dating methods. If you are unsure if your tree is binary,
rooted, or ultrametric CAFE will report this when you try to use it for an analysis. Alternatively, you can use the R package,
Ape with its included functions: is.ultrametric, is.rooted, and is.binary.  

Example: mammals_tree.txt
```
((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);
```
To get a list of commands just call CAFE with the -h or --help arguments:

    $ cagee -h

To estimate sigma with no among family rate variation issue the command:

    $ cagee -i mammal_gene_families.txt -t mammal_tree.txt

To estimate separate sigma values for different lineages in the tree, first identify the branches to which each sigma will apply.
This can be done by making a copy of your tree, and substituting the sigma identifier (1,2,3, etc.) for the branch length values.
For example, to apply a different sigma to the branches leading to human, chimp, and their ancestor, modify the branches as below.

Example chimphuman_separate_sigma.txt:
<pre>
((((cat:1,horse:1):1,cow:1):1,(((((<b>chimp:2,human:2):2</b>,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:1,mouse:1):1);
</pre>
For this tree, sigma #2 will be applied to branches leading to human, chimp, and their ancestor while sigma #1 will be applied to all other branches of the tree. 


To run this analysis with both sigmas estimated:


    $ cagee -i mammal_gene_families.txt -t mammal_tree.txt -y chimphuman_separate_sigma.txt 

***Caveats***

- **Always** perform multiple runs to ensure convergence, especially if multiple gamma rate categories or sigmas are used.
- More gamma rate categories (-k) does not always mean a better fit to the data. While -k=2 nearly always fits the data better than -k=1, it may be the case that -k=5 has a _worse_ likelihood than -k=3, and convergences between runs is more difficult with more categories. Try several and see what works.
- We recommend using the -o flag to assign a unique name to the output directory for each run so that results from previous runs are not overwritten.
- For all but the simplest of data sets, searching for multiple sigmas with multiple rate categories will result in a failure of convergence to a single optimum between runs. 

### Tutorial

A tutorial is provided in the _tutorial_ directory. It provides 
instructions on how to generate a reasonable gene family groups
in the correct format, dated ultrametric trees, and basic CAFE
analyses. The tutorial contains tutorial.md and some helper
scripts.

----
### Slow Start


CAGEE performs three different operations on either one or two
models. The operations are

-   Estimate Sigma - the traditional function of CAFE. Takes a tree and
    a file of gene family counts, and performs a maximum likelihood
    calculation to estimate the most likely rate of change across the
    entire tree.

-   Simulate - Given specified values, generate an artificial list of
    gene families that matches the values. To generate a simulation,
    pass the --simulate or -s parameter. Either pass a count of families
    to be simulated with the parameter (--simulate=1000) or pass a
    --rootdist (-f) parameter with a file containing the distribution to
    match (see \[rootdist\] for the file format).

The models are
-   Base - Perform computations as if no gamma function is available

-   Gamma - Perform computations as if each gene family can belong to a
    different evolutionary rate category. To use Gamma modelling, pass
    the -k parameter specifying the number of categories to use.

Unlike earlier versions, CAGEE does not require a script. All
options are given at once on the command line. Here is an example:

    cagee -t examples/mammals\_tree.txt -i examples/mammal\_gene\_families.txt -p -k 3

In this example, the -t parameter specifies a file containing the tree
that CAFE uses; and the -i parameter specifies a list of gene families.
The -p, in this instance given without a parameter, indicates that the
root equilibrium frequency will not be a uniform distribution. The -k 
parameter specifies how many gamma rate categories to use.  

Parameters
----------

-   **--fixed\_alpha, -a**

    Alpha value of the discrete gamma distribution to use in category
    calculations. If not specified, the alpha parameter will be
    estimated by maximum likelihood.

- 	**--sigma\_per\_family, -b**
	
	Estimate sigma by family (for testing purposes only).

- 	**--cores, -c**
	
	Number of processing cores to use, requires an integer argument. Default=All available 
	cores.

-   **--error_model, -e**

    Run with no file name to estimate the global error model file. This file can be 
    provided in subsequent runs by providing the path to the Error model file with no 
    spaces (e.g. -eBase\_error\_model.txt).

-   **--Expansion, -E**

    Expansion parameter for Nelder-Mead optimizer, Default=2.

-   **--rootdist, -f**

    Path to root distribution file for simulating datasets.

-   **--help, -h**

    Help menu with a list of all commands.

-   **--infile, -i**

    Path to tab delimited gene families file to be analyzed - Required for estimation.

-   **--Iterations, -I**

    Maximum number of iterations that will be performed in sigma search. 
    Default=300 (increase this number if likelihood is still improving when limit is hit).

-	**--n\_gamma\_cats, -k**

    Number of gamma categories to use. If specified, the Gamma model
    will be used to run calculations; otherwise the Base model will be
    used.

-   **--fixed\_sigma, -l**

    Value (between 0 and 1) for a single user provided sigma value, otherwise sigma is estimated.

-   **--log\_config, -L**

    Turn on logging, provide name of the configuration file for logging (see example log.config file).

-   **--fixed\_multiple\_sigmas, -m**

    Multiple sigma values, comma separated, must be used in conjunction with sigma tree (-y).

-   **--output\_prefix, -o**

    Output directory - Name of directory automatically created for output. Default=results.

-   **--poisson, -p**

    Use a Poisson distribution for the root frequency distribution.
  	If no -p flag is given, a uniform distribution will be used. A value
    can be specified (-p10, or --poisson=10); otherwise the distribution
    will be estimated from the gene families.  

-   **--pvalue, -P**

    P-value to use for determining significance of family size change, Default=0.05.

-   **--chisquare\_compare, -r**

    Chi square compare (not tested).

-   **--Reflection, -R**

    Reflection parameter for Nelder-Mead optimizer, Default=1.

-   **--simulate, -s**

    Simulate families. Either provide an argument of the number of families
	to simulate (-s100, or --simulate=100) or provide a rootdist file giving a set
	of root family sizes to match. Without such a file, the families will be generated
	with root sizes selected randomly between 0 and 100.

-   **--tree, -t**

    Path to file containing newick formatted tree - Required for estimation.

-   **--sigma\_tree, -y**

    Path to sigma tree, for use with multiple sigmas.
    
-   **--zero\_root, -z**

    Include gene families that don't exist at the root, not recommended.


Input files
-----------

- 
- Tree files

    A tree file is specified in Newick format.

        ((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);


    An example may be found in the examples/mammals\_tree.txt file.

-   Transcript files

    Family files can be specified in the CAFE input format:

	```
	Desc	Family ID	human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
	ATPase	ORTHOMCL1	 52	 55	 54	 57	 54	  56	  56	 53	 52	57	55	 54
	(null)	ORTHOMCL2	 76	 51	 41	 39	 45	  36	  37	 67	 79	37	41	 49
	HMG box	ORTHOMCL3	 50	 49	 48	 48	 46	  49	  48	 55	 52	51	47	 55
	(null)	ORTHOMCL4	 43	 43	 47	 53	 44	  47	  46	 59	 58	51	50	 55
	Dynamin	ORTHOMCL5	 43	 40	 43	 44	 31	  46	  33	 79	 70	43	49	 50
	......
	....
	..
	DnaJ	ORTHOMCL10016	 45	 46	 50	 46	 46 	  47	  46	 48	 49	45	44	 48
	``` 

    The file is tab-separated with a header line giving the order of the
    species. Each line thereafter consists of a description, a family
    ID, and counts for each species in the tree.

    Alternatively, the family file can be specified with a set of lines
    beginning with hashtags containing the species order:

        #human
        #chimp
        1       2     ORTHOMCL1
        2       1     ORTHOMCL2
        3       6     ORTHOMCL3
        6       3     ORTHOMCL4

    In this case, the family ID will be in the final column.

-   Root distributions

    A root distribution file takes the format “family\_size
    \[whitespace\] family\_count”, e.g.

        1 1
        2 5
        3 10
        4 15
        5 42

    An example may be found in the
    “examples/poisson\_root\_dist\_1000.txt” file.

-   Error models

    An error model consists of modifications of probabilities of moving
    from one family size to another through the tree. The file is
    structured as a series of lines containing the family size, the
    probability of moving to less than that size, the probability of
    that size staying the same, and the probabilities of the size
    becoming larger. Two header lines must be included: the maximum
    family size to process, and the differential of the probabilities.

        maxcnt:90
        cntdiff -1 0 1
        0 0.00 0.95 0.05
        1 0.05 0.9 0.05
        2 0.05 0.9 0.05
        3 0.05 0.9 0.05
        4 0.05 0.9 0.05

Output
------

All output will be stored to the "results" directory, unless another directory is specified with the "-o" parameter.


-   _model_\_asr.tre

    The file will be named Base\_asr.txt or Gamma\_asr.txt, based on
    which model is in play. It contains the reconstructed states of the
    families, in the Nexus file format
    (https://en.wikipedia.org/wiki/Nexus\_file). A tree is provided for
    each family,with the expected family size set off with an underscore
    from the node ID. 

    In the case of the Gamma reconstruction, the Sigma multipliers for
    each category are given their own section in this file. In this case,
	ony the fastest families are printed.

    IMPORTANT! If you want to view only gene families with significant changes 
    mapped onto the branch on which the change occurred, simple parsing can
    be employed using grep e.g.:

        echo $'#nexus\nbegin trees;'>Significant_trees.tre
        grep "*" _model_\_asr.tre >>Significant_trees.tre
        echo "end;">>Significant_trees.tre

    Open this file in Dendroscope or a similar program and you can view the tree with only families exhibiting a significant change mapped onto the branch on which the change occurred.

-   _model_\_family\_results.txt

    The file will be named Base\_family\_results.txt or
    Gamma\_family\_results.txt, based on which model is in play. It
    consists of a header line giving the name of each node in the tree,
    followed by a line consisting of the family ID, an estimate of
    whether the change is significant (’y’ or ’n’). The characters are 
    separated by tabs.

        #FamilyID   pvalue  Significant at 0.05
        0           0.436   n
        1           0.209   n
        2           0.002   y

    In the Gamma model, an additional set of probabilities are appended,
    representing the likelihood of the family belonging to each gamma
    category.

    IMPORTANT! If you want to know how many and which families underwent a significant expansion/contraction, you can parse this file using simple grep or awk commands. 
            
    To count significant families at the p=0.05 threshold:

            grep -c "\ty" Gammma_family_results.txt 
    To write the just significant families to a file:

            grep "y" Gamma_family_results.txt > Significant_families.txt

    If you used a the default p-value (0.05) in your analysis but have too many significant results, you can filter these to a lower p-value (0.01 in the example) using awk, e.g.:

        awk '$2 < .01 {print $0}' Gamma_family_results.txt > Sig_at_p.01.txt
        
    Count the number of significant families in this file using bash: 

        wc -l Sig_at_p.01.txt


-   _model_\_clade\_results.txt

    The file will be named Base\_clade\_results.txt or
    Gamma\_clade\_results.txt, based on which model is in play. It
    consists of a header line, “Taxon Increase Decrease”, and for each
    node in the tree, a tab-separated count of the number of families
    which have increased and decreased for that node.

-   _model_\_branch\_probabilities.txt

	Contains a tab-separated list of the probabilities calculated for each clade
	and significant family. Probabilities are displayed as "N/A" if the parent
	and child have the same value. In the case of the Gamma model, only
	contains significant families that are calculated to be rapidly changing.
	
-   _model_\_family\_likelihoods.txt

    Using the Base model, a tab-separated file consisting of the header
    line “\#FamilyID Likelihood of Family”, and additional tab-separated
    lines consisting of the family ID and the posterior probability of
    that family.

    Using the Gamma model, the file contains the following values:
    
        #FamilyID	Gamma Cat Mean	Lkhd of Category	Lkhd of Family	Posterior Probability	Significant
        6	        0.0459381	    3.72073e-33	        1.56448e-18	    2.37825e-15	            N/S
        6	        0.333538	    3.15818e-20	        1.56448e-18	    0.0201867	            N/S
        6	        0.969463	    1.37286e-18	        1.56448e-18	    0.877521	            N/S
        6	        2.65106	        1.60035e-19	        1.56448e-18	    0.102293	            N/S
        9	        0.0459381	    3.44912e-22	        2.5524e-16	    1.35133e-06	            N/S
        9	        0.333538	    1.53142e-16	        2.5524e-16	    0.599993	            N/S
        9	        0.969463	    9.99454e-17	        2.5524e-16	    0.391575	            N/S
        9	        2.65106	        2.1518e-18	        2.5524e-16	    0.00843051	            N/S
        ...
        ..
        .
    
    The values for each family are listed on each tab-separated line.

-   _model_\_results.txt

    A file giving the name of the model that was selected (“Base” or
    “Gamma”), the final likelihood of that model, the final value of the
    rate of change of families (Sigma) that was calculated, and, if an
    error model was specified, the final value of that value (Epsilon).

-	_model_\_change.txt - A tab-separated file listing, for each family 
        and clade, the difference between it and its parent clade in the 
		reconstruction that was performed.
		
-	_model_\_count.txt - A tab-separated file listing, for each family 
        and clade, the reconstructed value in that clade.
		
-   simulation.txt

    In the case of simulation, a family file is
    generated with simulated data based on the given input parameters.

-   simulation\_truth.txt
Lambda
    When simulating, an additional file is generated with simulated data
    for internal nodes included. The format is otherwise identical to
    the family file. Rather than node names, each internal node is
    assigned an integer ID.

Examples
========

Sigma Search
-------------

Search for a single sigma value using the mammal phylogeny and gene family set in the Examples directory:

    ../bin/cagee -t mammals_tree.txt -i mammal_gene_families.txt -p -o singlesigma

In earlier versions, the following script would have returned the same
values:

    tree ((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575)
    load -i filtered_cafe_input.txt -t 10 -l reports/run6_caferror_files/cafe_log.txt
    sigma -s -t ((((1,1)1,1)1,(((((1,1)1,1)1,1)1,(1,1)1)1,1)1)1,(1,1)1)

Sigma Search with Multiple Sigmas
-----------------------------------

Search for separate sigma values for the chimp/human clade and the rest of the tree separately, using the mammal phylogeny and gene family set in the Examples directory:

    ../bin/CAGEE -t mammals_tree.txt -i mammal_gene_families.txt -p -y chimphuman_separate_sigma.txt -o doublesigma

In earlier versions, the following script would have returned the same
values:

    ((((cat:68,horse:68):4,cow:73):20,(((((chimp:4,human:4):6,orang:11):2,gibbon:13):7,(macaque:4,baboon:4):16):16,marmoset:36):57):38,(rat:36,mouse:36):96) 
    load -i integral_test_families.txt -t 10
    sigma -s -t ((((1,1)1,1)1,(((((2,2)2,2)2,2)2,(1,1)1)1,1)1)1,(1,1)1)

Error Models
------------

Search for a single sigma value using the Newick tree of mammals in the
Samples folder, and the family files from the CAFE tutorial, applying an
error model:

    CAGEE -t data/mammals_integral_tree.txt -i data/filtered_cafe_input.txt -p -l 0.01 -e data/cafe_errormodel_0.0548828125.txt

In earlier versions, the following script would have returned the same
values:

    tree ((((cat:68,horse:68):4,cow:73):20,(((((chimp:4,human:4):6,orang:11):2,gibbon:13):7,(macaque:4,baboon:4):16):16,marmoset:36):57):38,(rat:36,mouse:36):96) load -i integral_test_families.txt -t 10
    errormodel -all -model cafe_errormodel_0.0548828125.txt
    sigma -l 0.01 -t ((((1,1)1,1)1,(((((1,1)1,1)1,1)1,(1,1)1)1,1)1)1,(1,1)1) -score

Error Model Estimation
----------------------

Estimate an error model:

    CAGEE -t mammals_integral_tree.txt -i filtered_cafe_input.txt -p -e errormodel.txt

In earlier versions, the following script would have returned the same
values:

    load -i filtered_cafe_input.txt -t 4 -l reports/log_run6.txt
    tree ((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575)
    sigma -s -t ((((1,1)1,1)1,(((((1,1)1,1)1,1)1,(1,1)1)1,1)1)1,(1,1)1)
    report reports/report_run6

Troubleshooting
===============

Known Limitations
-----------------

Because the random birth and death process assumes that each family has
at least one gene at the root of the tree, CAGEE will not provide
accurate results if included gene families were not present in the most
recent common ancestor (MRCA) of all taxa in the tree. For example, even
if all taxa have a gene family size of 0, CAFE will assign the MRCA a
gene family of size 1, and include the family in estimation of the birth
and death rate. This difficulty does not affect analyses containing
families that go extinct subsequent to the root node.

If a change in gene family size is very large on a single branch, CAGEE 
may fail to provide accurate λ estimation and/or die during
computation. To see if this is a problem, look at the likelihood scores
computed during the λ search (reported in the log file if the
job finishes). If ALL scores are “-inf” then there is a problem with
large size changes and CAGEE has calculated a probability of 0. Removing the
family with the largest difference in size among species and rerunning
CAGEE should allow λ to be estimated on the remaining data.
If the problem persists, remove the family with the next largest
difference and proceed in a like manner until CAGEE no longer finds
families with zero probability. However, if rapidly evolving families
are removed, care should be taken in interpretation of the estimated
average rate of evolution for the remaining data.

In very large phylogenetic trees there can be many independent sigma
parameters (2n - 2 in a rooted tree, where n is the number of taxa).
CAGEE does not always converge to a single global maximum with large
numbers of λ parameters, and therefore can give misleading
results. To check for this you should always run the λ search
multiple times to ensure that the same estimated values are found. Also,
the likelihood of models with more parameters should always be lower
than models with fewer parameters, which may not be true if [CAGEE
]{}has failed to find a global maximum. If CAFE does not converge over
multiple runs, then one should reduce the number of parameters estimated
and try again.

Logging
-------
More verbose logging can be provided with an EasyLogging log config file.
The file may look like this:

	* GLOBAL:
	   FORMAT               =  "%datetime %msg"
	   FILENAME             =  "cafe.log"
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

Pass the config file to CAFE with the --log_config flag. For example,

    CAGEE -t examples/mammals_tree.txt -i filtered_cafe_input.txt --log_config log.config	   

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

When optimizing for an alpha value with a set number of clusters, if the
largest multiplier in the longest branch is saturated, the scorer will
return an infinite value. This will be noted at the end of the run with
text like:

    90 values were attempted (10% rejected)

showing that 10

    The following families had failure rates >20% of the time:
    Family6 had 22 failures
    Family9 had 19 failures

Certain options are available at compile-time for the optimizer. If
OPTIMIZER\_STRATEGY\_INITIAL\_VARIANTS is defined, the optimizer will
take several shorter attempts at various values before settling on one
value to continue on with. This may cause the optimizer to take more
iterations to finish but may have greater accuracy. If
OPTIMIZER\_STRATEGY\_PERTURB\_WHEN\_CLOSE is defined, the optimizer will
begin searching a wider range of values when it is getting close to a
solution. This attempts to get the optimizer out of a local optima it
may have found.

How does the simulator choose what sigma to use?
-------------------------------------------------

Although the user specifies the sigma, in order to give more family
variety a multiplier is selected every 50 simulated families. So
if 10,000 families are being simulated, 200 different sigmas will be
used.

When simulating without the gamma model, the multiplier is a random
value based on a normal distribution with a mean of 1 and a standard
deviation of 0.3.

When simulating WITH the gamma model, the multiplier is drawn directly
from a gamma distribution based on the selected alpha and a mean of 1.
If clustering is requested via the -k parameter, the selected cluster
multiplier is modified by a normal distribution with the mean at the
value of the multiplier, and a standard deviation intended to reflect
the number of clusters requested.

Initial Guesses
---------------

One of the most important concerns when searching a parameter space is
what initial values to choose. Each of the three values that CAFE can
search for has a particular initial guess strategy. For sigma values,
the formula is 1 / (longest tree branch times a random number between 0 and 1).
For epsilon values, the initial guesses are taken directly from the
provided error model. For gamma values, a random value taken from an
exponential distribution is used.

In some situations the values may fail. In this case, the scorer will
return an infinite value and the optimizer will retry initialization,
up to a number of attempts determined at compile time. If, after this number
of attempts, the optimizer continues to fail, it will abort. The most
likely cause of failure is too wide a variety of species sizes inside
certain families, and a message will be shown giving the most likely
families to remove from the analysis for success. 	

Acknowledgements
================

Many people have contributed to the CAFE project, either by code or
ideas. Thanks to:

-   Ben Fulton

-   Matthew Hahn

-   Mira Han

-	Fabio Mendes

-   Gregg Thomas

-	Dan Vanderpool

CAFE uses the EasyLogging logging framework. https://github.com/amrayn/easyloggingpp

CAFE uses the DocTest testing framework. https://github.com/onqtam/doctest

