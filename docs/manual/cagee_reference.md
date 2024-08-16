## Command Reference

Here is a list of commands that can be used with the CAGEE tool:

### -c, --config <arg>
 Configuration file containing additional options. File should consist of key-value pairs:

```
   tree = example_tree.txt 
   infile = gene_transcripts.txt
```

### -@, --cores <arg>
 Number of processing cores to use, requires an integer argument. Default=All available cores.

### --count_all_changes [=arg(=1)]
 Reconstruction will count all changes rather than only credible changes

### -D, --discretization_size <arg> (=200)
 Size (length) of the discretization vector. Default=200. Can increase resolution at the cost of computation time.

### -e, --error [=arg(=true)] (=false)
 Run with no file name to estimate the global error model file. This file can be provided in subsequent runs by providing the path to the Error model file with no spaces (e.g. -eBase_error_model.txt).

### --fixed_alpha <arg>
 Value for a single user provided alpha value, otherwise alpha is estimated.

### -m, --fixed_multiple_sigmas <arg>
 Multiple sigma values, comma separated.

### -l, --fixed_sigma <arg>
 Value for a single user provided sigma value, otherwise sigma is estimated.

### - -h, --help
Produce help message

### -i, --infile <arg>
 Path to tab delimited gene families file to be analyzed (Required for estimation)

### -@, --n_gamma_cats <arg>
 Number of gamma categories to use, requires an integer argument. Default=1 (No gamma modelling)

### -E, --optimizer_expansion <arg>
 Expansion parameter for Nelder-Mead optimizer. Default=2.

### -I, --optimizer_iterations <arg>
 Maximum number of iterations that will be performed in sigma search. Default=300 (increase this number if likelihood is still improving when limit is hit).

### -R, --optimizer_reflection <arg>
 Reflection parameter for Nelder-Mead optimizer. Default=1.

### -o, --output_prefix <arg>
 Output directory - Name of directory automatically created for output. Default=results.

### --prior <arg>
 Expected distribution of the root in Inference Mode (mutually exclusive with --rootdist in Simulation Mode). Can be specified as one of three distributions:

 * ```gamma:k:theta```
 * ```fisher:d1:d2```
 * ```uniform:lower:upper```

Default is ```gamma:0.375:1600.0``` or ```fisher:0.75:0.75``` if the ratio flag is specified.

Note: The *uniform* distribution has no basis in biology but is available for testing and verifying simulated data.

### --replicate_map <arg>
 Filename of a file containing a list of specie replicates to be combined into a single species

### --rootdist <arg>
 Distribution of the root in Simulation Mode (mutually exclusive with --prior in Inference Mode). Can be gamma:[k]:[theta], fixed [count], or a path/to/tab_sep_file.txt with two columns: the root value to use and the frequency of it. Default=gamma:0.375:1600.0

### --sample_group <arg>
 Specifies sample groups (if any) for which to infer sigma^2. Each sample and sigma^2 estimate requires a --sample_group [your_sample_A] arg, or combine them with comma
 --sample_group [your_sample_A,your_sample_B,...]. Optional, no default.

### -y, --sigma_tree <arg>
 Path to sigma tree, for use with multiple sigmas

### -s, --simulate <arg> (=false)
 Simulate families. Optionally provide the number of simulations to generate

### -t, --tree <arg>
 Path to file containing newick formatted tree (Required for estimation)

### --unbounded [=arg(=1)]
 The input file contains ratios of gene expression values rather than absolute values

### -v, --version
Print version string

### --verbose <arg>
Extended logging. Arg must be an integer between 1 and 9.

### -z, --zero_root [=arg(=1)]
 Exclude gene families that don't exist at the root, not recommended.

