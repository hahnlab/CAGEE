# 16. Calculate a credible interval for each reconstructed value

Date: 2023-01-25

## Status

Accepted

## Context

Real world data are plagued by missing data. In particular, draft incomplete 
genomes and de novo transcriptome assemblies can lack many genes. This will
create a sparse expression matrix, where many cells are not available. CAGEE
should allow users to enter their data with incomplete values.

## Decision

We will add support for missing values, signified in the input file with the 
single characters "?", "-", or "N".

When computing probabilities, we will take the following steps: 

* leaf nodes with missing data will be ignored
* internal nodes where we have not calculated probabilities for either child will be ignored.
* internal nodes where a probability vector was calculated for a single child will
be calculated by multiplying the relevant matrix by that vector.

## Consequences

Currently we do an initial check for missing data before starting the optimizer.
This will have to be removed.

The user may want to signify missing data with different characters than we
expect. We will count any non-numeric data as missing, and add a warning if
the user has data with different characters.

Output files will need to be modified as follows: Ancestral_states.tab and 
change.tab will show an "N" where data is missing; ancestral_states.tre will 
show nodes with missing data as "species<node_id>_N:branch_length"; 
credible_intervals.tab will show [N-N] for missing data.


