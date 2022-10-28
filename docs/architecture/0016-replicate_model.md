# 16. Replicate model

Date: 2022-10-28

## Status

Candidate

## Context

We want to implement a replicate model to perform calculations on it.

A replicate model is one in which, rather than a species having a single value associated
with it, it has a series of values, or replicates.

One possible method of implementing this is to provide a mapping file as an argument. Each
replicate would be specified as a separate column (i.e. cow-1, cow-2, cow-3) and the mapping
file would identify which species are associated with which columns (cow => cow-1, cow => cow-2).

Another possiblity would be to allow replicates to be specified as comma-separated values in our
input file. We require the input to be in TSV format, so these values should not cause parsing errors.

### Implications
* The uppper bound should not need to reflect the exact values of each replicate - taking the average value should be sufficient
* When computing the CAGEE estimation, the replicates should be added together to create the initial value

## Decision



## Consequences
