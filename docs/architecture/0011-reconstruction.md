# 10. Reconstruction 

Date: 2022-01-12

## Status

Candidate

## Context

The application should be able to create the most likely set of expression values throughout the 
internal nodes of the tree. We could attempt to modify the Pupko algorithm used in CAFE 5 for
continuous values, or the CAFE 4 algorithm similarly, or come up with a new strategy. The CAFE 4
algorithm is difficult to glean from reading the code, however.

## Decision

We will attempt to recreate the CAFE 4 algorithm, modified for continuous values, as best we can.
Each node will be assigned a vector of likelihoods of size DISCRETIZATION_RANGE. Leaf nodes will
be initialized using VectorPos_bounds. Internal nodes will be calculated as follows: Compute
the relevant branch length matrix using ConvProp_bounds for each child, and multiply the matrix
times the child's likelihood vector. The internal node will then be the mean value of these two
vectors. To recreate a value, invert the discretization process at the maximum value of the vector -
that is (max/DISCRETIZATION_VALUE) * (the upper bound used for discretization).

## Consequences

None yet.
