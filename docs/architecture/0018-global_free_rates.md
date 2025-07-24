# 18. Global Free Rates

Date: 2025-07-07

## Status

Candidate

## Context

On a large tree, it may be more efficient to optimize a rate of evolution for each branch separately, rather
than to calculate a single rate across the entire tree. We can optimize a node using the following algorithm:

### Global Free-Rate Algorithm
1. Assign an initial value for one child node
2. Assign a value for the other child node by finding an optimal value for it, using the sibling value
3. Now assign a value for the initial child node by finding an optimal value for it, using the optimized value for its sibling that was just calculated
4. Store the likelihoods for the node based on these calculated likelihoods

### Assigning Initial Values

The above algorithm may be sensitive to the initial values selected - if every value in the tree is initialized
to the same value, the optimization may be limited in the values it finds. We compute initial values using
the following algorithm:

1. Compute a Spearman correlation for each pair of species in the tree 
2. The inverse of the correlation will be a distance matrix
3. Use the distance matrix to compute an unrooted tree using the Fitch-Margoliash algorithm
4. Normalize the values in this tree so the median value is 1
5. Compute the distribution mean across all of the species
6. For each pair of sibling nodes, compute the distance between them in the Fitch tree. Assign the initial value to this distance multiplied by the distribution mean.

## Decision

We will allow users to calculate a rate for each branch separately if they so desire.

* We will add a flag "free_rate=global" to allow this

## Consequences

The accuracy of GFR seems to drop as distance from the tips increases. Users may find this lack of
accuracy confusing.

The algorithm does not calculate a final value at the root. Users may find this confusing as well.

