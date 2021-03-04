# 8. Initial guess for Sigma

Date: 2021-03-04

## Status

Candidate

## Context

The Nelder-Mead optimizer performs better if given an initial guess that will approximate the final result.

## Decision

The initial guess will be calculated as follows:

* Take the variance of all of the leaf values of all of the gene transcripts.
* Multiply that by the length of the tree
* Take the square root of that value
* Select the initial guess from a normal distribution with the mean of that value and stddev of 0.2

## Consequences

We use the fact that the observed variance is proportional to sigma\*sqrt(t).  In principle we just need 
to get to the right order of magnitude otherwise we probably shouldn't be trusting the optimization anyway.
Note that we don't need a known root value to calculate the guess. That's just showing that no matter if we 
happen to start close to the lower bound (x=1) or away from the bounds (x=10), the result will be roughly 
the right order of magnitude.
