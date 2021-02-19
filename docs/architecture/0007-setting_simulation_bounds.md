# 7. Setting simulation bounds

Date: 2021-02-19

## Status

Accepted

## Context

When simulating gene families, a boundary needs to established in order to discretize the range of values.
A simulated gene famuly starts with one root value, which does not give, a priori, the upper limit to be
set in your state vector. 

## Decision

The boundary will be calculated based on the root size of the family being simulated. The formula will be:

root size + 4.5 * sigma / sqrt(t);


## Consequences

The thought is to do one calculation: multiplying the root value by t (the distance from the root to tips) and 
sigma^2, to get the upper 99% tail of this distribution (and then multiply that by 1.5). Using the root-tip 
distance is a good idea. If we assume unbounded Brownian motion then the variance is of a trait value about 
its root value after time t is equal to t \* Sigma^2. 3 times the square root of that is the 99% distance, so 
an upper bound of 1.5 \* 3 \* sigma / sqrt(t) should be safe. As a safety check, we should be sure that the transition 
matrix from root to tips at the rows corresponding to the highest allowed trait values is negligible.

It is not clear if the root value has to be in the right-side as well. Variance may be dependent on starting value.
We will dynamically vary the bounds for each gene family. But we may be assuming some sort of normalization across 
10,000 genes. 

One possibility is that variances are proportional to magnitudes which seems biologically plausible.
In that case we're really modelling log trait values. We wouldn't need to change the model, just take logs on the 
data. That is, things happen multiplicatively rather than additively (as in the stock market where percentage differences
are what matter, not absolute differences). In the birth-death model, variance is directly related to starting value,
so for the simulated data we just need to remember that the root value is actually the log of the actual trait value.


