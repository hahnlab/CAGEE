# 17. Unbounded Brownian Motion

Date: 2023-03-22

## Status

Candidate

## Context

Some users will want to run CAGEE-style analyses on ratios of data rather than absolute values.
This implies that the logs can fall below zero so we will have to provide a mode that does
not have a lower bound of zero. However, in order to make the algorithm work, a lower bound
will still have to be defined.

Some of the ratios may be 0, or 1/0. These will have to be handled. Currently we would handle
0 cases by adding a small epsilon value before taking the log, but we might not want to do 
that.

The user data may be in a variety of formats: straight ratios, logged ratios, log10 ratios.
We only need to accept a single format, but we will need to make sure that whatever input 
we accept, our reconstruction files will be in the same format.

## Decision

We will add a flag stating that the input file consists of ratios rather than absolute values.
We will define a lower bound that is the negative of the upper bound.
We will remove the epsilon value and...do what with 0 values?


## Consequences

Resolution of the grid will be half of what was in the bounded case, e.g. (-10,10) instead of (0,10)
Right now values are assumed to be absolute values and we take their logs. If we decide to 
do something different we will have to be careful that all values are in the same space.


