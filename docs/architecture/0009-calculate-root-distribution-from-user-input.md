# 9. Calculate a prior of root distributions based on user input

Date: 2021-08-31

## Status

Candidate

## Context

An expected root distribution is required in order to accurately calculate prior probabilities
and, in simulations, to select a root value for the simulated tree. The user may specify the
root distribution in one of several ways. Hopefully this decision reflects the least surprising
results for the user.

## Decision

A rootdist argument will be made available. The possible formats will be as follows:

- rootdist=gamma:alpha,beta
- rootdist=fixed:value
- rootdist=file:filename
- rootdist=estimate

If the user does not provide a rootdist argument, _estimate_ will be assumed.

The root distribution will be calculated from the sub-argument:

- gamma: a gamma probability will be used with the specified alpha and beta
- fixed: all root values will be identical. The correct prior probability is unclear.
- file: The user has specified a root distribution. Use that.
- estimate: estimate a gamma distribution from the --infile value. If no infile is specified, use a gamma distribution with defaults to be decided


## Consequences

Some users may be confused by the ordering. There may be some unforeseen circumstance where the combination of input flags is still surprising.


