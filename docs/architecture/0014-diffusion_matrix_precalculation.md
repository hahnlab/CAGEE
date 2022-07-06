# 9. Calculate a prior of root distributions based on user input

Date: 2021-08-31

## Status

Accepted

## Context

An expected root distribution is required in order to accurately calculate prior probabilities
and, in simulations, to select a root value for the simulated tree. The user may specify the
root distribution in one of several ways. Hopefully this decision reflects the least surprising
results for the user.

## Decision

For simulations, a rootdist argument will be made available. The possible formats will be as follows:

- rootdist=gamma:k:theta
- rootdist=fixed:value
- rootdist=file:filename

If the user does not provide a rootdist argument, the gamma distribution <del>(k=0.75, theta=30)</del>(k=0.375, theta=1600) will be used.
The root distribution will be calculated from the sub-argument:

- gamma: a gamma probability will be used with the specified k and theta
- fixed: all root values will be identical. The correct prior probability is unclear.
- file: The user has specified a root distribution. Use that.

For estimations, a prior argument will be made available. 
- prior=gamma:k:theta

If the user does not provide a prior argument, the gamma distribution <del>(k=0.75, theta=30)</del>(k=0.375, theta=1600) will be used.

Specifying a --rootdist argument without the --simulate argument will result in an error.
Specifying a --prior argument with the --simulate argument will result in an error.
Specifying an --infile argument with the --simulate argument will result in an error.

## Consequences

Some users may be confused by the ordering. There may be some unforeseen circumstance where the combination of input flags is still surprising.


