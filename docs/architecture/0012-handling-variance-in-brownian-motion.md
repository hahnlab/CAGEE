# 12. Handling variance in the Brownian motion algorithm

Date: 2022-01-26

## Status

Candidate

## Context

Brownian motion has equal variance (i.e. the average amount of change that occurs per unit time) at all values.
But we expect that the amount of evolutionary change in gene expression is related to the value - lower expressions
will vary less and larger expression values will vary more. Similar to the model underlying CAFE, the variance 
being proportional to the current value is more biologically realistic.

We could implement a stochastic model called Geometric Brownian Motion to model this proportional variance.
But it's much simpler to assume that we're not modeling gene expression itself, but rather log(gene expression)
It amounts to the same thing, except we can just stick to the familiar territory of regular Brownian motion.
We could have users convert their data to log space before sending it to CAGEE, or we could assume the data
is normal and convert the data ourselves.

## Decision

We will model log(gene expression) internally, but we will hide this from the user by converting the input data
to log space before running calculations. We will then un-log the output values before giving them to the user.
For simulations, after computing the root value of a simulated transcript, we will use the log of the computed
value. For the simulation output files, all of the values output should be un-logged before writing.

For estimations, we will take the log of the prior the user provides before using it; also, we will take the 
log of all of the values the user provides in the gene transcript input file. Most of the output files deal with
probabilities and relationships and thus don't need to be modified, but the nexus file named _model_asr.tre holds
reconstructed values, and these should be un-logged before being written. 

## Consequences

This is additional complexity in our application that may lead to bugs. But it is probably less complex than 
trying to explain variance issues to users.