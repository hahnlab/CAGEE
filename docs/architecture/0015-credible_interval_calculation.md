# 15. Calculate a credible interval for each reconstructed value

Date: 2022-09-06

## Status

Candidate

## Context

Reconstructed values, while representing the most likely values based on the 
discretized likelihoods calculated from the diffusion matrix, may still be 
more or less accurate based on the values. 

To avoid over-interpretation of small changes in inferred expression levels 
- especially when there is uncertainty in ancestral states - CAGEE will also 
compare the credible intervals at parent and child nodes to note if a change 
is “credible” (i.e. the intervals do not overlap). 

The credible interval will be calculated as follows:

Start from the reconstructed distribution, a vector of probabilities. To get
 the Qth percentile we first need to calculate the cumulative probabilities
(the sum of all values in the probability vector up to that position).  So 
that if the vector of probabilities is `[p1,p2,p3]` then the cumulative 
probabilities are `[p1,p1+p2,p1+p2+p3]`. These values must then be normalized 
such that the final value is equal to 1, so each value is divided by that 
final value. Each value in this vector is the total probability up to that 
point. From this vector, any percentile Q can be computed by finding the 
largest value in the vector that is smaller than Q. The expression value 
corresponding to that percentile is the value in the likelihood vector at 
the same index.

As with the most likely value, we will take the midpoint of the bin to be 
the final value.

## Decision

We will provide a credible interval reflecting the range of values for which 
we have a 95% confidence that the true value lies within. To compute the 95% 
credible interval we will take the expression value for Q=0.025 as the lower 
limit and the value for Q=0.975 as the upper limit.

We will provide a new output file, *base_credible_intervals.tab*, which will 
hold a table of values similar to the base_change and base_level files. 
For each node-family combination, a credible-interval will be output in 
the format `[lower-upper]`. Thus the file might look like

    TranscriptID	cat<1>	horse<2>	cow<3>
    transcript0	[45.67-92.73]	[564.23-893.54]	[12.09-44.6]	
    transcript1	[323.76-402.11]	[15.99-29.23]	[114.12-517.1]

Since the leaf values are known, their credible interval will simply be 
the known value repeated for the upper and lower bounds.

We will modify two of the other output files as well: *change.tab* and 
*clade_results.tab*. *Change.tab* will be modified to show an asterisk by
changes for which the credible intervals do not overlap. *clade_results.tab*
will no longer count results for which the credible intervals do overlap
(i.e. the predicted change is likely to represent a real change).
 
## Consequences

The word "Credible" is open to some interpretation. One might consider that
a change in value which lies inside the credible interval would be a "credible"
change, and a value outside that range would be a "non-credible" change.
We use the term in the opposite sense; a change that overlaps intervals is likely
to not represent a real change, thus, those changes that lie outside the interval
are considered "credible changes".
 
The user might not understand the reasoning for the new file. The format of the 
credible intervals file may be overly difficult to parse.


