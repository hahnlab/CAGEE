# 13. Setting estimation bounds

Date: 2022-06-28

## Status

Candidate

## Context

When estimating gene families, a boundary needs to established in order to discretize the range of values.
For simplicity's sake, it might be useful to use the same calculation as used in simulations, that is
the largest value plus the formula

4.5 * sigma * sqrt(t);

(We add this value rather than multiplying it since all of our values are in log space). If we take this
approach, we have to decide when the upper bound is calculated, since as the optimizer runs, it uses varying
values of sigma to find the most appropriate one. There are at least three options: We could use the _initial_ 
sigma that is selected; we could _recalculate_ the upper bound each time the sigma changes; or we could use
the initial sigma to start with and calculate a new upper bound once a final sigma has been found.

Alternatively, since we have the tip values available, we can calculate an upper bound based on these. It seems
reasonable that the upper bound will not ever be much larger than the largest tip value, so we can simply take
that value plus a constant.

## Decision

The upper bound will be selected as max(1.0, ceil(gt.get_max_expression_value() + 3.0));

## Consequences

In testing, changing the upper bound at each iteration can cause incorrect estimates. Since we increment in 
steps of 5, the bounds were suddenly jumping back and forth on each estimate, causing the optimization to
go in incorrect directions. It seems odd to use one upper bound at the beginning and another at the end, and
equally odd to use an upper bound calculated from a sigma that is not being used.

The 3.0 is an arbitrary value that will be defined in the program configuration so advanced users can change
it if they so desire.

