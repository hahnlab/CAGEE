# 10. Estimation of different sample types

Date: 2021-10-27

## Status

Candidate

## Context

A user should be able to provide several sample types in a single transcript file. Each sample type
will be included in a sample group and a separate sigma will be estimable for each sample group.

We need to specify what to do if the user has sample types that are not part of one of the sample
groups, and whether the groups need numerical identifiers or labels.

## Decision

The sample type will be the third column of the transcript file, after the description and ID fields,
but before the transcript values. The column header must be SAMPLETYPE; otherwise the column will be
read as a transcript value.

A parameter, --sample_group will be added. The user will be able to specify the parameter multiple
times. The parameter will support comma-separated samples, so the arguments should look like this:

$ cagee --sample_group heart,lungs --sample_group brain

If there are samples that the user does not specify in a sample group, they will be placed in a separate 
group, and a warning message will be issued.

An alternative would be to support a sample_group_count parameter and require that the user provide
the exact number of sample groups in the parameter. This seems to add additional cognitive load for the
user.

We might also support naming the groups, probably by making the parameter accept a value like 
"group1:heart,lungs". We will not support that at this time.

## Consequences

We now support two kinds of multiple sigmas: a lineage-specific sigma for different parts of the tree,
and a sample-specific sigma for different sample groups. These two kinds of sigma estimations may be
confusing.

Users may forget to place a sample in a specific group, and be surprised when those samples are not
added to any of their existing groups. The warning message will hopefully help with this.
