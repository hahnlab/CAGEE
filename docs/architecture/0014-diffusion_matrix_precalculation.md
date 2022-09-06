# 14. Precalculate a diffusion matrix to speed up processing time

Date: 2022-07-06

## Status

Accepted

## Context

Each time CAGEE starts, it calculates a diffusion matrix to work with, that is only dependent on the
discretization size specified by the user. If we precalculated this matrix it might save the user 
quite a bit of time.

Although the discretization size defaults to 200, we allow the user to set it to an arbitrary size so
we could not precalculate every size they might use. We could have CAGEE write out each matrix that
it calculates, but there is no guarantee that the application will be installed into a directory where
the user has write access. We could write it to an application-specific directory such as $HOME/.cagee
or similar, but since we don't require much in the way of configuration files, this seems like overkill.

## Decision

We will provide a separate application, diffmat_calc, which will calculate the diffusion matrix for a
given size and write it to a pair of files named "eigenvaluesx.bin" and "eigenvectorsx.bin" where the x
will be replaced by the discretization size. CAGEE will search for the file in the current directory and 
load it if available.

## Consequences

A user who does not read the directions may not take advantage of this optimization. It might be a good
idea to provide a precalculated matrix for the default size as many users probably will never change it.


