# xpam: Presence-Absence Matrix Construction from Species Range Maps

Constructs presence-absence matrices (PAMs) from species range polygons
or occurrence records. A streaming native engine rasterizes geometries
onto a user-defined grid or bins them in environmental space, and
returns the sites-by-species incidence as a sparse matrix plus
marginals; a native Jaccard backend computes species-by-species or
site-by-site similarity directly from the same incidence, at scale. A
small-input dense path is kept for convenience. Includes convenience
functions to retrieve example mammal range data and WorldClim
environmental layers projected to two dimensions via Principal Component
Analysis.

## See also

Useful links:

- <https://alrobles.github.io/xpam/>

- <https://github.com/alrobles/xpam>

- Report bugs at <https://github.com/alrobles/xpam/issues>

## Author

**Maintainer**: Angel Robles <a.l.robles.fernandez@gmail.com>
([ORCID](https://orcid.org/0000-0002-4559-1751))

Authors:

- Angel Robles <a.l.robles.fernandez@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-4559-1751))
