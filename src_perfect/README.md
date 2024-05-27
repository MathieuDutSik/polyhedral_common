Perfect form enumeration
========================

This is the same business as for IsoDelaunay domains,
IsoEdge domains and others.

However, due to the particular importance of perfect forms,
we ought to think carefully about what we want to do.

What we want from this is the following:
* Work in a T-space. That is super important. So that means
the space of the perfect domains and the Ryshkov space of the
perfect forms. This is like that in Opgenorth paper.
* The Ryshkov stuff allows to do the enumeration quite
efficiently. On the other hand, the perfect domain are better
for computing cohomology.
* The finding of initial form and the flipping are quite clear
from existing code. The equivalence and stabilizer can also
be done from the Tspace equivalence code. That part is clear.
* How to build the dual cone? Possible sources:
  --- The TechVor.pdf unpublished manuscript.
  --- The GL4(Z[i]) paper.
  --- The GAP code, there are several.
* What we want from the dual cone:
  --- That it embeds into S^n_{>0} so that the existing tools
  for stab/equi can be used.
  --- That the linear form of evaluation correspond to a vector
  in the dual.
* What to do:
  --- It is a little bit unclear whether we can have a general theory
  --- It would make sense to have the dual space put as argument. If
  there is no universal theory, then so be it.
  ---


