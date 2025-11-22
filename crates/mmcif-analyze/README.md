# mmcif-analyze

Spatial search helpers that sit on top of `mmcif-core`. The crate mirrors the
API style of pdbtbx/Biopython but keeps the implementation local so the
monorepo can stay dependency-light.

Current focus:
- Query for the closest atoms around an arbitrary (x, y, z) coordinate, honoring a
  distance cutoff and a maximum number of hits.
- Query for the closest amino-acid residues (based on `_atom_site.label_comp_id`)
  using minimum-image wrapping when unit-cell information is available.
