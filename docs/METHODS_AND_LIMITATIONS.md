# Methods and limitations

Adult cortical labels were transferred on each hemisphere using the verified existing surface correspondence. V4 uses local 1-mm graph closing and small-hole proposals, requiring adult-label support, at least two of three ipsilateral age-GM ribbon samples and no geometric caution. Mild area-weighted categorical voting protects small source components and limits per-ROI area removed during smoothing to 5%. Existing mesh geometry is retained.

Native labels are resampled to the paired 32k age surfaces and projected through the original native white/pial ribbons onto the original approximately 0.5-mm T1 grid. Bilateral overlap ownership uses the adult hemisphere prior transformed into age space, with age segmentation and reliable surface distances. All 112 overlap voxels have one selected owner, while 72 retain uncertainty flags.

The original 1,213 native and two 32k unassigned cortex vertices remain. Original pial and white/pial intersections persist. The conspicuous six-vertex Group04 right medial hole remains excluded in the primary atlas. Group02 left PHG.mi's single-vertex secondary source component was already absent at 32k and remains outside support. Source components represented in the v3 baseline were protected in v4.

Technical integrity and 248-label completeness passed, but anatomical status remains **HOLD_SOURCE_GEOMETRY_AND_LOCAL_REVIEW**. These atlases implement adult ROI definitions on age-template anatomy. They do not independently discover developmental boundaries or provide individual-level measurements. The cosmetic companion is explicitly excluded from quantitative use.
