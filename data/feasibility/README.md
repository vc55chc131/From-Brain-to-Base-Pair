# Spatial feasibility audit

This dataset records a completed audit of three proposed 2 mm radius ROIs using the spatial metadata for all six Allen Human Brain Atlas donors. It contains no gene-expression measurements. The proposed centres are left DLPFC at -45, 36, 28; left Caudate at -12, 15, 9; and right Caudate at 12, 15, 9 in MNI space. These are the coordinates under examination, without a claim that they reproduce a published fMRI contrast.

The inputs contain 3,702 tissue coordinates, Allen sample annotations and ontology information, and a saved API response that maps well IDs to donors. Every corrected-coordinate well ID maps uniquely to an API donor ID, and each donor's well-ID set agrees exactly with the corresponding SampleAnnot.csv file. Original and corrected MNI coordinates are evaluated separately.

## Reproduce the audit

Install the optional environment specified by `requirements-ahba.txt`. From the repository root, run:

```bash
python source_code/audit_roi_coverage.py --inputs data/feasibility/inputs --seeds config/roi_seeds_unverified.csv --output-dir outputs/roi_coverage_audit
```

The output directory must be empty. The script loads only coordinates and anatomical metadata. It does not load, download, process or infer gene-expression values. `results/audit_manifest.json` records the actual software versions, input and output checksums, and source-code checksums for the deposited run.

## Completed results

The audit distinguishes exact Euclidean spheres around the coordinates from the voxelized atlas assignment performed by abagen. The latter uses a 2 mm resolution MNI152 grid with labelled voxel centres inside the proposed spheres. Abagen projects sample coordinates onto its grid before tolerance matching. Consequently, its zero-tolerance result is not identical to a continuous 2 mm distance test.

| Coordinate set | Assignment | Left DLPFC | Left Caudate | Right Caudate | Paired donors |
|---|---|---:|---:|---:|---:|
| Corrected MNI | Exact closed 2 mm spheres | 0 | 0 | 0 | 0 |
| Corrected MNI | Voxel atlas, tolerance 0 mm | 0 | 0 | 0 | 0 |
| Corrected MNI | Voxel atlas, tolerance 2 mm | 0 | 1 | 0 | 0 |
| Original MNI | Exact closed 2 mm spheres | 0 | 1 | 0 | 0 |
| Original MNI | Voxel atlas, tolerance 0 mm | 0 | 0 | 0 | 0 |
| Original MNI | Voxel atlas, tolerance 2 mm | 0 | 4 | 0 | 0 |

These counts are identical before and after the implemented hemisphere and broad structural-class restrictions. The corrected-coordinate tolerance 2 mm assignment is well 11333 from donor H0351.1009, whose Allen anatomical annotation identifies the left head of the caudate. The nearest corrected tissue coordinates are 9.405998 mm from the proposed DLPFC centre, 2.788991 mm from the left Caudate centre, and 4.536933 mm from the right Caudate centre. These are distances to the proposed centres, not distances to labelled voxels.

No donor has an observed DLPFC–Caudate pair under any audited condition. The proposed regional expression comparison therefore cannot proceed under these specifications. This is a finding about the tested sampling geometry and assignment rules. It does not establish absence of DLPFC expression in the atlas or absence of a biological difference between the regions.

## Files and provenance

`results/donor_counts.csv` reports donor-specific counts for every condition. `coverage_summary.csv` aggregates those counts. `nearest_samples.csv` gives the nearest observed sample to each seed within each donor for both coordinate sets. `matched_samples.csv` lists every nonzero assignment and its original Allen well ID. `sample_totals.csv` reports the available tissue samples before and after the hemisphere mismatch filter. The atlas and its specification are retained in `seed_atlas.nii.gz` and `atlas_info.csv`.

The corrected coordinates were copied from the pinned upstream abagen commit and verified against the public source by SHA-256. The BSD-3 notice for their alleninf source is reproduced in `inputs/CORRECTED_COORDINATES_LICENSE.txt`. The original sample annotations and ontology were retrieved from the official Allen donor archives through bounded HTTP byte-range requests. Only SampleAnnot.csv and Ontology.csv were transferred and checked by ZIP CRC; no expression arrays were downloaded. The six ontology files were byte-identical, so one shared copy is retained. Input URLs, archive sizes, transferred byte counts and checksums appear in `inputs/input_provenance.json` and the donor-specific download provenance files. The code license does not replace the source-data terms or the separate corrected-coordinate license.
