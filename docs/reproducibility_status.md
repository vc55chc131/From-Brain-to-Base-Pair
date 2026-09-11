# Reproducibility status

Updated 11 September 2026.

## Completed spatial feasibility study

The executed audit in `source_code/audit_roi_coverage.py` is supported by the inputs, output tables, run metadata, source terms, and provenance under `data/feasibility/`. Independent coordinate and abagen matching checks agreed:

- 3,702 finite corrected sampling coordinates with unique well identifiers.
- Exact 2 mm sphere counts: zero at all three targets.
- Minimum center distances: DLPFC 9.405998 mm; left caudate 2.788991 mm; right caudate 4.536933 mm.
- Zero additional abagen tolerance: no assigned samples.
- Additional tolerance of 2 mm: one left caudate candidate (well 11333, H0351.1009); no DLPFC candidate and no right caudate candidate.
- Zero complete DLPFC–caudate donor pairs under both settings.

The completed result is a spatial feasibility finding. It is not a gene-expression null result or a reproduction of the old statistical summaries. Discrete atlas assignment and continuous seed distance have different definitions; abagen grid mapping must not be described as simply adding the tolerance to the sphere radius. No full microarray comparison or enrichment analysis was executed for these targets.

## What was added and checked

- A new executable processed-data analysis pipeline with donor-first regional aggregation, paired inference, Benjamini–Hochberg adjustment, descriptive rankings, optional heatmaps, optional g:Profiler calls, and input/run provenance.
- Eleven core tests passed. They cover identifier preservation, duplicate and missing data guards, donor pairing, unequal tissue-sample counts, known BH results, paired t statistics/probabilities/confidence limits, undefined zero-variance tests, arithmetic overflow, summary input, and invalid enrichment probabilities. Test fixtures are synthetic and are not study data.
- Core tests passed with both NumPy 2.3.5 / SciPy 1.17 and NumPy 1.26.4 / SciPy 1.12 environments. Exact installed versions are recorded by actual runs. Python 3.10 syntax was checked; a full Python 3.10 execution was not claimed.
- Five optional extraction tests passed in an isolated AHBA environment. They cover the unverified-seed guard, the pinned abagen aggregation API and missing regions, observed-count weighting, paired-file integration with the companion script, and rejection of zero-coverage data. Atlas preparation was also checked separately without microarray downloads. These checks do not validate biological results.
- Donor age, sex and postmortem interval were retrieved from the Allen Institute API. The saved response, derived CSV and checksum provenance are in `data/reference/`. Discrepancies between API PMI values and the abagen packaged donor table are documented.
- All twenty GO table rows and all nine visible selected-gene rows from the supplied table were transcribed to `data/author_reported/`. They are labeled unverified author-reported results and were not regenerated.
- The old README and legacy instruction files were corrected to describe actual paths, current availability and the distinction between examples and study data.

## What remains unexecuted or unavailable

The following are **not** established by this update:

- Full download, processing and extraction of the six AHBA microarray datasets.
- Biological/anatomical verification of the author-reported ROI coordinates.
- Reproduction of the original approximately 16,000-gene matrix, gene-level q values, 43-gene subset, standardized differences or enrichment probabilities.
- A live g:Profiler analysis of verified study inputs.
- Recovery of original source tissue mappings, donor ROI coverage, historical scripts or executed logs.
- Validation of the original heatmap, duplicate enrichment images, or claimed cross-pipeline/PCA/semantic analyses.

The new paired pipeline is a prospective reanalysis implementation. It does not establish which method produced the original paper's reported values. It deliberately does not tune settings to reproduce those values. The disposition of the original reporting issues is recorded in `docs/manuscript_reconciliation.md`.

## Scope of a future expression study

An empirical expression study requires a new, independently justified anatomical design with adequate observed coverage, or a recoverable original dataset with its authentic mapping and processing records. Do not enlarge regions merely to recover a desired significance result. Any changed coordinates or spatial rules define a new analysis and need explicit reporting. Preserve all coverage diagnostics, excluded donors, input lists, and executed results.

The revised manuscript is a methodology/feasibility paper based on the completed spatial audit. The earlier expression and enrichment claims have been removed rather than relabeled as reproduced findings.

Public GitHub availability, a software test suite, and this status document do not by themselves establish full FAIR compliance.
