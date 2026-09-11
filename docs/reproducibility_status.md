# Reproducibility status

Updated 11 September 2026.

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

The new paired pipeline is a prospective reanalysis implementation. It does not establish which method produced the original paper's reported values. It deliberately does not tune settings to reproduce those values. Read `docs/manuscript_reconciliation.md` before treating the manuscript as ready for submission.

## Next scientifically necessary step

Resolve the coordinate provenance and either supply the original processed dataset with its records or execute a clearly identified new extraction. Preserve the full outputs and all failed/insufficient-coverage diagnostics. Use the paired files exported by extraction for the companion analysis, and report any excluded donors. A completed, reviewable empirical run must precede any claim that the reported research findings are reproducible.

Public GitHub availability, a software test suite, and this status document do not by themselves establish full FAIR compliance.
