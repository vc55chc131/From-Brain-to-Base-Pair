# Input data schema

The analysis accepts two CSV files. These must contain actual processed study data, not the legacy five-gene example. No example values should be promoted to research observations.

## Expression matrix

`expression_matrix.csv` has a first column containing gene identifiers, preferably named `gene_id`, followed by one column per sample or explicitly documented donor-region summary. The script accepts any nonempty first-column header and normalizes the output name to `gene_id`. Gene identifiers must be unique. Entries must be finite numeric values on one consistent, documented normalized expression scale. Column names must match `sample_id` values in the metadata exactly, including any leading zeros.

The pipeline computes differences on the supplied scale. A difference between normalized values is not a log2 fold change. Probe selection and normalization belong to the preprocessing provenance and must be described separately. Use gene symbols or identifiers supported by g:Profiler if enrichment is requested; retain its mapping information and report unmapped identifiers.

## Sample metadata

`sample_metadata.csv` requires `sample_id`, `donor_id`, `roi`, and `hemisphere`. ROI values are `DLPFC` or `Caudate`. Each matrix column must have one corresponding metadata row. The same donor must contribute both regions to enter the paired analysis. A minimum of three complete donor pairs is required by the new implementation; this execution threshold does not establish adequate statistical power.

For tissue-level input, `sample_id` must identify the actual atlas tissue sample. `input_level` may be omitted or set to `tissue_sample`; `hemisphere` must be `L` or `R`. Each row represents one tissue sample. Retain additional columns for MNI coordinates, ROI center, matching distance, source hemisphere, and source-file identifiers whenever available. The script does not infer these fields from a synthetic numbering convention such as 1xxx/2xxx.

For data exported as donor-region summaries by the extraction script, IDs identify derived summaries rather than individual tissues. Set `input_level=donor_roi_summary` on every row and provide positive integer `n_tissue_samples` values from observed extraction counts. There must be exactly one summary row per donor and ROI. `hemisphere` may be `L`, `R`, or `B` for a bilateral summary. The implementation preserves `n_tissue_samples` separately from `n_input_rows`. The aggregation rule for bilateral caudate values must be stated; do not interpret hemisphere summaries as independent donors. Mixing tissue-level rows and summary rows in a single run is rejected.

## Required provenance

Retain the original atlas release or download manifest; input checksums; software versions; probe-selection and normalization settings; gene-filtering rules; donor and hemisphere handling; ROI definitions and their provenance; tissue matching and missing-data policies; and all gene sets submitted for enrichment. An output filename or a prose description is not a substitute for the corresponding file.

The files under `data/author_reported/` are transcriptions for reconciliation. They are not inputs to the gene-level analysis and must not be used to reconstruct unavailable donor measurements.
