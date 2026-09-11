# From Brain to Base Pair

Reproducibility materials for an exploratory comparison of dorsolateral prefrontal cortex and caudate transcriptomic profiles, motivated by simultaneous interpreting research.

## Current status

A **completed spatial feasibility audit** now precedes the conditional expression workflow. Across 3,702 corrected AHBA sampling coordinates, exact 2 mm spheres at the author-specified DLPFC and bilateral caudate centers contained no sampling centers. Matching the corresponding voxel atlas with no additional tolerance assigned no samples. A diagnostic 2 mm additional tolerance assigned only one left caudate candidate and no DLPFC candidates. Both configurations therefore yielded **zero complete regional donor pairs**. These settings do not support expression inference.

This result applies to the specified coordinate resource, centers, atlas grid, and tolerances. It does not mean that the anatomical DLPFC or caudate generally lack AHBA samples, and it does not test molecular equivalence between regions. Tissue centers do not represent the full physical extent of dissected samples.

- `source_code/audit_roi_coverage.py`: executed spatial coverage assessment.
- `data/feasibility/`: source inputs, provenance, exact-sphere and voxel-matching outputs.
- `source_code/DLPFC_vs_Caudate_Zscore.py`: conditional paired analysis of a real processed expression matrix and donor/ROI metadata.
- `source_code/extract_ahba.py`: conditional extraction from public AHBA data; it has not generated expression results for the audited targets.
- `requirements.txt` and `requirements-ahba.txt`: analysis and optional AHBA dependencies.
- `docs/data_schema.md` and `docs/extraction.md`: accepted inputs, parameters, and execution instructions.
- `docs/reproducibility_status.md`: completed checks and limits.
- `docs/manuscript_reconciliation.md`: disposition of the earlier manuscript discrepancies.
- `data/reference/`: official donor metadata with provenance.
- `data/author_reported/`: historical table transcriptions; **not validated results of the coverage study**.
- `tests/`: computational checks with synthetic fixtures, not research observations.

The original full processed expression matrix and executed enrichment outputs remain unavailable. The original coordinates are evaluated as exploratory author-specified inputs, not verified peaks from Hervais-Adelman et al. (2015). The clean manuscript now reports the spatial feasibility study and does not retain the unsupported gene-level or enrichment conclusions.

See `data/feasibility/README.md` for the exact audit command, input definitions, and source terms.

## Install and run

Use an isolated Python environment. The new analysis implementation uses donor-paired normalized expression differences. It does not reproduce the old example's across-gene standardization or calculate log2 fold changes from normalized scores.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python source_code/DLPFC_vs_Caudate_Zscore.py --help
python -m unittest discover -s tests -v
```

On Windows, activate the environment with `.venv\Scripts\activate`.

Provide actual analysis data in the documented schema, then run:

```bash
python source_code/DLPFC_vs_Caudate_Zscore.py \
  --expression data/processed/expression_matrix.csv \
  --metadata data/processed/sample_metadata.csv \
  --output outputs/run_20260911 \
  --normalization-description "REPLACE with the actual preprocessing and scale of your input matrix"
```

Replace the normalization description with a factual account of the input matrix. Add `--heatmap` for a visualization and `--enrich` to request functional enrichment through g:Profiler. The enrichment option sends selected gene symbols and a background gene list to the public g:Profiler service and requires network access. Its outputs depend on the service's annotation release. The response and provenance must be retained with the run.

The optional AHBA route is described in [docs/extraction.md](docs/extraction.md). Extraction downloads substantial public source data. It must not be mistaken for a completed analysis merely because the script exists.

## Outputs and interpretation

Only executed commands generate numerical analysis outputs. The workflow records input checksums and settings, aggregates within donors before a paired regional comparison, distinguishes descriptive top-ranked genes from statistically selected genes, and avoids unsupported fold-change terminology. Read the generated manifest and inspect donor representation before interpreting a result.

The legacy five-gene matrix `AHBA_expression_matrix_sample` is a formatting example. Its placeholder sample IDs and values are not study observations and cannot reproduce the manuscript. Legacy instructions under `Repository Contents/` are superseded by the documented workflow above.

## Scientific scope

AHBA donors were not selected as interpreters. A regional atlas contrast does not identify an interpreter-specific molecular signature, a training effect, an aptitude marker, or a causal mechanism of language control. Expression abundance also does not directly measure neurotransmitter activity. Anatomical matching, donor dependence, and gene-set definitions require explicit documentation.

## Data and code availability

Analysis code, documentation, source donor metadata, coordinate inputs, and executed spatial feasibility outputs are available in this repository. Historical author-reported table transcriptions are labeled separately. Original transcriptomic data are available through the [Allen Human Brain Atlas](https://human.brain-map.org/). The complete original processed expression dataset and original executed analysis outputs are not currently included. No claim of complete reproduction or full FAIR compliance is made.

## Sources and license

- AHBA: Hawrylycz et al. 2012, https://doi.org/10.1038/nature11405
- abagen: Markello et al. 2021, https://doi.org/10.7554/eLife.72129
- g:Profiler: Raudvere et al. 2019, https://doi.org/10.1093/nar/gkz369
- Imaging source discussed in the manuscript: Hervais-Adelman et al. 2015, https://doi.org/10.1093/cercor/bhu158

Repository code is covered by the existing MIT `License`. Upstream datasets and services remain subject to their own terms; the repository license does not relicense AHBA data. Cite the original atlas and methods alongside this repository. The research manuscript has not been assigned a publication DOI here.
