# From Brain to Base Pair

Reproducibility materials for an exploratory comparison of dorsolateral prefrontal cortex and caudate transcriptomic profiles, motivated by simultaneous interpreting research.

## Current status

This repository now provides a **new, documented analysis implementation** and supporting data specifications. It does **not** establish reproduction of the numerical results in the original manuscript. The original complete processed expression matrix, donor-to-tissue mapping, executed analysis logs, and full enrichment exports were not available when this implementation was added on 11 September 2026.

The author-reported ROI coordinates are also awaiting source verification. Their attribution to Hervais-Adelman et al. 2015 is not supported by the reported peak tables inspected in that paper. Do not silently replace the coordinates and describe the resulting analysis as a reproduction.

- `source_code/DLPFC_vs_Caudate_Zscore.py`: analysis of a real processed expression matrix with required donor/ROI metadata.
- `source_code/extract_ahba.py`: optional preparation from the public Allen Human Brain Atlas, subject to its documented data-download requirements.
- `requirements.txt`: analysis dependencies.
- `requirements-ahba.txt`: optional atlas-extraction dependencies.
- `docs/data_schema.md`: accepted input formats and distinctions between tissue samples and regional summaries.
- `docs/extraction.md`: extraction parameters, coordinate provenance, and execution instructions.
- `docs/reproducibility_status.md`: what has and has not been verified.
- `docs/manuscript_reconciliation.md`: discrepancies that must be resolved before manuscript submission.
- `data/author_reported/`: transcriptions of supplied manuscript tables, explicitly unverified as computational results.
- `data/reference/`: source metadata from the Allen Institute, with provenance.
- `tests/`: computational checks using clearly identified synthetic fixtures, not research data.

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

Analysis code, documentation, source donor metadata, and clearly labeled transcriptions of author-reported tables are available in this repository. Original transcriptomic data are available through the [Allen Human Brain Atlas](https://human.brain-map.org/). The complete original processed expression dataset and original executed analysis outputs are not currently included. No claim of complete reproduction or full FAIR compliance is made.

## Sources and license

- AHBA: Hawrylycz et al. 2012, https://doi.org/10.1038/nature11405
- abagen: Markello et al. 2021, https://doi.org/10.7554/eLife.72129
- g:Profiler: Raudvere et al. 2019, https://doi.org/10.1093/nar/gkz369
- Imaging source discussed in the manuscript: Hervais-Adelman et al. 2015, https://doi.org/10.1093/cercor/bhu158

Repository code is covered by the existing MIT `License`. Upstream datasets and services remain subject to their own terms; the repository license does not relicense AHBA data. Cite the original atlas and methods alongside this repository. The research manuscript has not been assigned a publication DOI here.
