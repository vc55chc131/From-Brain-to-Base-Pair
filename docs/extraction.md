# Optional extraction from the Allen Human Brain Atlas

This is a newly specified reanalysis pipeline. It does not recover the original study dataset or establish that the original figures, gene lists, statistical tests or enrichment results are reproducible. Full microarray extraction has not been executed or validated. No scientific results are supplied by the preparation check. The optional environment pins abagen to upstream commit `dc4a007e4e902e51f97251390c8d1bbf7e58c6d3`, which includes pandas 2 compatibility fixes missing from the PyPI 0.1.3 release.

The source microarrays are available from the [Allen Human Brain Atlas](https://human.brain-map.org/static/download). The script calls the documented [abagen expression workflow](https://abagen.readthedocs.io/en/stable/generated/abagen.get_expression_data.html). Read the source terms before redistributing Allen data; this repository's code license does not replace them.

## Coordinate provenance and preparation

`config/roi_seeds_unverified.csv` records the coordinates stated in the supplied manuscript. Their correspondence to the cited primary fMRI study has not been established. The file marks this explicitly. The extractor stops by default for unverified seeds. To create the atlas and inspect the configuration without downloading expression data:

```bash
python -m venv .venv-ahba
source .venv-ahba/bin/activate
python -m pip install -r requirements-ahba.txt
python source_code/extract_ahba.py --seeds config/roi_seeds_unverified.csv --output-dir outputs/atlas_preparation --accept-unverified-seeds --prepare-only
```

The unverified-seed flag acknowledges an unresolved issue for exploratory work; it does not validate the anatomical labels or the source attribution. Before a publishable reanalysis, check the original coordinates, contrast, participant group and coordinate space. Save corrected/reviewed seeds to a new CSV, set their provenance status to `reviewed` and record the precise source table or location. Anatomical checks should also confirm that the spheres lie in the intended structures.

## Full extraction after coordinate review

```bash
python source_code/extract_ahba.py --seeds config/roi_seeds_reviewed.csv --data-dir data/ahba_cache --output-dir outputs/ahba_extraction
```

The reviewed CSV in this command is an author-created input, not a file supplied here. Full execution may download all six large microarray archives and require substantial memory and disk space. Use an empty output directory for each run.

The pipeline uses a 2 mm resolution MNI152 template distributed with nilearn. Sphere radii come from the seed CSV. The separate `--tolerance-mm` setting defaults to zero, so it does not expand the atlas matching search. Setting it to 2 adds up to 2 mm of abagen matching tolerance beyond the labelled atlas voxels; that is not equivalent to a 2 mm radius criterion around the coordinate centre. Record and justify any change before examining results.

The parameters are explicit: intensity-based probe filtering at 0.5, maximum-intensity probe selection across donors, reannotated probes, corrected MNI coordinates, scaled robust sigmoid normalization within samples and then within genes across all retained donor samples. Normalization is not restricted to the three small ROIs and is not separated by structural class. Tissue matching is constrained by hemisphere and broad structure. Missing regions remain missing. No mirroring, interpolation, nearest-centroid filling or inferred tissue IDs are used.

## Outputs and companion analysis

The extractor writes the atlas and its voxel counts, explicit parameters and versions, a checksum for the installed abagen workflow source, a manifest of the downloaded source files and their checksums, the abagen methods report, exact tissue counts for each donor and atlas region, per-donor region expression matrices, and a gene-ID inclusion audit. The gene identifiers are the gene symbols returned by abagen. They are preserved verbatim. This high-level interface does not expose the final selected probe IDs; the audit is not a probe-selection map. An archive installation can report the abagen version as `0+unknown`; the pinned upstream commit and workflow source checksum establish the intended version more precisely.

`expression.csv` has genes in rows with first column `gene_id`. `sample_metadata.csv` identifies each generated donor-region summary with `sample_id`, `donor_id`, `roi`, `hemisphere`, `input_level`, and `n_tissue_samples`. These generated sample IDs are explicitly not Allen tissue sample IDs. ROIs are `DLPFC` and `Caudate`. Separate hemisphere means are preserved in `donor_region_expression/`. The combined Caudate summary weights hemisphere means by their actual matched tissue-sample counts within the same donor. A single observed hemisphere is retained and marked L or R; B means both were observed. This does not establish bilateral coverage in all donors.

No donor or gene count is predetermined. Genes with nonfinite expression in any observed donor-region summary are excluded and logged before complete donors are selected; there is no further gene filtering or reintroduction after this selection. Absent donor regions do not become zero-valued observations. The full observed outputs remain in `expression.csv` and `sample_metadata.csv`. For the companion analysis, `paired_expression.csv` and `paired_sample_metadata.csv` contain only donors with both DLPFC and Caudate observations. Coverage diagnostics report the actual paired donors and identify every excluded donor, their observed ROIs and the reason for exclusion. If no donor has both regions, paired files are not created. The script returns an unsuccessful status when fewer than three paired donors remain, without fabricating data or automatically enlarging ROIs. Coverage failure is a reason to reconsider the design and anatomical specification.

After successful extraction with at least three paired donors, run the companion analysis on the paired exports. The command below reads the actual normalization description without shell interpolation:

```python
from pathlib import Path
import subprocess
import sys

inputs = Path("outputs/ahba_extraction")
subprocess.run([
    sys.executable, "source_code/DLPFC_vs_Caudate_Zscore.py",
    "--expression", str(inputs / "paired_expression.csv"),
    "--metadata", str(inputs / "paired_sample_metadata.csv"),
    "--output-dir", "outputs/paired_reanalysis",
    "--normalization-description", (inputs / "normalization_description.txt").read_text(),
], check=True)
```

Expression values are normalized units, not log2 fold changes, z scores, or absolute transcript counts. Compare regions at the donor level and report the number of paired donors retained for each analysis. Full data execution, anatomical validation and evaluation of analytical sensitivity remain to be completed before scientific conclusions can be drawn.

## Verified donor reference records

`data/reference/allen_donors_api.json` is the retrieved Allen API response. `allen_donors.csv` transcribes donor ID, age, sex, postmortem interval and handedness. `provenance.json` contains the URL, access timestamp and response checksum. These are whole-donor reference records, not evidence of ROI sampling coverage. No cause of death is inferred. The current API reports PMI of 26 hours for H0351.1009 and 17 hours for H0351.1012; abagen's packaged donor table lists 25.5 and 17.5 hours respectively. The CSV follows the saved current Allen API response and preserves this discrepancy in its provenance note.
