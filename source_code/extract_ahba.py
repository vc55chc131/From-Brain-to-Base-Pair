#!/usr/bin/env python3
"""Prepare or run a NEW, explicitly specified AHBA extraction.

This does not recover the missing input data or results of the original paper.
No downloads occur in --prepare-only mode. Full execution can download all six
AHBA microarray archives through abagen and requires substantial RAM/disk space.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import sys

DONORS = ["9861", "10021", "12876", "14380", "15496", "15697"]
DONOR_NAMES = dict(zip(DONORS, ["H0351.2001", "H0351.2002", "H0351.1009",
                              "H0351.1012", "H0351.1015", "H0351.1016"]))
NORMALIZATION = (
    "AHBA microarray data processed with abagen: intensity-based filtering at "
    "0.5; maximum-intensity representative probe selected across donors; "
    "reannotated probes; corrected MNI coordinates; scaled robust sigmoid "
    "normalization across genes within samples and across all retained samples "
    "within each donor for each gene, without separate structural-class normalization; "
    "arithmetic tissue means within each donor and ROI; Caudate hemispheres "
    "combined within each donor using observed tissue-sample counts as weights. "
    "Values are normalized expression units, not log2 fold changes or z scores."
)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_seeds(path, accept_unverified=False):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    required = {"id", "roi", "hemisphere", "structure", "x", "y", "z",
                "radius_mm", "provenance_status", "source"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Seed CSV must contain: " + ", ".join(sorted(required)))
    for row in rows:
        if row["provenance_status"] != "reviewed" and not accept_unverified:
            raise ValueError("Seed provenance is unverified. Review the source coordinates "
                             "or explicitly use --accept-unverified-seeds for exploratory preparation.")
        if row["provenance_status"] == "reviewed" and not row["source"].strip():
            raise ValueError("Reviewed seeds require a source and coordinate location.")
        for key in ("x", "y", "z", "radius_mm"):
            row[key] = float(row[key])
        row["id"] = int(row["id"])
        if row["id"] < 1 or not 0 < row["radius_mm"] <= 20:
            raise ValueError("Seed IDs must be positive and radii must be in (0, 20] mm.")
        if row["roi"] not in ("DLPFC", "Caudate") or row["hemisphere"] not in ("L", "R"):
            raise ValueError("ROI must be DLPFC or Caudate; hemisphere must be L or R.")
        if (row["hemisphere"] == "L" and row["x"] >= 0) or (row["hemisphere"] == "R" and row["x"] <= 0):
            raise ValueError("MNI x coordinate must agree with the hemisphere designation.")
        expected = "cortex" if row["roi"] == "DLPFC" else "subcortex/brainstem"
        if row["structure"] != expected:
            raise ValueError("Structural class must be cortex for DLPFC and subcortex/brainstem for Caudate.")
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("Seed IDs must be unique.")
    if {(r["roi"], r["hemisphere"]) for r in rows} != {("DLPFC", "L"), ("Caudate", "L"), ("Caudate", "R")} or len(rows) != 3:
        raise ValueError("Exactly one left DLPFC and one Caudate seed per hemisphere are required.")
    return rows


def make_atlas(seeds, np, pd, nib, load_mni152_template):
    # This template is distributed with nilearn; its affine defines MNI space.
    template = load_mni152_template(resolution=2)
    voxels = np.indices(template.shape).reshape(3, -1).T
    world = nib.affines.apply_affine(template.affine, voxels)
    atlas_flat = np.zeros(len(world), dtype=np.int16)
    info = []
    for row in seeds:
        centre = np.array([row["x"], row["y"], row["z"]])
        selected = np.linalg.norm(world - centre, axis=1) <= row["radius_mm"]
        if not selected.any():
            raise ValueError(f"Seed {row['id']} contains no template voxel centres.")
        if (atlas_flat[selected] != 0).any():
            raise ValueError("Seed spheres overlap; revise the seeds before extraction.")
        atlas_flat[selected] = row["id"]
        info.append({**row, "voxel_count": int(selected.sum())})
    atlas = nib.Nifti1Image(atlas_flat.reshape(template.shape), template.affine)
    return atlas, pd.DataFrame(info)


def export_expression(expression, counts, seeds, out, np, pd):
    """Convert actual abagen donor-region means; never label them tissue IDs."""
    if not isinstance(expression, dict):
        raise ValueError("Expected an abagen donor-keyed dictionary; incompatible API output.")
    counts = counts.copy()
    counts.columns = counts.columns.map(str)
    counts.index = counts.index.astype(int)
    counts.to_csv(out / "roi_tissue_counts.csv", index_label="atlas_roi_id")
    region_dir = out / "donor_region_expression"
    region_dir.mkdir()
    columns, metadata = {}, []
    for donor, frame in expression.items():
        donor = str(donor)
        frame = frame.copy()
        frame.index = frame.index.astype(int)
        frame.to_csv(region_dir / f"{donor}.csv", index_label="atlas_roi_id")
        for roi in ("DLPFC", "Caudate"):
            used, weights, hemispheres = [], [], []
            for seed in (s for s in seeds if s["roi"] == roi):
                n = int(counts.loc[seed["id"], donor])
                if n > 0:
                    used.append(frame.loc[seed["id"]])
                    weights.append(n)
                    hemispheres.append(seed["hemisphere"])
            if not used:
                continue
            values = np.average(np.vstack(used), axis=0, weights=weights)
            sample_id = f"donor_{donor}_{roi}"
            columns[sample_id] = pd.Series(values, index=used[0].index)
            metadata.append({"sample_id": sample_id, "donor_id": DONOR_NAMES[donor],
                             "roi": roi, "hemisphere": "B" if len(hemispheres) == 2 else hemispheres[0],
                             "input_level": "donor_roi_summary", "n_tissue_samples": sum(weights),
                             "hemispheres_observed": ";".join(hemispheres)})
    if not columns:
        raise ValueError("No tissue samples matched these seeds. Counts were saved; no expression dataset was invented.")
    matrix = pd.DataFrame(columns)
    kept = np.isfinite(matrix.to_numpy()).all(axis=1)
    gene_audit = pd.DataFrame({"gene_id": matrix.index, "identifier_type": "gene_symbol", "retained_complete": kept})
    gene_audit.to_csv(out / "gene_id_audit.csv", index=False)
    matrix = matrix.loc[kept]
    if matrix.empty:
        raise ValueError("No genes have finite values across the observed donor-region summaries.")
    meta = pd.DataFrame(metadata)
    matrix.to_csv(out / "expression.csv", index_label="gene_id")
    meta.to_csv(out / "sample_metadata.csv", index=False)
    pairs = meta.groupby("donor_id")["roi"].nunique()
    paired = pairs.index[pairs == 2].tolist()
    paired_meta = meta.loc[meta["donor_id"].isin(paired)].copy()
    if paired:
        matrix.loc[:, paired_meta["sample_id"]].to_csv(out / "paired_expression.csv", index_label="gene_id")
        paired_meta.to_csv(out / "paired_sample_metadata.csv", index=False)
    excluded = []
    for donor in expression:
        donor_name = DONOR_NAMES[str(donor)]
        if donor_name not in paired:
            observed = sorted(meta.loc[meta["donor_id"] == donor_name, "roi"].tolist())
            excluded.append({"donor_id": donor_name, "observed_rois": observed,
                             "reason": "incomplete_ROI_pair" if observed else "no_ROI_coverage"})
    diagnostics = {"n_retained_genes": len(matrix), "n_donor_region_summaries": len(metadata),
                   "n_paired_donors": len(paired), "paired_donors": paired,
                   "paired_analysis_exclusions": excluded,
                   "paired_exports_created": bool(paired),
                   "gene_filter_policy": "Finite expression in all observed donor-region summaries before complete-donor subsetting; no subsequent gene refiltering.",
                   "note": "Donor-region summary IDs are generated identifiers, not Allen tissue sample IDs."}
    (out / "coverage_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    return diagnostics


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=Path("data/ahba_cache"))
    parser.add_argument("--accept-unverified-seeds", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--tolerance-mm", type=int, choices=range(0, 11), default=0,
                        help="Extra abagen voxel-matching tolerance; default 0, distinct from sphere radius.")
    args = parser.parse_args(argv)
    try:
        seeds = read_seeds(args.seeds, args.accept_unverified_seeds)
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            raise ValueError("Output directory must be empty to avoid mixing runs.")
        try:
            import numpy as np
            import pandas as pd
            import nibabel as nib
            import abagen
            from nilearn.datasets import load_mni152_template
        except ImportError as error:
            raise ValueError("Missing AHBA extraction dependency. Install requirements-ahba.txt in a separate environment. " + str(error)) from error
        atlas, info = make_atlas(seeds, np, pd, nib, load_mni152_template)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        atlas_file = args.output_dir / "seed_atlas.nii.gz"
        nib.save(atlas, atlas_file)
        info.to_csv(args.output_dir / "atlas_info.csv", index=False)
        parameters = dict(ibf_threshold=0.5, probe_selection="max_intensity", donor_probes="aggregate",
                          sim_threshold=None, lr_mirror=None, missing=None, tolerance=args.tolerance_mm,
                          sample_norm="srs", gene_norm="srs", norm_matched=False, norm_structures=False,
                          region_agg="donors", agg_metric="mean", corrected_mni=True, reannotated=True,
                          return_donors=True, return_counts=True, return_report=True, donors=DONORS,
                          data_dir=str(args.data_dir.resolve()), verbose=1, n_proc=1)
        manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                    "status": "prepared_only" if args.prepare_only else "started",
                    "purpose": "new_reanalysis_not_recovery_of_original_results",
                    "seed_sha256": sha256(args.seeds), "atlas_sha256": sha256(atlas_file),
                    "parameters": parameters, "seeds": seeds,
                    "seed_provenance_unverified": any(s["provenance_status"] != "reviewed" for s in seeds),
                    "versions": {p: importlib.metadata.version(p) for p in ("abagen", "nilearn", "nibabel", "numpy", "pandas", "scipy")},
                    "abagen_reported_version": abagen.__version__,
                    "abagen_workflow_source_sha256": sha256(inspect.getsourcefile(abagen.get_expression_data)),
                    "requirements_abagen_commit": "dc4a007e4e902e51f97251390c8d1bbf7e58c6d3",
                    "normalization": NORMALIZATION,
                    "sample_assignment": "Spherical atlas in MNI152 space, hemisphere and broad structural-class constraints; no mirroring or imputation."}
        manifest_file = args.output_dir / "extraction_manifest.json"
        manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
        (args.output_dir / "normalization_description.txt").write_text(NORMALIZATION + "\n")
        if args.prepare_only:
            print("Atlas and configuration prepared. AHBA expression data were not downloaded or analysed.")
            return 0
        try:
            expression, counts, report = abagen.get_expression_data(atlas, info, **parameters)
            (args.output_dir / "abagen_methods.txt").write_text(report + "\n")
            # The cached raw files identify the exact public inputs used. This
            # performs no new downloads when the preceding extraction succeeded.
            files = abagen.fetch_microarray(data_dir=str(args.data_dir.resolve()),
                                            donors=DONORS, convert=False, verbose=0)
            source_rows = []
            for donor, donor_files in files.items():
                for role, filename in donor_files.items():
                    source_rows.append({"donor_id": DONOR_NAMES[str(donor)], "role": role,
                                        "filename": Path(filename).name,
                                        "size_bytes": Path(filename).stat().st_size,
                                        "sha256": sha256(filename)})
            pd.DataFrame(source_rows).to_csv(args.output_dir / "source_file_manifest.csv", index=False)
            diagnostics = export_expression(expression, counts, seeds, args.output_dir, np, pd)
            manifest.update(status="completed_extraction", diagnostics=diagnostics,
                            sufficient_paired_coverage=diagnostics["n_paired_donors"] >= 3)
        except Exception as error:
            manifest.update(status="failed", error=str(error))
            manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
            raise
        manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
        print(json.dumps(diagnostics, indent=2))
        if diagnostics["n_paired_donors"] < 3:
            print("Insufficient paired donors for the companion inferential analysis. Do not expand seeds to obtain significance.", file=sys.stderr)
            return 2
        return 0
    except (ValueError, OSError) as error:
        parser.exit(2, "Error: " + str(error) + "\n")


if __name__ == "__main__":
    sys.exit(main())
