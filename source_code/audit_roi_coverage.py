#!/usr/bin/env python3
"""Audit proposed ROI coverage using public spatial metadata only.

No gene-expression files are loaded. Exact coordinate spheres are evaluated
separately from abagen's voxel-atlas assignment at tolerances 0 and 2 mm.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
import platform

import abagen
import nibabel as nib
from nilearn.datasets import load_mni152_template
import nilearn
import numpy as np
import pandas as pd
import scipy

DONORS = {"9861": "H0351.2001", "10021": "H0351.2002", "12876": "H0351.1009",
          "14380": "H0351.1012", "15496": "H0351.1015", "15697": "H0351.1016"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_atlas(seeds):
    template = load_mni152_template(resolution=2)
    voxels = np.indices(template.shape).reshape(3, -1).T
    coords = nib.affines.apply_affine(template.affine, voxels)
    labels = np.zeros(len(coords), dtype=np.int16)
    info = []
    for row in seeds.to_dict("records"):
        within = np.linalg.norm(coords - np.array([row["x"], row["y"], row["z"]]), axis=1) <= row["radius_mm"]
        if not within.any() or np.any(labels[within]):
            raise ValueError("Empty or overlapping voxelized ROI")
        labels[within] = int(row["id"])
        info.append({**row, "n_atlas_voxels": int(within.sum())})
    return nib.Nifti1Image(labels.reshape(template.shape), template.affine), pd.DataFrame(info)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path("data/feasibility/inputs"))
    parser.add_argument("--seeds", type=Path, default=Path("config/roi_seeds_unverified.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/feasibility/results"))
    args = parser.parse_args()
    out = args.output_dir
    if out.exists() and any(out.iterdir()):
        parser.error("Output directory must be empty to preserve run provenance")
    out.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(args.seeds)
    required = {"id", "roi", "hemisphere", "structure", "x", "y", "z", "radius_mm"}
    if not required.issubset(seeds) or seeds["id"].duplicated().any():
        parser.error("Missing seed fields or duplicate ROI IDs")
    atlas, atlas_info = build_atlas(seeds)
    nib.save(atlas, out / "seed_atlas.nii.gz")
    atlas_info.to_csv(out / "atlas_info.csv", index=False)
    trees = {"unconstrained": abagen.images.check_atlas(atlas),
             "hemisphere_and_structure": abagen.images.check_atlas(atlas, atlas_info)}
    corrected_path = args.inputs / "corrected_mni_coordinates.csv.gz"
    corrected = pd.read_csv(corrected_path).set_index("well_id", verify_integrity=True)
    response = json.loads((args.inputs / "well_donor_mapping.json").read_text())
    if not response["success"] or response["num_rows"] != response["total_rows"]:
        raise ValueError("Incomplete or unsuccessful donor mapping response")
    wells = pd.DataFrame(response["msg"]).set_index("id", verify_integrity=True)
    matched = wells.loc[corrected.index]
    input_files = [corrected_path, args.inputs / "well_donor_mapping.json", args.seeds]
    counts, nearest, assignments, sample_totals = [], [], [], []
    for donor_api_id, donor_id in DONORS.items():
        donor_dir = args.inputs / ("donor_" + donor_api_id)
        annotation_path, ontology_path = donor_dir / "SampleAnnot.csv", args.inputs / "Ontology.csv"
        input_files.extend([annotation_path, ontology_path])
        original = abagen.io.read_annotation(str(annotation_path))
        if original["well_id"].duplicated().any():
            raise ValueError("Duplicate well IDs within donor annotations")
        expected_wells = matched.index[matched["donor_id"] == int(donor_api_id)]
        if set(original["well_id"]) != set(expected_wells):
            raise ValueError("Corrected coordinate, annotation and API donor mapping disagree")
        modes = {"original_mni": original.copy(), "corrected_mni": original.copy()}
        modes["corrected_mni"][["mni_x", "mni_y", "mni_z"]] = corrected.loc[original["well_id"]].to_numpy()
        for mode, annotation in modes.items():
            xyz = annotation[["mni_x", "mni_y", "mni_z"]].to_numpy()
            if not np.isfinite(xyz).all():
                raise ValueError("Nonfinite sample coordinates")
            constrained = abagen.samples_.drop_mismatch_samples(annotation, str(ontology_path))
            sample_totals.append({"coordinate_set": mode, "donor_id": donor_id, "donor_api_id": donor_api_id,
                                  "n_samples": len(annotation), "n_after_hemisphere_mismatch_filter": len(constrained)})
            for seed in seeds.to_dict("records"):
                distances = np.linalg.norm(xyz - np.array([seed["x"], seed["y"], seed["z"]]), axis=1)
                i = int(np.argmin(distances))
                sample = annotation.iloc[i]
                nearest.append({"coordinate_set": mode, "donor_id": donor_id, "donor_api_id": donor_api_id,
                                "roi_id": seed["id"], "roi": seed["roi"], "hemisphere": seed["hemisphere"],
                                "seed_x": seed["x"], "seed_y": seed["y"], "seed_z": seed["z"],
                                "well_id": int(sample["well_id"]), "sample_x": sample["mni_x"],
                                "sample_y": sample["mni_y"], "sample_z": sample["mni_z"],
                                "distance_to_seed_mm": float(distances[i]), "filtering": "unconstrained"})
            for filtering, samples in [("unconstrained", annotation), ("hemisphere_and_structure", constrained)]:
                sample_xyz = samples[["mni_x", "mni_y", "mni_z"]].to_numpy()
                exact = np.zeros(len(samples), dtype=int)
                for seed in seeds.to_dict("records"):
                    keep = np.linalg.norm(sample_xyz - np.array([seed["x"], seed["y"], seed["z"]]), axis=1) <= seed["radius_mm"]
                    if filtering != "unconstrained":
                        keep &= (samples["hemisphere"].to_numpy() == seed["hemisphere"])
                        keep &= (samples["structure"].to_numpy() == seed["structure"])
                    exact[keep] = int(seed["id"])
                outputs = [("exact_euclidean_sphere", "not_applicable", exact)]
                for tolerance in (0, 2):
                    labels = trees[filtering].label_samples(samples, tolerance=tolerance)["label"].to_numpy()
                    outputs.append(("abagen_voxel_atlas", tolerance, labels))
                for method, tolerance, labels in outputs:
                    row = {"coordinate_set": mode, "filtering": filtering, "matching_method": method,
                           "tolerance_mm": tolerance, "donor_id": donor_id, "donor_api_id": donor_api_id,
                           "n_samples_evaluated": len(samples)}
                    for seed in seeds.to_dict("records"):
                        row[f"n_{seed['roi']}_{seed['hemisphere']}"] = int((labels == seed["id"]).sum())
                    row["has_DLPFC_Caudate_pair"] = bool(row["n_DLPFC_L"] > 0 and row["n_Caudate_L"] + row["n_Caudate_R"] > 0)
                    counts.append(row)
                    for idx in np.flatnonzero(labels):
                        sample = samples.iloc[idx]
                        assignments.append({**{k: row[k] for k in ("coordinate_set", "filtering", "matching_method", "tolerance_mm", "donor_id")},
                                            "well_id": int(sample["well_id"]), "roi_id": int(labels[idx]),
                                            "sample_structure": sample.get("structure_acronym", ""),
                                            "mni_x": sample["mni_x"], "mni_y": sample["mni_y"], "mni_z": sample["mni_z"]})
    counts = pd.DataFrame(counts)
    counts.to_csv(out / "donor_counts.csv", index=False)
    pd.DataFrame(nearest).to_csv(out / "nearest_samples.csv", index=False)
    pd.DataFrame(sample_totals).to_csv(out / "sample_totals.csv", index=False)
    pd.DataFrame(assignments).to_csv(out / "matched_samples.csv", index=False)
    groups = ["coordinate_set", "filtering", "matching_method", "tolerance_mm"]
    summary = counts.groupby(groups, dropna=False, sort=False).agg(
        n_DLPFC_L=("n_DLPFC_L", "sum"), n_Caudate_L=("n_Caudate_L", "sum"), n_Caudate_R=("n_Caudate_R", "sum"),
        n_paired_donors=("has_DLPFC_Caudate_pair", "sum")).reset_index()
    summary.to_csv(out / "coverage_summary.csv", index=False)
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "analysis": "spatial_coverage_only_no_gene_expression",
                "n_corrected_coordinate_wells": len(corrected), "n_donors": len(DONORS),
                "coordinate_source": "Pinned abagen corrected MNI table derived from alleninf; original MNI coordinates from Allen SampleAnnot.csv",
                "donor_mapping": "All corrected-coordinate well IDs uniquely joined to saved Allen Well API donor_id records and exact per-donor SampleAnnot.csv well-ID sets",
                "seed_provenance": "Proposed manuscript coordinates; no claim that these coordinates reproduce a published fMRI contrast",
                "atlas": "2 mm nilearn MNI152 template; labelled voxels have centres within the proposed 2 mm coordinate spheres",
                "matching_note": "Exact continuous Euclidean spheres and abagen voxel-atlas matching are distinct procedures. abagen first projects coordinates onto its voxel grid before incremental tolerance matching. Tolerance 0 is not a continuous 2 mm distance test.",
                "filtering_note": "Unconstrained counts are spatial upper bounds. Constrained counts additionally remove hemisphere-coordinate mismatches and constrain assignments by hemisphere and broad structural class using Allen ontology metadata.",
                "scope_limit": "This audit cannot test gene expression differences, enrichment, or interpreter-specific biology. It evaluates only these coordinates, radii, coordinate sets and matching rules.",
                "software": {"python": platform.python_version(), "abagen_module_version": abagen.__version__,
                             "abagen_pinned_commit": "dc4a007e4e902e51f97251390c8d1bbf7e58c6d3", "nilearn": nilearn.__version__,
                             "nibabel": nib.__version__, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
                "source_sha256": {"audit_roi_coverage.py": sha256(__file__), "abagen_matching.py": sha256(inspect.getsourcefile(abagen.AtlasTree)),
                                  "abagen_samples.py": sha256(inspect.getsourcefile(abagen.samples_.drop_mismatch_samples))},
                "input_sha256": {str(p.relative_to(args.inputs)) if p.is_relative_to(args.inputs) else str(p): sha256(p) for p in input_files},
                "output_sha256": {p.name: sha256(p) for p in out.iterdir() if p.is_file()}}
    (out / "audit_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
