#!/usr/bin/env python3
"""New donor-paired analysis implementation, not the unavailable original pipeline.

Input is an already processed gene-by-tissue-sample CSV, with genes in the first
column, plus sample metadata containing sample_id, donor_id, roi, and hemisphere.
Preprocessing and ROI assignment must be completed and documented upstream.
The legacy filename is retained for discoverability. This script DOES NOT compute
log2 fold changes or an across-gene z-score. Its effect is the within-donor DLPFC
minus Caudate difference on the supplied normalized expression scale.

All donors must contribute both ROIs. Tissue samples are averaged within each
donor and ROI before inference, so a donor with more tissue samples receives no
additional weight. Three donors is a minimum execution guard, not a claim of
adequate power. With the small AHBA donor sample, inference remains exploratory.

Primary API references:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html
https://biit.cs.ut.ee/gprofiler/page/apis
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from scipy import stats


IMPLEMENTATION_VERSION = "1.0.0-new-reanalysis"
REGIONS = ("DLPFC", "Caudate")
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
GO_SOURCES = ["GO:BP", "GO:MF", "GO:CC"]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)
                    + "\n", encoding="utf-8")


def read_csv_checked(path: Path) -> pd.DataFrame:
    """Reject duplicate headers before pandas can silently rename them."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle), None)
    if not header or any(not field.strip() for field in header):
        raise ValueError(f"{path.name}: missing or blank column headers")
    if len(header) != len(set(header)):
        raise ValueError(f"{path.name}: duplicate column headers")
    if any(field != field.strip() for field in header):
        raise ValueError(f"{path.name}: column headers contain surrounding whitespace")
    # Read identifiers as strings, preserving numeric IDs and leading zeroes.
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")


def require_ids(series: pd.Series, description: str, unique: bool = True) -> None:
    if series.empty or series.eq("").any() or series.str.strip().ne(series).any():
        raise ValueError(f"{description}: identifiers must be nonempty and trimmed")
    if unique and series.duplicated().any():
        raise ValueError(f"{description}: duplicate identifiers are not allowed")


def load_inputs(expression_path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = read_csv_checked(expression_path)
    if raw.shape[1] < 2 or raw.empty:
        raise ValueError("Expression matrix needs at least one gene and one sample")
    require_ids(raw.iloc[:, 0], "Gene identifiers")
    expression = raw.set_index(raw.columns[0])
    expression.index.name = "gene_id"
    try:
        expression = expression.apply(pd.to_numeric, errors="raise").astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Expression values must all be numeric") from exc
    if not np.isfinite(expression.to_numpy()).all():
        raise ValueError("Expression values must all be finite, without missing values")

    metadata = read_csv_checked(metadata_path)
    required = {"sample_id", "donor_id", "roi", "hemisphere"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Metadata is missing columns: {sorted(required - set(metadata.columns))}")
    require_ids(metadata["sample_id"], "Sample identifiers")
    require_ids(metadata["donor_id"], "Donor identifiers", unique=False)
    if not metadata["roi"].isin(REGIONS).all():
        raise ValueError("Metadata roi must be exactly DLPFC or Caudate")
    if "input_level" not in metadata:
        metadata["input_level"] = "tissue_sample"
    input_levels = metadata["input_level"].unique().tolist()
    if len(input_levels) != 1 or input_levels[0] not in ("tissue_sample", "donor_roi_summary"):
        raise ValueError("input_level must consistently be tissue_sample or donor_roi_summary")
    if input_levels[0] == "tissue_sample":
        if not metadata["hemisphere"].isin(["L", "R"]).all():
            raise ValueError("Tissue sample hemisphere must be L or R")
        if "n_tissue_samples" in metadata and not metadata["n_tissue_samples"].eq("1").all():
            raise ValueError("Each tissue_sample input row must represent exactly one tissue sample")
        metadata["n_tissue_samples"] = 1
    else:
        if not metadata["hemisphere"].isin(["L", "R", "B"]).all():
            raise ValueError("Summary hemisphere must be L, R, or B for bilateral")
        if metadata.duplicated(["donor_id", "roi"]).any():
            raise ValueError("donor_roi_summary requires exactly one input row per donor and ROI")
        if "n_tissue_samples" not in metadata:
            raise ValueError("donor_roi_summary requires observed n_tissue_samples")
        try:
            counts = pd.to_numeric(metadata["n_tissue_samples"], errors="raise")
        except ValueError as exc:
            raise ValueError("n_tissue_samples must be positive integers") from exc
        if (not np.isfinite(counts).all() or (counts < 1).any()
                or not np.equal(counts, np.floor(counts)).all()):
            raise ValueError("n_tissue_samples must be positive integers")
        metadata["n_tissue_samples"] = counts.astype(int)
    expression_samples = set(expression.columns)
    metadata_samples = set(metadata["sample_id"])
    if expression_samples != metadata_samples:
        raise ValueError(
            "Expression and metadata sample IDs must match exactly; "
            f"without metadata={sorted(expression_samples - metadata_samples)}; "
            f"without expression={sorted(metadata_samples - expression_samples)}"
        )
    metadata = metadata.set_index("sample_id").loc[expression.columns].rename_axis("sample_id").reset_index()
    pairing = metadata.groupby("donor_id")["roi"].agg(set)
    incomplete = [str(donor) for donor, regions in pairing.items()
                  if regions != set(REGIONS)]
    if incomplete:
        raise ValueError(f"Every donor must contribute both ROIs; incomplete donors={incomplete}")
    if len(pairing) < 3:
        raise ValueError("At least three independent donors with both ROIs are required")
    return expression, metadata


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; missing tests remain missing."""
    pvalues = np.asarray(pvalues, dtype=float)
    finite = np.isfinite(pvalues)
    if np.any((pvalues[finite] < 0) | (pvalues[finite] > 1)):
        raise ValueError("Finite p-values must be within [0, 1]")
    if np.isinf(pvalues).any():
        raise ValueError("P-values may not be infinite")
    result = np.full(pvalues.shape, np.nan, dtype=float)
    positions = np.flatnonzero(finite)
    if not positions.size:
        return result
    order = np.argsort(pvalues[positions], kind="stable")
    sorted_positions = positions[order]
    n_tests = len(sorted_positions)
    adjusted = pvalues[sorted_positions] * n_tests / np.arange(1, n_tests + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[sorted_positions] = np.clip(adjusted, 0, 1)
    return result


def analyze(expression: pd.DataFrame, metadata: pd.DataFrame,
            alpha: float = 0.05) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return all-gene statistics, donor means, counts, and donor differences."""
    donors = sorted(metadata["donor_id"].unique())
    means: dict[tuple[str, str], pd.Series] = {}
    for donor in donors:
        for roi in REGIONS:
            sample_ids = metadata.loc[
                (metadata["donor_id"] == donor) & (metadata["roi"] == roi), "sample_id"
            ]
            if sample_ids.empty:
                raise ValueError(f"Missing {roi} for donor {donor}")
            with np.errstate(over="ignore", invalid="ignore"):
                means[(donor, roi)] = expression.loc[:, sample_ids].mean(axis=1)
            if not np.isfinite(means[(donor, roi)].to_numpy()).all():
                raise ValueError("Nonfinite donor-ROI means after aggregation; check input scale for overflow")
    donor_means = pd.DataFrame(means)
    donor_means.columns = pd.MultiIndex.from_tuples(donor_means.columns,
                                                   names=["donor_id", "roi"])
    dlpfc = np.column_stack([means[(donor, "DLPFC")] for donor in donors])
    caudate = np.column_stack([means[(donor, "Caudate")] for donor in donors])
    with np.errstate(over="ignore", invalid="ignore"):
        difference = dlpfc - caudate
    if not np.isfinite(difference).all():
        raise ValueError("Nonfinite donor differences after subtraction; check input scale for overflow")
    with np.errstate(over="ignore", invalid="ignore"):
        effect = difference.mean(axis=1)
        standard_deviation = difference.std(axis=1, ddof=1)
        dlpfc_mean = dlpfc.mean(axis=1)
        caudate_mean = caudate.mean(axis=1)
    if not all(np.isfinite(values).all() for values in
               [effect, standard_deviation, dlpfc_mean, caudate_mean]):
        raise ValueError("Nonfinite effect or moments; check input scale for numerical overflow")
    # A t-test is undefined when the donor differences have zero variance.
    # Use a roundoff-scale tolerance, without converting it into a p=0 result.
    tolerance = np.finfo(float).eps * np.maximum(1.0, np.abs(difference).max(axis=1)) * 10
    valid = standard_deviation > tolerance
    tstat = np.full(len(expression), np.nan)
    pvalue = np.full(len(expression), np.nan)
    ci_low = np.full(len(expression), np.nan)
    ci_high = np.full(len(expression), np.nan)
    if valid.any():
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            paired = stats.ttest_rel(dlpfc[valid], caudate[valid], axis=1,
                                    nan_policy="raise", alternative="two-sided")
        tstat[valid] = paired.statistic
        pvalue[valid] = paired.pvalue
        critical = stats.t.ppf(0.975, df=len(donors) - 1)
        with np.errstate(over="ignore", invalid="ignore"):
            halfwidth = critical * standard_deviation[valid] / np.sqrt(len(donors))
            ci_low[valid] = effect[valid] - halfwidth
            ci_high[valid] = effect[valid] + halfwidth
        if not all(np.isfinite(values[valid]).all() for values in
                   [tstat, pvalue, ci_low, ci_high]):
            raise ValueError("Nonfinite test or confidence interval; check input scale for numerical overflow")
    qvalue = bh_adjust(pvalue)
    result = pd.DataFrame({
        "gene_id": expression.index,
        "n_paired_donors": len(donors),
        "DLPFC_mean_across_donors": dlpfc_mean,
        "Caudate_mean_across_donors": caudate_mean,
        "mean_paired_difference_DLPFC_minus_Caudate": effect,
        "sd_paired_difference": standard_deviation,
        "ci95_low_unadjusted": ci_low,
        "ci95_high_unadjusted": ci_high,
        "t_statistic": tstat,
        "degrees_of_freedom": len(donors) - 1,
        "p_value_two_sided": pvalue,
        "q_value_BH": qvalue,
        "test_status": np.where(valid, "estimated", "undefined_near_zero_difference_variance"),
        "significant_BH": np.isfinite(qvalue) & (qvalue < alpha),
        "direction": np.where(effect > 0, "DLPFC_higher",
                              np.where(effect < 0, "Caudate_higher", "equal")),
        "absolute_mean_paired_difference": np.abs(effect),
    })
    # A stable tie break by gene ID makes the descriptive ranking deterministic.
    result = result.sort_values(["absolute_mean_paired_difference", "gene_id"],
                                ascending=[False, True], kind="stable").reset_index(drop=True)
    result.insert(1, "descriptive_rank", np.arange(1, len(result) + 1))
    long_means = donor_means.T.rename_axis(columns="gene_id").stack().rename(
        "mean_expression").reset_index()
    counts = metadata.groupby(["donor_id", "roi", "hemisphere", "input_level"], sort=True).agg(
        n_tissue_samples=("n_tissue_samples", "sum"), n_input_rows=("sample_id", "size")
    ).reset_index()
    counts["donor_roi_total_samples"] = counts.groupby(["donor_id", "roi"])[
        "n_tissue_samples"].transform("sum")
    differences = pd.DataFrame(difference, index=expression.index, columns=donors)
    differences.index.name = "gene_id"
    return result, long_means, counts, differences


def save_heatmap(result: pd.DataFrame, output: Path, top_n: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = result.head(top_n)
    values = selected[["DLPFC_mean_across_donors", "Caudate_mean_across_donors"]].to_numpy()
    fig, ax = plt.subplots(figsize=(7, max(4, 0.22 * len(selected))))
    plot = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks([0, 1], ["DLPFC", "Caudate"])
    ax.set_yticks(np.arange(len(selected)), selected["gene_id"], fontsize=7)
    ax.set_title(f"Top {len(selected)} by absolute mean donor-paired difference\nDescriptive ranking")
    fig.colorbar(plot, ax=ax, label="Mean expression on supplied normalized scale")
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def run_enrichment(result: pd.DataFrame, output_dir: Path, alpha: float,
                   numeric_namespace: str = "ENTREZGENE") -> dict:
    """Query GO separately by direction, using only successfully tested genes.

    Uses the official HTTPS JSON API. A query with no significant genes is
    recorded as skipped. Raw requests and responses preserve API metadata.
    Returned g:Profiler p_value is already adjusted using the requested FDR
    method; it must not be presented as an uncorrected hypergeometric p-value.
    The correction supplied by g:Profiler is not a correction across the two
    directional analyses jointly.
    """
    background = sorted(result.loc[result["test_status"] == "estimated", "gene_id"].tolist())
    write_json(output_dir / "GO_background_genes.json", background)
    summary = {}
    for direction in ("DLPFC_higher", "Caudate_higher"):
        query = sorted(result.loc[result["significant_BH"] &
                                  result["direction"].eq(direction), "gene_id"].tolist())
        request_body = {
            "organism": "hsapiens", "query": query, "sources": GO_SOURCES,
            "user_threshold": alpha, "all_results": True, "ordered": False,
            "combined": False, "measure_underrepresentation": False,
            "domain_scope": "custom", "background": background,
            "significance_threshold_method": "fdr", "no_evidences": False,
            "numeric_ns": numeric_namespace,
        }
        prefix = output_dir / f"GO_{direction}"
        write_json(prefix.with_name(prefix.name + "_request.json"), request_body)
        if not query:
            summary[direction] = {"status": "skipped_empty_significant_gene_set",
                                  "query_genes": 0, "background_genes": len(background)}
            write_json(prefix.with_suffix(".json"), summary[direction])
            pd.DataFrame(columns=["source", "native", "name", "gprofiler_adjusted_p_value"]
                         ).to_csv(prefix.with_suffix(".csv"), index=False)
            continue
        request = Request(GPROFILER_URL,
                          data=json.dumps(request_body).encode("utf-8"),
                          headers={"Content-Type": "application/json",
                                   "User-Agent": f"BrainBasePairReanalysis/{IMPLEMENTATION_VERSION}"},
                          method="POST")
        try:
            with urlopen(request, timeout=45) as response:
                raw = response.read()
            prefix.with_suffix(".json").write_bytes(raw)
            parsed = json.loads(raw)
            if not isinstance(parsed, dict) or not isinstance(parsed.get("result"), list):
                raise ValueError("g:Profiler returned an unexpected response structure")
            rows = parsed["result"]
            table = pd.json_normalize(rows)
            if rows:
                if "p_value" not in table:
                    raise ValueError("g:Profiler response is missing adjusted p_value")
                try:
                    probabilities = pd.to_numeric(table["p_value"], errors="raise").astype(float)
                except (TypeError, ValueError) as exc:
                    raise ValueError("g:Profiler adjusted p_value must be numeric") from exc
                if (not np.isfinite(probabilities).all()
                        or not probabilities.between(0, 1, inclusive="both").all()
                        or table["p_value"].map(lambda value: isinstance(value, bool)).any()):
                    raise ValueError("g:Profiler adjusted p_value must be finite and within [0, 1]")
                table["p_value"] = probabilities
                table = table.rename(columns={"p_value": "gprofiler_adjusted_p_value"})
                for column in table.columns:
                    table[column] = table[column].map(
                        lambda value: json.dumps(value, ensure_ascii=False)
                        if isinstance(value, (dict, list)) else value)
            else:
                table = pd.DataFrame(columns=["source", "native", "name", "gprofiler_adjusted_p_value"])
            table.to_csv(prefix.with_suffix(".csv"), index=False)
            summary[direction] = {"status": "complete", "query_genes": len(query),
                                  "background_genes": len(background), "returned_terms": len(rows)}
        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            summary[direction] = {"status": "failed", "error": str(exc),
                                  "query_genes": len(query), "background_genes": len(background)}
            write_json(prefix.with_name(prefix.name + "_error.json"), summary[direction])
    write_json(output_dir / "GO_status.json", summary)
    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--expression", type=Path, required=True,
                        help="Real, already normalized gene-by-sample CSV; first column gene_id")
    parser.add_argument("--metadata", type=Path, required=True,
                        help="CSV with sample_id, donor_id, roi, hemisphere")
    parser.add_argument("--output-dir", "--output", dest="output_dir", type=Path, required=True,
                        help="New or empty directory; existing outputs are never overwritten")
    parser.add_argument("--normalization-description", required=True,
                        help="Document the actual upstream expression scale and preprocessing")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Strict q < alpha threshold for gene and GO results; default 0.05")
    parser.add_argument("--top-n", type=int, default=50,
                        help="Descriptive gene ranking length, independent of significance")
    parser.add_argument("--heatmap", action="store_true", help="Write descriptive heatmap")
    parser.add_argument("--enrich", action="store_true",
                        help="Send significant directional gene IDs and tested background to g:Profiler")
    parser.add_argument("--numeric-namespace", default="ENTREZGENE",
                        help="g:Profiler namespace for numeric gene IDs; recorded in request")
    args = parser.parse_args(argv)
    if not 0 < args.alpha < 1:
        parser.error("--alpha must be between 0 and 1")
    if args.top_n < 1:
        parser.error("--top-n must be positive")
    if not args.normalization_description.strip():
        parser.error("--normalization-description must not be blank")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    expression, metadata = load_inputs(args.expression, args.metadata)
    if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
        raise ValueError("Output directory must be new or empty; refusing to overwrite existing files")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = args.output_dir
    started = datetime.now(timezone.utc).isoformat()
    result, means, counts, differences = analyze(expression, metadata, args.alpha)
    result.to_csv(output_dir / "all_gene_statistics.csv", index=False, na_rep="NA")
    result.head(args.top_n).to_csv(output_dir / f"top{args.top_n}_descriptive_genes.csv",
                                 index=False, na_rep="NA")
    for direction in ("DLPFC_higher", "Caudate_higher"):
        result.loc[result["significant_BH"] & result["direction"].eq(direction)].to_csv(
            output_dir / f"significant_{direction}.csv", index=False, na_rep="NA")
    means.to_csv(output_dir / "donor_roi_gene_means.csv", index=False)
    counts.to_csv(output_dir / "donor_roi_sample_counts.csv", index=False)
    differences.to_csv(output_dir / "donor_paired_differences.csv")
    metadata.to_csv(output_dir / "sample_metadata_used.csv", index=False)
    status = "complete"
    optional_status = {}
    if args.heatmap:
        try:
            save_heatmap(result, output_dir / f"top{args.top_n}_descriptive_heatmap.png", args.top_n)
            optional_status["heatmap"] = "complete"
        except Exception as exc:
            optional_status["heatmap"] = {"status": "failed", "error": str(exc)}
            status = "partial_optional_output_failure"
    if args.enrich:
        optional_status["GO_enrichment"] = run_enrichment(result, output_dir, args.alpha,
                                                         args.numeric_namespace)
        if any(value["status"] == "failed" for value in optional_status["GO_enrichment"].values()):
            status = "partial_optional_output_failure"
    versions = {}
    for package in ["numpy", "pandas", "scipy", "matplotlib"]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    provenance = {
        "implementation_version": IMPLEMENTATION_VERSION,
        "status": status,
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": "New analysis implementation; does not reconstruct or validate original manuscript results",
        "python": platform.python_version(), "platform": platform.platform(),
        "packages": versions,
        "arguments": {key: str(value) if isinstance(value, Path) else value
                      for key, value in vars(args).items()},
        "inputs": {
            "expression": {"filename": args.expression.name, "sha256": file_sha256(args.expression)},
            "metadata": {"filename": args.metadata.name, "sha256": file_sha256(args.metadata)},
            "analysis_script": {"filename": Path(__file__).name, "sha256": file_sha256(Path(__file__))},
        },
        "n_genes_supplied": len(result),
        "n_genes_with_estimable_test": int(result["test_status"].eq("estimated").sum()),
        "n_paired_donors": int(metadata["donor_id"].nunique()),
        "donor_ids": sorted(metadata["donor_id"].unique().tolist()),
        "effect_definition": "Mean within-donor DLPFC minus Caudate on supplied normalized scale; not log2 fold change",
        "input_level": str(metadata["input_level"].iloc[0]),
        "aggregation": ("Arithmetic mean over tissue samples within donor and ROI, then equal donor weighting"
                        if metadata["input_level"].iloc[0] == "tissue_sample" else
                        "One supplied summary per donor and ROI, then equal donor weighting; upstream aggregation must be documented"),
        "gene_testing": "Two-sided paired t-test; BH across estimable genes; strict q < alpha",
        "undefined_tests": "Near-zero donor difference variance gives NA p and q and excludes gene from GO background",
        "confidence_intervals": "95 percent pointwise intervals, unadjusted for multiple comparisons",
        "GO_scope": "Separate directional lists from BH-significant genes; custom background of estimable genes; FDR returned by g:Profiler per query, not joint correction across directions",
        "GO_endpoint": GPROFILER_URL if args.enrich else None,
        "optional_outputs": optional_status,
        "limitations": [
            "Input donor identities, gene identity, upstream normalization, and ROI assignment are not independently verified by this script",
            "Hemisphere labels are reported; hemisphere balance is not automatically corrected or modeled",
            "Minimum three donors is an execution guard, not a sample-size justification",
            "Regional AHBA comparisons cannot establish interpreter-specific effects or training-induced change",
            "No spatial null model or independent cohort validation is implemented",
        ],
        "output_sha256": {path.name: file_sha256(path)
                          for path in sorted(output_dir.iterdir()) if path.is_file()},
    }
    write_json(output_dir / "run_manifest.json", provenance)
    print(json.dumps({"status": status, "output_directory": str(output_dir),
                      "genes": len(result), "paired_donors": metadata["donor_id"].nunique(),
                      "significant_genes": int(result["significant_BH"].sum())}))
    return 0 if status == "complete" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
