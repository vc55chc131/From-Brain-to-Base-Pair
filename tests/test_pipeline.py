"""Unit tests use artificial numbers only, never manuscript evidence or outputs."""
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "source_code" / "DLPFC_vs_Caudate_Zscore.py"
SPEC = importlib.util.spec_from_file_location("pipeline", SCRIPT)
pipeline = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pipeline)


class PipelineTests(unittest.TestCase):
    def fixture(self, folder):
        # Donor A has two cortical samples. Donor averaging must precede inference.
        expression = pd.DataFrame({
            "gene_id": ["TEST_ONLY_A", "TEST_ONLY_B"],
            "a1": [0, 5], "a2": [6, 7], "ac": [1, 2],
            "b1": [6, 4], "bc": [2, 2], "c1": [4, 6], "cc": [3, 3],
        })
        metadata = pd.DataFrame({
            "sample_id": ["a1", "a2", "ac", "b1", "bc", "c1", "cc"],
            "donor_id": ["A", "A", "A", "B", "B", "C", "C"],
            "roi": ["DLPFC", "DLPFC", "Caudate", "DLPFC", "Caudate", "DLPFC", "Caudate"],
            "hemisphere": ["L"] * 7,
        })
        expression_path = folder / "test_expression.csv"
        metadata_path = folder / "test_metadata.csv"
        expression.to_csv(expression_path, index=False)
        metadata.to_csv(metadata_path, index=False)
        return expression_path, metadata_path

    def test_bh_known_values_with_missing(self):
        actual = pipeline.bh_adjust(np.array([0.01, 0.04, 0.03, 0.8, np.nan]))
        np.testing.assert_allclose(actual[:4], [0.04, 0.053333333333, 0.053333333333, 0.8])
        self.assertTrue(np.isnan(actual[4]))

    def test_donor_pairing_and_unequal_sample_weight(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = self.fixture(Path(folder))
            expression, metadata = pipeline.load_inputs(*paths)
            result, means, counts, difference = pipeline.analyze(expression, metadata)
            row = result.set_index("gene_id").loc["TEST_ONLY_A"]
            self.assertAlmostEqual(row["mean_paired_difference_DLPFC_minus_Caudate"], 7 / 3)
            np.testing.assert_allclose(difference.loc["TEST_ONLY_A"], [2, 4, 1])
            self.assertEqual(row["degrees_of_freedom"], 2)
            self.assertAlmostEqual(row["t_statistic"], np.sqrt(7))
            self.assertAlmostEqual(row["p_value_two_sided"], 0.11808289631180315)
            self.assertAlmostEqual(row["ci95_low_unadjusted"], -1.4612497002634264)
            self.assertAlmostEqual(row["ci95_high_unadjusted"], 6.127916366930093)
            self.assertEqual(len(means), 12)
            self.assertEqual(int(counts["n_tissue_samples"].sum()), 7)
            # Reordering metadata does not change pairing or the estimates.
            metadata.sample(frac=1, random_state=42).to_csv(paths[1], index=False)
            expression2, metadata2 = pipeline.load_inputs(*paths)
            result2, *_ = pipeline.analyze(expression2, metadata2)
            pd.testing.assert_frame_equal(result, result2)

    def test_missing_pair_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = self.fixture(Path(folder))
            expression = pd.read_csv(paths[0]).drop(columns="cc")
            metadata = pd.read_csv(paths[1]).query("sample_id != 'cc'")
            expression.to_csv(paths[0], index=False)
            metadata.to_csv(paths[1], index=False)
            with self.assertRaisesRegex(ValueError, "both ROIs"):
                pipeline.load_inputs(*paths)

    def test_nonfinite_duplicate_and_unmatched_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = self.fixture(Path(folder))
            frame = pd.read_csv(paths[0])
            frame["a1"] = frame["a1"].astype(float)
            frame.loc[0, "a1"] = float("inf")
            frame.to_csv(paths[0], index=False)
            with self.assertRaisesRegex(ValueError, "finite"):
                pipeline.load_inputs(*paths)
            paths = self.fixture(Path(folder))
            raw = paths[0].read_text().replace("gene_id,a1,a2,", "gene_id,a1,a1,")
            paths[0].write_text(raw)
            with self.assertRaisesRegex(ValueError, "duplicate column"):
                pipeline.load_inputs(*paths)
            paths = self.fixture(Path(folder))
            metadata = pd.read_csv(paths[1])
            metadata.loc[0, "sample_id"] = "unknown"
            metadata.to_csv(paths[1], index=False)
            with self.assertRaisesRegex(ValueError, "match exactly"):
                pipeline.load_inputs(*paths)

    def test_zero_variance_not_falsely_significant(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = self.fixture(Path(folder))
            expression, metadata = pipeline.load_inputs(*paths)
            expression.loc[:, metadata.loc[metadata.roi.eq("DLPFC"), "sample_id"]] = 5
            expression.loc[:, metadata.loc[metadata.roi.eq("Caudate"), "sample_id"]] = 2
            result, *_ = pipeline.analyze(expression, metadata)
            self.assertTrue(result["p_value_two_sided"].isna().all())
            self.assertFalse(result["significant_BH"].any())

    def test_finite_input_overflow_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = self.fixture(Path(folder))
            expression, metadata = pipeline.load_inputs(*paths)
            expression.loc[:, metadata.loc[metadata.roi.eq("DLPFC"), "sample_id"]] = 1e308
            expression.loc[:, metadata.loc[metadata.roi.eq("Caudate"), "sample_id"]] = -1e308
            with self.assertRaisesRegex(ValueError, "[Nn]onfinite"):
                pipeline.analyze(expression, metadata)

    def test_no_query_for_empty_significant_set(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = self.fixture(Path(folder))
            expression, metadata = pipeline.load_inputs(*paths)
            result, *_ = pipeline.analyze(expression, metadata)
            result["significant_BH"] = False
            with patch.object(pipeline, "urlopen") as network:
                summary = pipeline.run_enrichment(result, Path(folder), 0.05)
                network.assert_not_called()
            self.assertTrue(all(item["status"] == "skipped_empty_significant_gene_set"
                                for item in summary.values()))

    def test_summary_level_preserves_effect_and_observed_tissue_counts(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            paths = self.fixture(folder)
            expression, metadata = pipeline.load_inputs(*paths)
            expected, means, *_ = pipeline.analyze(expression, metadata)
            summary = {}
            summary_metadata = []
            for (donor, roi), group in metadata.groupby(["donor_id", "roi"]):
                sample_id = f"{donor}_{roi}"
                summary[sample_id] = expression[group["sample_id"]].mean(axis=1)
                summary_metadata.append({"sample_id": sample_id, "donor_id": donor,
                                         "roi": roi, "hemisphere": "B" if roi == "Caudate" else "L",
                                         "input_level": "donor_roi_summary",
                                         "n_tissue_samples": len(group)})
            pd.DataFrame(summary).to_csv(paths[0])
            pd.DataFrame(summary_metadata).to_csv(paths[1], index=False)
            summary_expression, summary_meta = pipeline.load_inputs(*paths)
            actual, _, counts, _ = pipeline.analyze(summary_expression, summary_meta)
            pd.testing.assert_frame_equal(expected, actual)
            self.assertEqual(int(counts.n_tissue_samples.sum()), 7)
            self.assertEqual(int(counts.n_input_rows.sum()), 6)

    def test_enrichment_records_explicit_background_and_adjusted_p(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            paths = self.fixture(folder)
            expression, metadata = pipeline.load_inputs(*paths)
            result, *_ = pipeline.analyze(expression, metadata)
            result["significant_BH"] = [True, False]
            response = {"result": [{"source": "GO:BP", "native": "GO:TEST_ONLY",
                                      "name": "ARTIFICIAL TEST RESPONSE", "p_value": 0.03}],
                        "meta": {"version": "TEST_ONLY"}}
            with patch.object(pipeline, "urlopen", return_value=io.BytesIO(json.dumps(response).encode())) as network:
                status = pipeline.run_enrichment(result, folder, 0.05)
            self.assertEqual(network.call_count, 1)
            sent = json.loads(network.call_args.args[0].data)
            self.assertEqual(sent["domain_scope"], "custom")
            self.assertEqual(sent["background"], ["TEST_ONLY_A", "TEST_ONLY_B"])
            self.assertEqual(sent["significance_threshold_method"], "fdr")
            self.assertEqual(status["DLPFC_higher"]["status"], "complete")
            output = pd.read_csv(folder / "GO_DLPFC_higher.csv")
            self.assertAlmostEqual(output.loc[0, "gprofiler_adjusted_p_value"], 0.03)

    def test_enrichment_invalid_probability_is_failed(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            paths = self.fixture(folder)
            expression, metadata = pipeline.load_inputs(*paths)
            result, *_ = pipeline.analyze(expression, metadata)
            result["significant_BH"] = [True, False]
            for probability in [float("nan"), float("inf"), -0.1, 1.1, "invalid", True]:
                with self.subTest(probability=probability):
                    response = {"result": [{"p_value": probability}]}
                    with patch.object(pipeline, "urlopen", return_value=io.BytesIO(json.dumps(response).encode())):
                        status = pipeline.run_enrichment(result, folder, 0.05)
                    self.assertEqual(status["DLPFC_higher"]["status"], "failed")

    def test_end_to_end_manifest_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            paths = self.fixture(folder)
            args = ["--expression", str(paths[0]), "--metadata", str(paths[1]),
                    "--output-dir", str(folder / "test_outputs"),
                    "--normalization-description", "ARTIFICIAL UNIT TEST FIXTURE ONLY"]
            self.assertEqual(pipeline.main(args), 0)
            manifest = json.loads((folder / "test_outputs/run_manifest.json").read_text())
            self.assertEqual(manifest["n_paired_donors"], 3)
            self.assertEqual(manifest["n_genes_supplied"], 2)
            self.assertIn("all_gene_statistics.csv", manifest["output_sha256"])
            self.assertEqual(len(manifest["inputs"]["expression"]["sha256"]), 64)
            with self.assertRaisesRegex(ValueError, "overwrite"):
                pipeline.main(args)


if __name__ == "__main__":
    unittest.main()
