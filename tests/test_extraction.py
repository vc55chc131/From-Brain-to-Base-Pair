"""Optional extractor checks using synthetic fixtures, never scientific results.

Run in the requirements-ahba.txt environment:
python -m unittest discover -s tests -p test_extraction.py
"""
import importlib.util
from pathlib import Path
import tempfile
import unittest

MODULE_PATH = Path(__file__).resolve().parents[1] / "source_code" / "extract_ahba.py"
spec = importlib.util.spec_from_file_location("extract_ahba", MODULE_PATH)
extract = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extract)

try:
    import abagen
    import numpy as np
    import pandas as pd
    AHBA_AVAILABLE = True
except ImportError:
    AHBA_AVAILABLE = False


class SeedTests(unittest.TestCase):
    def test_unverified_seeds_fail_before_data_access(self):
        path = MODULE_PATH.parent.parent / "config" / "roi_seeds_unverified.csv"
        with self.assertRaisesRegex(ValueError, "provenance is unverified"):
            extract.read_seeds(path)


@unittest.skipUnless(AHBA_AVAILABLE, "Optional AHBA dependencies are not installed")
class ExtractionTests(unittest.TestCase):
    def test_abagen_aggregation_api_and_missing_regions(self):
        frame = pd.DataFrame({"GENE_A": [2.0, 4.0, 8.0]},
                             index=pd.Index([1, 1, 2], name="label"))
        result = abagen.samples_.aggregate_samples({"9861": frame}, labels=[1, 2, 3],
                                                  region_agg="donors", agg_metric="mean",
                                                  return_donors=True)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["9861"].loc[1, "GENE_A"], 3.0)
        self.assertTrue(np.isnan(result["9861"].loc[3, "GENE_A"]))

    def test_weighted_caudate_preserves_donor_and_count(self):
        frame = pd.DataFrame({"GENE_A": [1.0, 2.0, 8.0], "GENE_B": [2.0, 4.0, 6.0]}, index=[1, 2, 3])
        counts = pd.DataFrame({"9861": [2, 3, 1], "10021": [0, 2, 0]}, index=[1, 2, 3])
        seeds = [{"id": 1, "roi": "DLPFC", "hemisphere": "L"},
                 {"id": 2, "roi": "Caudate", "hemisphere": "L"},
                 {"id": 3, "roi": "Caudate", "hemisphere": "R"}]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            diagnostics = extract.export_expression({"9861": frame, "10021": frame}, counts, seeds, out, np, pd)
            matrix = pd.read_csv(out / "expression.csv", index_col=0)
            meta = pd.read_csv(out / "sample_metadata.csv", index_col="sample_id")
            self.assertEqual(matrix.loc["GENE_A", "donor_9861_Caudate"], 3.5)
            self.assertEqual(meta.loc["donor_9861_Caudate", "n_tissue_samples"], 4)
            self.assertEqual(meta.loc["donor_9861_Caudate", "hemisphere"], "B")
            self.assertEqual(meta.loc["donor_10021_Caudate", "hemisphere"], "L")
            self.assertNotIn("donor_10021_DLPFC", matrix.columns)
            self.assertEqual(diagnostics["n_paired_donors"], 1)
            paired_matrix = pd.read_csv(out / "paired_expression.csv", index_col=0)
            paired_meta = pd.read_csv(out / "paired_sample_metadata.csv")
            self.assertEqual(paired_matrix.columns.tolist(), ["donor_9861_DLPFC", "donor_9861_Caudate"])
            self.assertEqual(paired_meta["donor_id"].unique().tolist(), ["H0351.2001"])
            self.assertEqual(diagnostics["paired_analysis_exclusions"], [{"donor_id": "H0351.2002", "observed_rois": ["Caudate"], "reason": "incomplete_ROI_pair"}])

    def test_paired_exports_load_in_companion_with_incomplete_donor(self):
        companion_spec = importlib.util.spec_from_file_location("companion", MODULE_PATH.parent / "DLPFC_vs_Caudate_Zscore.py")
        companion = importlib.util.module_from_spec(companion_spec)
        companion_spec.loader.exec_module(companion)
        frame = pd.DataFrame({"GENE_A": [1.0, 2.0, 8.0], "GENE_B": [2.0, 4.0, 6.0]}, index=[1, 2, 3])
        donors = ["9861", "10021", "12876", "14380"]
        counts = pd.DataFrame({donor: [2, 3, 1] for donor in donors}, index=[1, 2, 3])
        counts.loc[1, "14380"] = 0
        seeds = [{"id": 1, "roi": "DLPFC", "hemisphere": "L"},
                 {"id": 2, "roi": "Caudate", "hemisphere": "L"},
                 {"id": 3, "roi": "Caudate", "hemisphere": "R"}]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            diagnostics = extract.export_expression({donor: frame for donor in donors}, counts, seeds, out, np, pd)
            matrix, metadata = companion.load_inputs(out / "paired_expression.csv", out / "paired_sample_metadata.csv")
            self.assertEqual(diagnostics["n_paired_donors"], 3)
            self.assertEqual(metadata["donor_id"].nunique(), 3)
            self.assertEqual(matrix.shape[1], 6)
            with self.assertRaisesRegex(ValueError, "incomplete donors"):
                companion.load_inputs(out / "expression.csv", out / "sample_metadata.csv")

    def test_no_coverage_never_produces_expression_matrix(self):
        frame = pd.DataFrame({"GENE_A": [np.nan, np.nan, np.nan]}, index=[1, 2, 3])
        counts = pd.DataFrame({"9861": [0, 0, 0]}, index=[1, 2, 3])
        seeds = [{"id": 1, "roi": "DLPFC", "hemisphere": "L"},
                 {"id": 2, "roi": "Caudate", "hemisphere": "L"},
                 {"id": 3, "roi": "Caudate", "hemisphere": "R"}]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            with self.assertRaisesRegex(ValueError, "No tissue samples"):
                extract.export_expression({"9861": frame}, counts, seeds, out, np, pd)
            self.assertFalse((out / "expression.csv").exists())


if __name__ == "__main__":
    unittest.main()
