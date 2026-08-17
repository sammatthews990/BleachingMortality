import datetime as dt
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.fetch_patmosx_cloud import (
    PatmosFile,
    aggregate_year,
    nanmean_stack,
    sample_reefs,
)


class PatmosAggregationTests(unittest.TestCase):
    def test_nanmean_stack_preserves_all_missing_cells(self):
        first = np.array([[0.2, np.nan], [0.4, 0.6]], dtype=np.float32)
        second = np.array([[0.6, np.nan], [np.nan, 0.8]], dtype=np.float32)
        result = nanmean_stack([first, second])
        np.testing.assert_allclose(result[[0, 1, 1], [0, 0, 1]], [0.4, 0.4, 0.7])
        self.assertTrue(np.isnan(result[0, 1]))

    def test_aggregate_year_keeps_orbit_and_day_weighting_explicit(self):
        day1 = dt.date(2020, 1, 1)
        day2 = dt.date(2020, 1, 2)
        records = [
            PatmosFile("a", "d1_asc.npz", day1, "asc", day1, False),
            PatmosFile("b", "d1_des.npz", day1, "des", day1, False),
            PatmosFile("c", "d2_asc.npz", day2, "asc", day2, False),
        ]
        arrays = {
            "d1_asc.npz": np.array([[0.2, np.nan]], dtype=np.float32),
            "d1_des.npz": np.array([[0.4, 0.6]], dtype=np.float32),
            "d2_asc.npz": np.array([[0.6, 0.8]], dtype=np.float32),
        }
        temporary = Path("tests/.patmosx_cloud_unit_test")
        temporary.mkdir(exist_ok=True)
        try:
            paths = {}
            for name, values in arrays.items():
                path = temporary / name
                np.savez(path, cloud_fraction=values, latitude=[-10.0], longitude=[145.0, 145.1])
                paths[name] = path
            grids, _, _, metadata = aggregate_year(records, paths)
        finally:
            for path in temporary.glob("*.npz"):
                path.unlink()
            temporary.rmdir()

        np.testing.assert_allclose(grids["cloud_fraction_asc_jfm"], [[0.4, 0.8]])
        np.testing.assert_allclose(grids["cloud_fraction_passmean_jfm"], [[0.4, 0.7]])
        np.testing.assert_allclose(grids["cloud_fraction_dailymean_jfm"], [[0.45, 0.7]])
        self.assertEqual(metadata["n_days"], 2)
        self.assertEqual(metadata["n_paired_days"], 1)

    def test_sample_reefs_retains_nearest_and_paper_grid_conventions(self):
        reefs = pd.DataFrame({"LABEL_ID": ["reef-1"], "lon": [142.04], "lat": [-10.04]})
        latitude = np.array([-10.1, -10.0])
        longitude = np.array([142.0, 142.1, 142.2])
        base = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32) / 10
        grids = {
            "cloud_fraction_asc_jfm": base,
            "cloud_fraction_des_jfm": base,
            "cloud_fraction_passmean_jfm": base,
            "cloud_fraction_dailymean_jfm": base,
        }
        metadata = {
            "n_days": 80,
            "n_asc_files": 80,
            "n_des_files": 80,
            "n_paired_days": 80,
            "n_preliminary_files": 160,
        }
        result = sample_reefs(reefs, 2024, grids, latitude, longitude, metadata)
        self.assertAlmostEqual(result.loc[0, "cloud_fraction_asc_jfm"], 0.4)
        self.assertAlmostEqual(result.loc[0, "cloud_fraction_asc_jfm_paper_cell"], 0.5)
        self.assertAlmostEqual(result.loc[0, "cloudp_90"], 0.5)
        self.assertEqual(result.loc[0, "cloud_n_preliminary_files"], 160)


if __name__ == "__main__":
    unittest.main()
