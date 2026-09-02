import numpy as np
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'src' / 'data'))
from extract_sst_chla_features import (  # noqa: E402
    chlorophyll_summary,
    distribution_summary,
    nearest_valid_pixel,
)


class SstChlaFeatureTests(unittest.TestCase):
    def test_distribution_summary_matches_scipy_conventions(self):
        values = np.column_stack([
            np.arange(1.0, 11.0),
            np.array([1, 1, 1, 1, 1, 1, 1, 2, 4, 9], dtype=float),
        ])
        summary = distribution_summary(values, 'sst')
        self.assertTrue(np.all(summary['sst_n'] == 10))
        self.assertTrue(np.isclose(summary['sst_median'][0], 5.5))
        self.assertLess(abs(summary['sst_skewness'][0]), 1e-12)
        self.assertGreater(summary['sst_skewness'][1], 2)
        self.assertGreater(summary['sst_excess_kurtosis'][1], 4)

    def test_shape_metrics_require_eight_observations(self):
        values = np.array(
            [[1], [2], [3], [4], [5], [6], [7], [np.nan]], dtype=float
        )
        summary = distribution_summary(values, 'sst')
        self.assertEqual(summary['sst_n'][0], 7)
        self.assertTrue(np.isnan(summary['sst_skewness'][0]))
        self.assertTrue(np.isnan(summary['sst_excess_kurtosis'][0]))

    def test_chlorophyll_summary_ignores_nonpositive_values(self):
        values = np.array([[0, 1], [-1, 2], [1, 3], [3, np.nan]], dtype=float)
        summary = chlorophyll_summary(values, 'chla')
        self.assertTrue(np.array_equal(summary['chla_n'], [2, 3]))
        self.assertTrue(np.allclose(summary['chla_median'], [2, 2]))
        self.assertTrue(np.allclose(summary['chla_q10'], [1.2, 1.2]))
        self.assertTrue(np.allclose(summary['chla_q90'], [2.8, 2.8]))

    def test_nearest_valid_pixel_recovers_masked_centroid(self):
        lats = np.array([-20.00, -20.01, -20.02])
        lons = np.array([150.00, 150.01, 150.02])
        valid = np.zeros((3, 3), dtype=bool)
        valid[1, 2] = True
        pixel = nearest_valid_pixel(
            valid, lats, lons, -20.01, 150.01, radius_km=5
        )
        self.assertIsNotNone(pixel)
        self.assertEqual(pixel[:2], (1, 2))
        self.assertGreater(pixel[2], 0)
        self.assertLess(pixel[2], 2)

    def test_nearest_valid_pixel_respects_radius(self):
        lats = np.array([-20.0, -20.1])
        lons = np.array([150.0, 150.1])
        valid = np.zeros((2, 2), dtype=bool)
        valid[1, 1] = True
        pixel = nearest_valid_pixel(
            valid, lats, lons, -20.0, 150.0, radius_km=5
        )
        self.assertIsNone(pixel)


if __name__ == '__main__':
    unittest.main()
