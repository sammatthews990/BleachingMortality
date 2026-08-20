import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from secchi_extraction import extract_q1_secchi  # noqa: E402


class SecchiExtractionTests(unittest.TestCase):
    def test_acute_optical_metric_uses_q90_k490(self):
        times = pd.to_datetime(['2024-01-01', '2024-01-02'])
        values = np.array([0.17, 0.51]).reshape(2, 1, 1)
        dataset = xr.Dataset(
            {'K_490': (('time', 'latitude', 'longitude'), values)},
            coords={'time': times, 'latitude': [0.0], 'longitude': [0.0]},
        )
        reefs = pd.DataFrame({'lon': [0.0], 'lat': [0.0]})

        result = extract_q1_secchi(dataset, reefs, 2024)

        expected_q90 = np.quantile([0.17, 0.51], 0.90)
        self.assertAlmostEqual(result.loc[0, 'k490_q90'], expected_q90)
        self.assertAlmostEqual(result.loc[0, 'secc3m_p10'], 1.7 / expected_q90)

    def test_masked_centroid_uses_nearby_water_without_changing_valid_point(self):
        times = pd.to_datetime(["2024-01-01", "2024-01-02"])
        latitudes = np.array([0.01, 0.00, -0.01])
        longitudes = np.array([0.00, 0.01, 0.02])
        values = np.full((2, 3, 3), np.nan)

        # The first reef has a valid centroid cell: 1.7 / 0.34 = 5 m.
        values[:, 0, 0] = 0.34
        # The second centroid is masked. Its 2 km neighbourhood contains four
        # observations at Kd490 0.17 and two at 0.34.
        values[:, 1, 0] = 0.17
        values[:, 1, 2] = 0.17

        dataset = xr.Dataset(
            {"K_490": (("time", "latitude", "longitude"), values)},
            coords={
                "time": times,
                "latitude": latitudes,
                "longitude": longitudes,
            },
        )
        reefs = pd.DataFrame(
            {
                "lon": [0.00, 0.01],
                "lat": [0.01, 0.00],
            }
        )

        result = extract_q1_secchi(dataset, reefs, 2024)

        self.assertEqual(result.loc[0, "secc_match_method"], "centroid_grid_cell")
        self.assertAlmostEqual(result.loc[0, "secc3m"], 5.0)
        self.assertEqual(result.loc[0, "secc_radius_km"], 0.0)

        self.assertEqual(result.loc[1, "secc_match_method"], "water_neighbourhood")
        self.assertTrue(np.isclose(result.loc[1, "secc3m"], 7.5))
        self.assertEqual(result.loc[1, "secc_radius_km"], 2.0)
        self.assertEqual(result.loc[1, "secc_n_observations"], 6)


if __name__ == "__main__":
    unittest.main()
