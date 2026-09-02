import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src' / 'data'))

from build_local_first_dhw_correction import (  # noqa: E402
    direct_site_predict,
    estimate_uncertainty_calibration,
)


class LocalDhwUncertaintyTests(unittest.TestCase):
    def test_uncertainty_scales_are_heldout_and_source_specific(self):
        site_events = pd.DataFrame({
            'discrepancy_sd': [0.4, 0.8, np.nan],
        })
        validation = pd.DataFrame({
            'local_radius_km': [25.0, 25.0, 25.0, 25.0],
            'correction_source': [
                'current_local', 'current_local', 'none', 'none'
            ],
            'prediction': [1.0, -1.0, 0.0, 0.0],
            'discrepancy': [0.0, 0.0, 3.0, -3.0],
            'local_discrepancy_sd': [0.6, 0.6, np.nan, np.nan],
        })
        result = estimate_uncertainty_calibration(site_events, validation)

        self.assertAlmostEqual(result['direct_measurement_sd_floor'], 0.6)
        self.assertAlmostEqual(result['local_heldout_rmse'], 1.0)
        self.assertAlmostEqual(
            result['local_interpolation_sd_floor'], 0.8
        )
        self.assertAlmostEqual(result['unsupported_correction_sd'], 3.0)

    def test_direct_name_match_cannot_jump_to_distant_duplicate(self):
        sources = pd.DataFrame({
            'site': ['Snake Reef'],
            'lat': [-10.0],
            'lon': [145.0],
            'discrepancy': [5.0],
            'discrepancy_sd': [0.5],
        })
        targets = pd.DataFrame({
            'ReefName': ['Snake Reef (14-087)'],
            'lat': [-20.0],
            'lon': [152.0],
        })
        result = direct_site_predict(sources, targets, 'discrepancy')

        self.assertTrue(np.isnan(result.iloc[0]['prediction']))
        self.assertTrue(np.isnan(result.iloc[0]['measurement_sd']))


if __name__ == '__main__':
    unittest.main()
