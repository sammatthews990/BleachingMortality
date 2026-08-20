import unittest
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from scripts.fetch_dms_environmental_data import compute_ten_year_dhw_metrics


class ThermalHistoryTests(unittest.TestCase):
    def test_window_excludes_current_year_and_keeps_novelty_signed(self):
        years = list(range(2010, 2022))
        history = np.array(
            [[float(year - 2010), float(2020 - year)] for year in years]
        )

        result = compute_ten_year_dhw_metrics(2021, years, history)
        expected_past = history[1:11, :]

        np.testing.assert_allclose(
            result['dhw10_load4'],
            np.maximum(expected_past - 4.0, 0).sum(axis=0),
        )
        np.testing.assert_allclose(
            result['dhw_novelty10'],
            history[11, :] - expected_past.max(axis=0),
        )
        self.assertGreater(result['dhw_novelty10'][0], 0)
        self.assertLess(result['dhw_novelty10'][1], 0)


if __name__ == '__main__':
    unittest.main()
