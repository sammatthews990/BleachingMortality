import unittest
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / 'src' / 'data'))
from src.data.fetch_dms_environmental_data import (
    compute_repeated_exposure_metrics,
    compute_ten_year_dhw_metrics,
)


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

    def test_recent_history_is_prior_only_and_uses_strict_threshold(self):
        years = list(range(2015, 2025))
        history = np.zeros((len(years), 3), dtype=float)
        history[years.index(2016), :] = [7.0, 6.0, 2.0]
        history[years.index(2020), :] = [8.0, 7.0, 3.0]
        history[years.index(2022), :] = [9.0, 5.0, 4.0]
        history[years.index(2024), :] = [20.0, 20.0, 20.0]

        result = compute_repeated_exposure_metrics(2024, years, history)

        np.testing.assert_array_equal(
            result['dhw_events_since2016_n6'], [3, 1, 0]
        )
        np.testing.assert_array_equal(
            result['dhw_years_since_last_n6'][:2], [2, 4]
        )
        np.testing.assert_array_equal(
            result['dhw_years_since_last_n6_capped8'], [2, 4, 8]
        )
        self.assertTrue(np.isnan(result['dhw_years_since_last_n6'][2]))
        np.testing.assert_array_equal(
            result['dhw_no_prior_n6'], [0, 0, 1]
        )
        np.testing.assert_array_equal(
            result['dhw_history_years_since2016'], [8, 8, 8]
        )

    def test_recent_history_for_2016_has_zero_length_prior_window(self):
        years = list(range(2015, 2018))
        history = np.ones((len(years), 2), dtype=float) * 10

        result = compute_repeated_exposure_metrics(2016, years, history)

        np.testing.assert_array_equal(
            result['dhw_events_since2016_n6'], [0, 0]
        )
        np.testing.assert_array_equal(result['dhw_no_prior_n6'], [0, 0])
        np.testing.assert_array_equal(
            result['dhw_years_since_last_n6'], [1, 1]
        )


if __name__ == '__main__':
    unittest.main()
