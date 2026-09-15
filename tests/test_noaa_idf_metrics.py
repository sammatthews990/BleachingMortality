import unittest

import numpy as np
import pandas as pd

from src.data.extract_noaa_idf_metrics import (
    empirical_return_period,
    summarise_return_curve,
    summer_block_maxima,
)


class NoaaIdfMetricTests(unittest.TestCase):
    def test_empirical_return_period_is_smoothed_and_upper_tailed(self):
        reference = [1, 2, 3, 4]
        self.assertEqual(empirical_return_period(reference, 5), 5.0)
        self.assertEqual(empirical_return_period(reference, 3), 5 / 3)

    def test_summer_block_maxima_capture_short_and_persistent_heat(self):
        dates = pd.date_range('2019-11-01', '2020-04-30', freq='D')
        values = np.zeros(len(dates))
        values[:3] = 6
        maxima = summer_block_maxima(dates, values, 2020, durations=(1, 3, 7))
        self.assertEqual(maxima[1], 6)
        self.assertEqual(maxima[3], 6)
        self.assertAlmostEqual(maxima[7], 18 / 7)

    def test_return_curve_retains_breadth_of_rarity(self):
        from src.data.extract_noaa_idf_metrics import DURATIONS
        reference = {d: np.array([1, 2, 3, 4]) for d in DURATIONS}
        event = {d: (5 if d <= 7 else 3) for d in DURATIONS}
        result = summarise_return_curve(reference, event)
        self.assertGreater(result['noaa_idf_mean_log_return_period'], 0)
        self.assertLess(
            result['noaa_idf_persistence_logrp_difference'], 0
        )


if __name__ == '__main__':
    unittest.main()
