import unittest

import numpy as np
import pandas as pd

from src.data.extract_era5_weather import (
    longest_true_run,
    nearest_grid,
    summarise_hourly_weather,
)


class Era5WeatherTests(unittest.TestCase):
    def test_nearest_grid_avoids_bankers_rounding(self):
        result = nearest_grid(np.array([-18.125, -18.12, 147.125, 147.12]))
        np.testing.assert_allclose(result, [-18.0, -18.0, 147.25, 147.0])

    def test_longest_true_run(self):
        self.assertEqual(
            longest_true_run([False, True, True, False, True]), 2
        )
        self.assertEqual(longest_true_run([False, False]), 0)

    def test_weather_summary_uses_q1_and_antecedent_december(self):
        time = pd.date_range('2023-12-01', '2024-03-31 23:00', freq='h')
        hourly = pd.DataFrame(
            {
                'time': time,
                'precipitation': 1.0,
                'wind_speed_10m': 2.0,
                'temperature_2m': 28.0,
            }
        )
        result = summarise_hourly_weather(hourly)
        self.assertEqual(result['rain_december_total'], 31 * 24)
        self.assertEqual(result['rain_q1_total'], 91 * 24)
        self.assertEqual(result['q1_n_days'], 91)
        self.assertEqual(result['wind_q1_fraction_below_3'], 1)
        self.assertEqual(result['wind_q1_longest_spell_below_3'], 91)
        self.assertEqual(result['temperature_q1_mean'], 28)


if __name__ == '__main__':
    unittest.main()
