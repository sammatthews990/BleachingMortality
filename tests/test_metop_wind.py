import unittest

import pandas as pd

from scripts.extract_metop_wind import summarise_daily_wind, tile_details


class MetopWindTests(unittest.TestCase):
    def test_tile_names_span_northern_and_southern_gbr(self):
        self.assertEqual(
            tile_details(-17.5, 146.3),
            (
                '020S_140E',
                'IMOS_SRS-Surface-Waves_M_Wind-METOP-B_FV02_018S-146E-DM00.nc',
            ),
        )
        self.assertEqual(
            tile_details(-20.2, 149.1),
            (
                '040S_140E',
                'IMOS_SRS-Surface-Waves_M_Wind-METOP-B_FV02_021S-149E-DM00.nc',
            ),
        )

    def test_daily_summary_does_not_bridge_missing_days(self):
        dates = pd.date_range('2024-01-01', '2024-01-06', freq='D')
        daily = pd.Series(
            [2.0, 2.5, 2.0, 4.0, 6.0],
            index=dates[[0, 1, 3, 4, 5]],
        )
        summary = summarise_daily_wind(daily, dates)
        self.assertEqual(summary['wind_n_days'], 5)
        self.assertEqual(summary['wind_days_below_3'], 3)
        self.assertEqual(summary['wind_longest_spell_below_3'], 2)
        self.assertAlmostEqual(summary['wind_fraction_below_3'], 3 / 5)


if __name__ == '__main__':
    unittest.main()
