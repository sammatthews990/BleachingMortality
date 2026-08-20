import unittest

import numpy as np

from scripts.extract_ereefs_wind import summarise_wind_matrix


class EreefsWindTests(unittest.TestCase):
    def test_calm_days_and_spells_are_columnwise(self):
        wind = np.array([
            [2.0, 6.0],
            [2.5, 4.0],
            [4.0, np.nan],
            [2.0, 4.5],
        ])
        summary = summarise_wind_matrix(wind)
        np.testing.assert_array_equal(summary['ereefs_wind_n_days'], [4, 3])
        np.testing.assert_array_equal(
            summary['ereefs_wind_days_below_3'], [3, 0]
        )
        np.testing.assert_array_equal(
            summary['ereefs_wind_longest_spell_below_3'], [2, 0]
        )
        np.testing.assert_allclose(
            summary['ereefs_wind_fraction_below_5'], [1, 2 / 3]
        )


if __name__ == '__main__':
    unittest.main()
