import unittest

import numpy as np

from scripts.extract_ereefs_salinity import summarise_salinity


class EReefsSalinityTests(unittest.TestCase):
    def test_freshwater_exposure_accumulates_daily_deficit(self):
        daily = np.array(
            [
                [29.0, 27.0, np.nan],
                [27.0, 25.0, np.nan],
                [25.0, 23.0, np.nan],
            ]
        )

        result = summarise_salinity(daily)

        np.testing.assert_allclose(
            result['freshwater_exposure_28'][:2], [4.0, 9.0]
        )
        np.testing.assert_allclose(
            result['freshwater_exposure_26'][:2], [1.0, 4.0]
        )
        np.testing.assert_allclose(result['days_below_28'][:2], [2.0, 3.0])
        self.assertTrue(np.isnan(result['freshwater_exposure_28'][2]))
        self.assertEqual(result['salinity_n_days'][2], 0)


if __name__ == '__main__':
    unittest.main()
