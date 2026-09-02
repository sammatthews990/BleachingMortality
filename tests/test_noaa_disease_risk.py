import numpy as np
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'src' / 'data'))
from extract_noaa_disease_risk import summarise_risk


class DiseaseRiskTests(unittest.TestCase):
    def test_negative_masks_are_not_treated_as_risk(self):
        values = np.array([
            [-10.0, 0.0, 0.5],
            [-10.0, 2.0, 1.5],
            [-10.0, 4.0, np.nan],
        ])
        result = summarise_risk(values)
        np.testing.assert_array_equal(
            result['disease_risk_applicable'], [0, 1, 1]
        )
        np.testing.assert_allclose(result['disease_risk_max'], [0, 4, 1.5])
        np.testing.assert_array_equal(
            result['disease_risk_days_ge1'], [0, 2, 1]
        )
        np.testing.assert_allclose(
            result['disease_risk_burden'], [0, 6 / 7, 2 / 7]
        )


if __name__ == '__main__':
    unittest.main()
