import unittest

import pandas as pd

from src.data.extract_enso_indices import (
    agreement_label,
    enso_category,
    soi_phase,
    summarise_events,
)


class EnsoIndexTests(unittest.TestCase):
    def test_enso_categories(self):
        self.assertEqual(enso_category(0.49), ("Neutral", "Neutral", "Neutral"))
        self.assertEqual(
            enso_category(0.5), ("El Nino", "Weak", "Weak El Nino")
        )
        self.assertEqual(
            enso_category(-1.2), ("La Nina", "Moderate", "Moderate La Nina")
        )
        self.assertEqual(
            enso_category(1.5), ("El Nino", "Strong", "Strong El Nino")
        )

    def test_soi_direction(self):
        self.assertEqual(soi_phase(-7.1), "El Nino-like")
        self.assertEqual(soi_phase(7.1), "La Nina-like")
        self.assertEqual(soi_phase(0), "Neutral/mixed")

    def test_ocean_atmosphere_agreement(self):
        self.assertEqual(agreement_label("El Nino", "El Nino-like"), "coupled")
        self.assertEqual(agreement_label("La Nina", "Neutral/mixed"), "mixed")

    def test_event_window_uses_djf_jfm_fma_and_previous_december(self):
        roni = pd.DataFrame(
            {
                "seas": ["DJF", "JFM", "FMA"],
                "year": [2022, 2022, 2022],
                "roni": [-1.0, -1.2, -1.4],
            }
        )
        soi = pd.DataFrame(
            {
                "date": pd.to_datetime(
                    ["2021-12-01", "2022-01-01", "2022-02-01", "2022-03-01"]
                ),
                "soi": [8.0, 10.0, 12.0, 14.0],
            }
        )
        result = summarise_events(roni, soi, (2022,)).iloc[0]
        self.assertAlmostEqual(result["roni_bleaching_summer_mean"], -1.2)
        self.assertAlmostEqual(result["soi_dec_mar_mean"], 11.0)
        self.assertEqual(result["enso_category"], "Moderate La Nina")


if __name__ == "__main__":
    unittest.main()
