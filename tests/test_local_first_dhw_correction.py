import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_local_first_dhw_correction import (  # noqa: E402
    aggregate_site_events,
    direct_site_predict,
    enso_historical_effects,
    local_predict,
)


class LocalFirstDhwCorrectionTests(unittest.TestCase):
    def test_local_predict_has_a_hard_spatial_boundary(self):
        sources = pd.DataFrame({
            "lat": [-14.0, -18.0],
            "lon": [145.0, 149.0],
            "value": [4.0, -8.0],
        })
        targets = pd.DataFrame({"lat": [-14.02, -20.0], "lon": [145.02, 152.0]})
        result = local_predict(sources, targets, "value", radius_km=25.0)

        self.assertAlmostEqual(result.loc[0, "prediction"], 4.0)
        self.assertEqual(result.loc[0, "sources_used"], 1)
        self.assertTrue(np.isnan(result.loc[1, "prediction"]))
        self.assertEqual(result.loc[1, "sources_used"], 0)

    def test_local_predict_uses_only_the_nearest_requested_sources(self):
        sources = pd.DataFrame({
            "lat": [-14.00, -14.01, -14.02],
            "lon": [145.00, 145.01, 145.02],
            "value": [1.0, 2.0, 100.0],
        })
        targets = pd.DataFrame({"lat": [-14.0], "lon": [145.0]})
        result = local_predict(
            sources, targets, "value", radius_km=25.0, max_sources=2
        )

        self.assertEqual(result.loc[0, "sources_used"], 2)
        self.assertLess(result.loc[0, "prediction"], 2.0)

    def test_historical_fallback_requires_repeated_phase_and_sign(self):
        history = pd.DataFrame({
            "site": ["stable", "stable", "single", "mixed", "mixed"],
            "event_year": [2016, 2024, 2024, 2016, 2024],
            "enso_phase": ["El Nino"] * 5,
            "lat": [-14.0] * 5,
            "lon": [145.0] * 5,
            "discrepancy": [3.0, 5.0, 6.0, 3.0, -3.0],
        })
        result = enso_historical_effects(history, "El Nino")

        self.assertEqual(result["site"].tolist(), ["stable"])
        self.assertEqual(result.iloc[0]["events"], 2)
        self.assertEqual(result.iloc[0]["effect"], 4.0)

    def test_reef_logger_precedes_automated_station(self):
        rows = pd.DataFrame({
            "event_year": [2024, 2024],
            "site": ["Lizard Island", "Lizard Island"],
            "source": ["temperature_logger", "automated_weather"],
            "depth_class": ["shallow_0_5m", "shallow_0_5m"],
            "habitat_position": ["flat", "unknown"],
            "lat": [-14.67, -14.67],
            "lon": [145.46, 145.46],
            "depth_m": [2.0, 5.0],
            "dhw_discrepancy_logger_minus_noaa": [4.26, -2.37],
        })
        result = aggregate_site_events(rows)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["measurement_source"], "reef_logger")
        self.assertAlmostEqual(result.iloc[0]["discrepancy"], 4.26)


if __name__ == "__main__":
    unittest.main()
