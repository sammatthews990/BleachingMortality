import csv
import unittest
from pathlib import Path


class AerialNowcastContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]

    def test_crosswalk_has_both_scales_and_unique_native_categories(self):
        path = cls_path = self.root / 'config' / 'aerial_bleaching_score_crosswalk.csv'
        with cls_path.open(newline='', encoding='utf-8-sig') as handle:
            rows = list(csv.DictReader(handle))
        keys = [(row['source_scale'], row['native_category']) for row in rows]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(
            {row['source_scale'] for row in rows},
            {'historical_hughes', 'contemporary_aims_sop11'},
        )
        for row in rows:
            midpoint = float(row['midpoint_bleached_cover'])
            self.assertGreaterEqual(midpoint, 0.0, path)
            self.assertLessEqual(midpoint, 1.0, path)

    def test_unexpected_historical_score_is_explicitly_top_coded(self):
        path = self.root / 'config' / 'aerial_bleaching_score_crosswalk.csv'
        with path.open(newline='', encoding='utf-8-sig') as handle:
            rows = list(csv.DictReader(handle))
        row = next(
            item for item in rows
            if item['source_scale'] == 'historical_hughes'
            and item['native_category'] == '5'
        )
        self.assertEqual(float(row['midpoint_bleached_cover']), 0.8)
        self.assertIn('topcoded', row['review_action'])

    def test_locked_forecast_code_excludes_update_inputs(self):
        script = (
            self.root / 'src' / 'evaluation' /
            'build_locked_2025_initial_forecast.R'
        ).read_text(encoding='utf-8')
        self.assertIn('maximum_training_event', script)
        self.assertIn('assessment_responses_enter_fit = FALSE', script)
        self.assertIn('aerial_or_rhis_used = FALSE', script)

    def test_promotion_code_requires_locked_2025_guard(self):
        script = (
            self.root / 'src' / 'evaluation' /
            'test_spatial_early_bleaching_update.R'
        ).read_text(encoding='utf-8')
        self.assertIn('event_year == 2025', script)
        self.assertIn('prospective_2025_guard_both_designs', script)

    def test_consensus_update_is_constrained_and_event_centred(self):
        script = (
            self.root / 'src' / 'evaluation' /
            'test_spatial_early_bleaching_update.R'
        ).read_text(encoding='utf-8')
        self.assertIn("consensus_candidate <- 'Reliability-weighted", script)
        self.assertIn('standardise_event_signal', script)
        self.assertIn('consensus_centre', script)
        self.assertIn('fit_monotone_consensus_update', script)
        self.assertIn("optimize(objective, interval = c(0, upper_bound)", script)
        self.assertIn("'event-centred signal; magnitude fixed'", script)
        self.assertIn('diagnostic_direction_guard_both_designs', script)


if __name__ == '__main__':
    unittest.main()
