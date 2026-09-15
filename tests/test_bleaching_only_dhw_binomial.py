import unittest
from pathlib import Path


class BleachingOnlyDhwBinomialContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.script = (
            cls.root / 'src' / 'evaluation' /
            'test_bleaching_only_dhw_binomial.R'
        ).read_text(encoding='utf-8')

    def test_response_is_occurrence_without_invented_trials(self):
        self.assertIn('has_mortality = mortality_prop > 0', self.script)
        self.assertIn('has_mortality ~ dhw_input', self.script)
        self.assertNotIn('cbind(', self.script)

    def test_explicit_nonthermal_rows_are_excluded(self):
        for flag in (
            'disturbance_has_cots',
            'disturbance_has_cyclone',
            'disturbance_has_flood',
        ):
            self.assertIn(flag, self.script)
        self.assertIn('filter(!explicit_nonthermal)', self.script)

    def test_transfer_schemes_and_six_events_are_explicit(self):
        for scheme in (
            'leave_one_event_out',
            'leave_one_sector_out',
            'reef_blocked_5fold',
            'leave_one_programme_out',
        ):
            self.assertIn(scheme, self.script)
        self.assertIn('2025L', self.script)

    def test_2025_missing_adjustment_is_labelled(self):
        self.assertIn('unavailable_2025_no_adjustment', self.script)


if __name__ == '__main__':
    unittest.main()
