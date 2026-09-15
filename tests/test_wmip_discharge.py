import unittest

import numpy as np
import pandas as pd

from src.data.fetch_wmip_discharge import (
    parse_wmip_traces,
    quality_status,
    summarise_catchment_events,
)


class WmipDischargeTests(unittest.TestCase):
    def test_unavailable_zero_becomes_missing(self):
        response = {
            'error_num': 0,
            'licence': 'CC_40_BY',
            'return': {'traces': [{
                'error_num': 0,
                'site': 'A',
                'site_details': {
                    'name': 'Test River', 'longitude': '145.0',
                    'latitude': '-16.0', 'timezone': '10.0'
                },
                'quality_codes': {'151': 'Data not yet available'},
                'varto_details': {'units': 'Cumecs'},
                'trace': [{'v': '0.0000', 't': 20240102000000, 'q': 151}],
            }]},
        }
        result = parse_wmip_traces(
            response, 'snapshot', '2024-05-01T00:00:00+00:00', 2024
        )
        self.assertEqual(result.loc[0, 'date_local'], '2024-01-01')
        self.assertTrue(np.isnan(result.loc[0, 'discharge_mean_m3_s']))
        self.assertEqual(result.loc[0, 'value_status'], 'unavailable')

    def test_quality_mapping_is_conservative(self):
        self.assertEqual(quality_status('Good'), 'validated')
        self.assertEqual(quality_status('Provisional'), 'provisional')
        self.assertEqual(quality_status('Fair'), 'flagged')

    def test_event_summary_uses_declared_window(self):
        dates = pd.date_range('2023-12-01', '2024-04-30', freq='D')
        gauge = pd.DataFrame({
            'snapshot_id': 's', 'event_year': 2024, 'station_id': 'A',
            'date_local': dates.strftime('%Y-%m-%d'),
            'retrieved_at_utc': '2026-09-07T00:00:00+00:00',
            'discharge_mean_m3_s': 2.0,
            'value_status': 'validated',
        })
        crosswalk = pd.DataFrame([{
            'station_id': 'A', 'station_weight': 1.0, 'source_id': 'river',
            'basin_id': 'B', 'catchment_name': 'River', 'routing_lon': 145.0,
            'routing_lat': -16.0, 'kernel_scale_km': 45,
            'max_distance_km': 150, 'pilot_role': 'test',
            'selection_status': 'pilot_only', 'regulation_flag': 'unknown',
            'rating_audit_status': 'requires_review',
        }])
        result = summarise_catchment_events(gauge, crosswalk)
        march = result.loc[result['cutoff_name'].eq('march')].iloc[0]
        april = result.loc[result['cutoff_name'].eq('april')].iloc[0]
        self.assertEqual(march['expected_days'], 122)
        self.assertEqual(april['expected_days'], 152)
        self.assertAlmostEqual(march['discharge_total_ml'], 122 * 2 * 86.4)
        self.assertAlmostEqual(april['coverage_fraction'], 1.0)
        self.assertEqual(april['product_mode'], 'environmental_hindcast')


if __name__ == '__main__':
    unittest.main()
