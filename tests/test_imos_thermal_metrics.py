import unittest

import numpy as np
import pandas as pd

from src.data.extract_imos_thermal_metrics import (
    fill_short_gaps,
    interpolate_monthly_climatology,
    longest_true_run,
    reconstruct_sstaars_daily,
    sstaars_mmm,
    summarise_marine_heatwaves,
    summarise_percentile_heat_dose,
    summarise_day_night_relief,
    summarise_intensity_duration_profile,
    summarise_site_year,
)


def constant_coefficients(sites=1, temperature=28.0):
    coefficients = {'tm': np.full(sites, temperature)}
    for harmonic in range(1, 5):
        stem = 'ta' if harmonic == 1 else f't{harmonic}a'
        coefficients[f'{stem}_real'] = np.zeros(sites)
        coefficients[f'{stem}_imag'] = np.zeros(sites)
    return coefficients


class ImosThermalMetricTests(unittest.TestCase):
    def test_sstaars_reconstruction_uses_documented_conjugate_sign(self):
        coefficients = constant_coefficients()
        coefficients['ta_real'][0] = 2.0
        coefficients['ta_imag'][0] = 3.0
        quarter_cycle = 365.25 / 4
        reconstructed = reconstruct_sstaars_daily(
            coefficients, quarter_cycle
        )
        self.assertAlmostEqual(reconstructed[0], 31.0, places=10)

    def test_mmm_of_constant_climatology_is_constant(self):
        coefficients = constant_coefficients(sites=2, temperature=27.5)
        np.testing.assert_allclose(sstaars_mmm(coefficients), [27.5, 27.5])

    def test_short_gap_fill_does_not_partially_fill_long_gaps(self):
        values = np.array([1.0, np.nan, np.nan, 4.0, np.nan, np.nan, np.nan, 8.0])
        filled = fill_short_gaps(values, max_gap=2)
        np.testing.assert_allclose(filled[:4], [1.0, 2.0, 3.0, 4.0])
        self.assertTrue(np.isnan(filled[4:7]).all())
        self.assertEqual(filled[7], 8.0)

    def test_longest_true_run(self):
        self.assertEqual(
            longest_true_run([False, True, True, False, True, True, True]),
            3,
        )

    def test_noaa_equivalent_dhw_keeps_full_hotspot_above_one(self):
        dates = pd.date_range('2019-08-01', '2020-12-31', freq='D')
        sst = np.full(len(dates), 28.0)
        hot = (dates >= '2020-01-01') & (dates <= '2020-03-24')
        sst[hot] = 30.0
        result = summarise_site_year(dates, sst, 28.0, 2020)
        self.assertAlmostEqual(result['imos_dhw1_84_max'], 24.0, places=8)
        self.assertAlmostEqual(
            result['imos_dhw1_84_max_observed'], 24.0, places=8
        )

    def test_subthreshold_hotspot_is_excluded_from_dhw(self):
        dates = pd.date_range('2019-08-01', '2020-12-31', freq='D')
        sst = np.full(len(dates), 28.9)
        result = summarise_site_year(dates, sst, 28.0, 2020)
        self.assertAlmostEqual(result['imos_dhw1_84_max'], 0.0, places=8)

    def test_event_metrics_capture_duration_peak_and_hot_nights(self):
        dates = pd.date_range('2019-08-01', '2020-12-31', freq='D')
        sst = np.full(len(dates), 28.0)
        spell = (dates >= '2020-02-01') & (dates <= '2020-02-10')
        sst[spell] = 30.5
        result = summarise_site_year(dates, sst, 28.0, 2020)
        self.assertEqual(result['imos_hotspell_mmm1_max_days'], 10)
        self.assertAlmostEqual(result['imos_hotspot_3d_max'], 2.5)
        self.assertAlmostEqual(result['imos_extreme_dhd_mmm2'], 0.5 * 10)
        self.assertEqual(result['imos_hot_nights_mmm1_count'], 10)

    def test_monthly_threshold_interpolation_is_periodic(self):
        monthly = np.arange(1, 13, dtype=float)
        dates = pd.DatetimeIndex(['2005-01-01', '2005-07-16', '2005-12-31'])
        daily = interpolate_monthly_climatology(monthly, dates)
        self.assertEqual(daily.shape, (3,))
        self.assertGreater(daily[0], 1.0)
        self.assertLess(daily[0], 12.0)
        self.assertAlmostEqual(daily[1], 7.0, places=1)
        self.assertGreater(daily[2], 1.0)

    def test_hobday_events_require_five_days_and_join_two_day_gap(self):
        sst = np.zeros(20)
        sst[1:6] = 2.0
        sst[8:13] = 3.0
        sst[15:19] = 4.0  # four days: not a separate event
        result = summarise_marine_heatwaves(
            sst, np.zeros(20), np.zeros(20), np.ones(20)
        )
        self.assertEqual(result['imos_mhw_event_count'], 1)
        self.assertEqual(result['imos_mhw_max_duration_days'], 12)
        self.assertEqual(result['imos_mhw_max_intensity_c'], 3.0)
        self.assertEqual(result['imos_mhw_category2_days'], 10)

    def test_percentile_heat_dose_counts_isolated_extreme_days(self):
        dates = pd.date_range('2019-11-01', '2020-04-30', freq='D')
        peak = pd.Timestamp('2020-02-15')
        offset = np.abs((dates - peak).days.to_numpy())
        climatology = 30.0 - offset / 365.0
        threshold90 = climatology + 1.0
        sst = climatology.copy()
        for delta in (-10, 0, 10):
            sst[dates == peak + pd.Timedelta(days=delta)] += 2.0
        result = summarise_percentile_heat_dose(
            dates, sst, climatology, threshold90, 2020
        )
        self.assertEqual(
            result['imos_pbd_peak_climatology_date'], '2020-02-15'
        )
        self.assertAlmostEqual(result['imos_pbd12_c_weeks'], 6 / 7)
        self.assertEqual(result['imos_pbd12_exceedance_days'], 3)
        self.assertAlmostEqual(result['imos_pbd12_coverage'], 1.0)

    def test_daytime_pairs_with_following_local_night(self):
        day_dates = pd.date_range('2019-11-01', periods=6, freq='D')
        night_dates = day_dates + pd.Timedelta(days=1)
        result = summarise_day_night_relief(
            day_dates, np.full(6, 30.0),
            night_dates, np.array([28.5, 29.5, 28.0, 30.0, 28.2, 28.9]),
            28.0, 2020,
        )
        self.assertEqual(result['imos_hot_day_night_pair_count'], 6)
        self.assertAlmostEqual(result['imos_night_relief_fraction'], 4 / 6)
        self.assertAlmostEqual(
            result['imos_day_to_following_night_drop_c'], 1.3
        )

    def test_intensity_duration_profile_distinguishes_persistence(self):
        sharp = np.zeros(100)
        sharp[:3] = 3.0
        persistent = np.full(100, 1.5)
        sharp_result = summarise_intensity_duration_profile(sharp)
        persistent_result = summarise_intensity_duration_profile(persistent)
        self.assertAlmostEqual(
            persistent_result[
                'imos_id_hotspot_persistence_slope_c_per_log_day'
            ], 0.0, places=10
        )
        self.assertLess(
            sharp_result['imos_id_hotspot_persistence_slope_c_per_log_day'],
            -0.3,
        )


if __name__ == '__main__':
    unittest.main()
