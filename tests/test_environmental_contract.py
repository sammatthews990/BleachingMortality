import unittest

import pandas as pd

from scripts.validate_environmental_data import (
    EXPECTED_YEARS,
    validate_environmental_dataframe,
)


def valid_frame() -> pd.DataFrame:
    days = {2016: 70, 2017: 89, 2020: 89, 2022: 82, 2024: 82}
    preliminary = {2016: 0, 2017: 0, 2020: 176, 2022: 0, 2024: 164}
    rows = []
    for year in EXPECTED_YEARS:
        rows.append(
            {
                "year": year,
                "LABEL_ID": "reef-1",
                "LOC_NAME_S": "Test Reef",
                "lon": 145.0,
                "lat": -18.0,
                "ann_maxdhw": 8.0,
                "histmDHW6": 7.0,
                "yrsince6": 3,
                "histmDHW4": 5.0,
                "yrsince4": 2,
                "ann_maxsst": 30.0,
                "winyear_mean": 25.0,
                "winyear_sd": 1.5,
                "mcur_90": 0.4,
                "dist_to_er_km": 2.0,
                "secc3m": 15.0,
                "cloudp_90": 0.7,
                "cloud_n_days": days[year],
                "cloud_n_asc": days[year],
                "cloud_n_des": days[year],
                "cloud_n_paired_days": days[year],
                "cloud_n_preliminary_files": preliminary[year],
                "cloud_platform": "NOAA-18",
                "cloud_product_version": "v06r00",
                "cloud_aggregation": "mean(Jan-Mar NOAA-18 ascending pass)",
                "cloud_spatial_match": "paper convention",
            }
        )
    frame = pd.DataFrame(rows)
    frame['dhw10_mean'] = 3.0
    frame['dhw10_max'] = 7.0
    frame['dhw10_load4'] = 8.0
    frame['dhw10_n4'] = 3
    frame['dhw10_n6'] = 1
    frame['dhw_novelty10'] = 1.0
    frame['k490_q90'] = 0.2
    frame['secc3m_p10'] = 8.5
    return frame


class EnvironmentalContractTests(unittest.TestCase):
    def test_valid_table_passes(self):
        report = validate_environmental_dataframe(valid_frame())
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["reef_features"], 1)

    def test_duplicate_feature_year_fails(self):
        data = valid_frame()
        data = pd.concat([data, data.iloc[[0]]], ignore_index=True)
        report = validate_environmental_dataframe(data)
        self.assertFalse(report["valid"])
        self.assertTrue(any("duplicate reef-feature-year" in error for error in report["errors"]))

    def test_out_of_range_cloud_fails(self):
        data = valid_frame()
        data.loc[0, "cloudp_90"] = 1.2
        report = validate_environmental_dataframe(data)
        self.assertFalse(report["valid"])
        self.assertTrue(any("cloudp_90" in error for error in report["errors"]))


if __name__ == "__main__":
    unittest.main()
