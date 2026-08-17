# BleachingMortality

## Reconstructing PATMOS-x cloud cover

Install the isolated cloud-pipeline dependencies and reconstruct all five
bleaching-event years before assembling the environmental table:

```powershell
python -m pip install -r scripts/requirements-patmosx.txt
python scripts/fetch_patmosx_cloud.py
python scripts/fetch_dms_environmental_data.py
```

The cloud step streams only GBR chunks from NOAA's public PATMOS-x v6 archive,
caches them under `data/cache/patmosx_s3_gbr`, and writes the reef-year predictor
and provenance metadata under `data/processed`. It does not use the archived
Cheung cloud values as inputs or fill later years from a historical climatology.

See `reports/research_pipeline_implementation_plan.md` for the ordered pipeline
improvements and cloud validation summary.
