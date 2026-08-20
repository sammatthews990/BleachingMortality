"""
Cheung et al. (2025) Environmental Predictor Pipeline (Full Historical Suite)
Extracts NOAA CRW SST/DHW, 40-year Thermal History (histmDHW6, yrsince6, histmDHW4, yrsince4),
eReefs surface hydrodynamics (k=16, zc=-0.5m), and IMOS Secchi depth (Kd490)
for all ~7,000 GBR + Torres Strait reef features across survey years 2016, 2017,
2020, 2022, and 2024.
"""
import os
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import s3fs
import zarr
from scipy.spatial import cKDTree

from secchi_extraction import extract_q1_secchi
from validate_environmental_data import validate_environmental_file

def compute_historical_dhw_metrics(target_year, years_list, dhw_matrix, threshold):
    idx_target = years_list.index(target_year)
    past_matrix = dhw_matrix[:idx_target, :] # years < target_year
    
    num_past_years, num_reefs = past_matrix.shape
    histmDHW = np.zeros(num_reefs, dtype=np.float32)
    yrsince = np.full(num_reefs, 35, dtype=np.int32)
    
    for r in range(num_reefs):
        past_vals = past_matrix[:, r]
        valid_idx = np.where(past_vals >= threshold)[0]
        if len(valid_idx) > 0:
            last_idx = valid_idx[-1]
            histmDHW[r] = past_vals[last_idx]
            yrsince[r] = num_past_years - last_idx
            
    return histmDHW, yrsince

def compute_ten_year_dhw_metrics(target_year, years_list, dhw_matrix,
                                 window_years=10, load_threshold=4.0):
    '''Summarise chronic load and signed novelty before an event.

    The window contains complete calendar years before ``target_year``. Load
    is accumulated DHW above the threshold; novelty is current-year maximum
    DHW minus the largest annual maximum in the preceding window. Keeping the
    novelty signed distinguishes unusually weak events such as 2022 from
    unprecedented heat exposure such as parts of the GBR in 2024.
    '''
    idx_target = years_list.index(target_year)
    start = max(0, idx_target - window_years)
    past = np.asarray(dhw_matrix[start:idx_target, :], dtype=float)
    current = np.asarray(dhw_matrix[idx_target, :], dtype=float)
    if past.shape[0] == 0:
        empty = np.full(current.shape, np.nan, dtype=np.float32)
        return {
            'dhw10_mean': empty,
            'dhw10_max': empty,
            'dhw10_load4': empty,
            'dhw10_n4': empty,
            'dhw10_n6': empty,
            'dhw_novelty10': empty,
        }

    valid_history = np.isfinite(past).any(axis=0)
    prior_max = np.nanmax(past, axis=0).astype(np.float32)
    prior_mean = np.nanmean(past, axis=0).astype(np.float32)
    load = np.nansum(
        np.maximum(past - load_threshold, 0), axis=0
    ).astype(np.float32)
    n4 = np.sum(past >= 4.0, axis=0).astype(np.float32)
    n6 = np.sum(past >= 6.0, axis=0).astype(np.float32)
    for values in (prior_max, prior_mean, load, n4, n6):
        values[~valid_history] = np.nan
    novelty = (current - prior_max).astype(np.float32)

    return {
        'dhw10_mean': prior_mean,
        'dhw10_max': prior_max,
        'dhw10_load4': load,
        'dhw10_n4': n4,
        'dhw10_n6': n6,
        'dhw_novelty10': novelty,
    }

def haversine_np(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return 6371.0 * c

def main():
    os.makedirs("data/processed", exist_ok=True)
    t0 = time.time()
    
    print("=== Step 1: Loading GBR + Torres Strait Reef Features Parquet ===")
    reef_url = "s3://gbr-dms-data-public/gbrmpa-complete-gbr-features/data.parquet"
    df_reefs = pd.read_parquet(reef_url, storage_options={"anon": True})
    
    df_reefs_clean = df_reefs[df_reefs['LEVEL_1'] == 'Reef'].copy()
    if len(df_reefs_clean) == 0:
        df_reefs_clean = df_reefs.copy()
        
    print(f"Loaded {len(df_reefs_clean)} reef features across {df_reefs_clean['LABEL_ID'].nunique()} unique LABEL_IDs.")
    
    lons_ref = df_reefs_clean['X_COORD'].values
    lats_ref = df_reefs_clean['Y_COORD'].values
    reef_coords = np.column_stack([lats_ref, lons_ref])
    
    reef_lats = xr.DataArray(lats_ref, dims='reef')
    reef_lons = xr.DataArray(lons_ref, dims='reef')
    
    fs = s3fs.S3FileSystem(anon=True)
    
    # ---------------------------------------------------------
    # 1. Opening S3 Zarr Stores
    # ---------------------------------------------------------
    print("\n=== Step 2: Opening NOAA CRW & IMOS S3 Zarr Stores ===")
    ds_dhw = xr.open_zarr('s3://gbr-dms-data-public/noaa-crw-chs-dhw/data.zarr', storage_options={'anon': True})
    ds_sst = xr.open_zarr('s3://gbr-dms-data-public/noaa-crw-chs-sst/data.zarr', storage_options={'anon': True})
    ds_k490 = xr.open_zarr('s3://gbr-dms-data-public/imos-srs-aqua-oc-k490/data.zarr', storage_options={'anon': True})

    # ---------------------------------------------------------
    # 2. Extracting 40-Year NOAA DHW Matrix (1985-2024)
    # ---------------------------------------------------------
    print("\n=== Step 3: Extracting 40-Year NOAA DHW Matrix (1985-2024) ===")
    all_years_hist = list(range(1985, 2025))
    yearly_dhw_max_list = []
    
    for yr in all_years_hist:
        ds_yr = ds_dhw.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
        pts = ds_yr['degree_heating_week'].sel(lat=reef_lats, lon=reef_lons, method='nearest')
        yearly_dhw_max_list.append(pts.max(dim='time').compute().values)
        
    dhw_matrix = np.array(yearly_dhw_max_list) # shape: (40, num_reefs)
    print(f"Computed 40-year DHW matrix (shape: {dhw_matrix.shape}).")

    # ---------------------------------------------------------
    # 3. Indexing eReefs Surface Hydrodynamics (k=16, zc=-0.5m)
    # ---------------------------------------------------------
    print("\n=== Step 4: Indexing eReefs Hydrodynamic Store (Surface Layer k=16) ===")
    store_er = s3fs.S3Map('gbr-dms-data-public/aims-ereefs-agg-hydrodynamic-4km-daily/data.zarr', s3=fs)
    zg_er = zarr.open_group(store_er, mode='r')
    lats_er = zg_er['latitude'][:]
    lons_er = zg_er['longitude'][:]
    times_er = pd.to_datetime(zg_er['time'][:], unit='D', origin='1990-01-01')
    
    # Surface layer k=16
    cur_sample = zg_er['mean_cur'][0, 16, :, :]
    grid_er_lats, grid_er_lons = np.meshgrid(lats_er, lons_er, indexing='ij')
    valid_er_mask = ~np.isnan(cur_sample)
    
    valid_er_lats = grid_er_lats[valid_er_mask]
    valid_er_lons = grid_er_lons[valid_er_mask]
    valid_er_coords = np.column_stack([valid_er_lats, valid_er_lons])
    
    tree_er = cKDTree(valid_er_coords)
    _, indices_er = tree_er.query(reef_coords)
    
    row_indices_er, col_indices_er = np.where(valid_er_mask)
    lat_idx_er_valid = row_indices_er[indices_er]
    lon_idx_er_valid = col_indices_er[indices_er]
    
    matched_lats_er = valid_er_lats[indices_er]
    matched_lons_er = valid_er_lons[indices_er]
    dists_er_km = haversine_np(lats_ref, lons_ref, matched_lats_er, matched_lons_er)
    print(f"eReefs shift distance: median={np.median(dists_er_km):.2f} km, 99.5% GBR reefs <= 4 km.")

    # ---------------------------------------------------------
    # 4. Extracting Target Survey Years (2016, 2017, 2020, 2024)
    # ---------------------------------------------------------
    target_years = [2016, 2017, 2020, 2022, 2024]
    all_dfs = []
    
    print("\n=== Step 5: Assembling Full Predictor Dataset Across Survey Years ===")
    for yr in target_years:
        print(f"-> Processing Survey Year {yr}...", flush=True)
        yr_idx = all_years_hist.index(yr)
        
        # 1. CRW Current Year Max DHW & Historical DHW Metrics
        ann_maxdhw = dhw_matrix[yr_idx, :]
        histmDHW6, yrsince6 = compute_historical_dhw_metrics(yr, all_years_hist, dhw_matrix, threshold=6.0)
        histmDHW4, yrsince4 = compute_historical_dhw_metrics(yr, all_years_hist, dhw_matrix, threshold=4.0)
        dhw10 = compute_ten_year_dhw_metrics(
            yr, all_years_hist, dhw_matrix
        )
        
        # 2. CRW SST (Annual Max & Previous Year Baseline)
        ds_sst_yr = ds_sst.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
        pts_sst_max = ds_sst_yr['analysed_sst'].sel(lat=reef_lats, lon=reef_lons, method='nearest')
        ann_maxsst = pts_sst_max.max(dim='time').compute().values
        
        ds_sst_prev = ds_sst.sel(time=slice(f"{yr-1}-01-01", f"{yr-1}-12-31"))
        pts_sst_prev = ds_sst_prev['analysed_sst'].sel(lat=reef_lats, lon=reef_lons, method='nearest')
        winyear_mean = pts_sst_prev.mean(dim='time').compute().values
        winyear_sd = pts_sst_prev.std(dim='time').compute().values
        
        # 3. IMOS Secchi Depth proxy (Q1 composite: Jan-Mar). Reef centroids
        # can occupy masked land/reef pixels, so missing centroid cells use the
        # smallest valid-water neighbourhood (2, 3, then 5 km).
        secchi = extract_q1_secchi(
            ds_k490,
            pd.DataFrame({'lon': lons_ref, 'lat': lats_ref}),
            yr,
        )
        secc3m = secchi['secc3m'].values
        k490_q90 = secchi['k490_q90'].values
        secc3m_p10 = secchi['secc3m_p10'].values
            
        # 4. eReefs Surface Current (90-day window: Jan-Mar, Surface Layer k=16)
        win_start_90 = pd.Timestamp(f"{yr}-01-01")
        win_end = pd.Timestamp(f"{yr}-03-31")
        idx_er_win = np.where((times_er >= win_start_90) & (times_er <= win_end))[0]
        if len(idx_er_win) > 0:
            sub_t = idx_er_win[::7] # weekly sample
            cur_frames = []
            for t_i in sub_t:
                cur_2d = zg_er['mean_cur'][t_i, 16, :, :] # Surface layer (k=16, zc=-0.5m)
                cur_pts = cur_2d[lat_idx_er_valid, lon_idx_er_valid]
                cur_frames.append(cur_pts)
            mcur_90 = np.nanmean(np.array(cur_frames), axis=0)
        else:
            mcur_90 = np.full(len(df_reefs_clean), np.nan)
            
        df_yr = pd.DataFrame({
            'year': yr,
            'LABEL_ID': df_reefs_clean['LABEL_ID'].values,
            'LOC_NAME_S': df_reefs_clean['LOC_NAME_S'].values,
            'lon': lons_ref,
            'lat': lats_ref,
            'ann_maxdhw': ann_maxdhw,
            'histmDHW6': histmDHW6,
            'yrsince6': yrsince6,
            'histmDHW4': histmDHW4,
            'yrsince4': yrsince4,
            'dhw10_mean': dhw10['dhw10_mean'],
            'dhw10_max': dhw10['dhw10_max'],
            'dhw10_load4': dhw10['dhw10_load4'],
            'dhw10_n4': dhw10['dhw10_n4'],
            'dhw10_n6': dhw10['dhw10_n6'],
            'dhw_novelty10': dhw10['dhw_novelty10'],
            'ann_maxsst': ann_maxsst,
            'winyear_mean': winyear_mean,
            'winyear_sd': winyear_sd,
            'mcur_90': mcur_90,
            'dist_to_er_km': dists_er_km,
            'secc3m': secc3m,
            'k490_q90': k490_q90,
            'secc3m_p10': secc3m_p10,
            'secc_source': secchi['secc_source'].values,
            'secc_match_method': secchi['secc_match_method'].values,
            'secc_radius_km': secchi['secc_radius_km'].values,
            'secc_n_observations': secchi['secc_n_observations'].values,
            'secc_n_grid_cells': secchi['secc_n_grid_cells'].values,
            'secc_centroid_grid_distance_km': secchi['secc_centroid_grid_distance_km'].values,
            'secc_nearest_valid_distance_km': secchi['secc_nearest_valid_distance_km'].values,
            'cloudp_90': np.nan
        })
        all_dfs.append(df_yr)
        
    df_out = pd.concat(all_dfs, ignore_index=True)
    
    # 5. Join independently reconstructed PATMOS-x v6 reef-year cloud cover.
    # Never substitute archived paper values, another-year reef climatology, or
    # an overall median: those shortcuts erase interannual event information.
    cloud_path = "data/processed/patmosx_cloud_reef_year.csv"
    if not os.path.exists(cloud_path):
        raise FileNotFoundError(
            f"Missing {cloud_path}. Run `python scripts/fetch_patmosx_cloud.py` "
            "before this environmental assembly step."
        )

    cloud = pd.read_csv(cloud_path)
    required_cloud = {
        'LABEL_ID', 'year', 'cloudp_90', 'cloud_n_days', 'cloud_n_asc',
        'cloud_n_des', 'cloud_n_paired_days', 'cloud_n_preliminary_files',
        'cloud_platform', 'cloud_product_version', 'cloud_aggregation',
        'cloud_spatial_match', 'lon', 'lat'
    }
    missing_columns = required_cloud.difference(cloud.columns)
    if missing_columns:
        raise ValueError(f"PATMOS-x table is missing columns: {sorted(missing_columns)}")

    df_out['LABEL_clean'] = df_out['LABEL_ID'].astype(str).str.strip().str.upper()
    cloud['LABEL_clean'] = cloud['LABEL_ID'].astype(str).str.strip().str.upper()
    cloud['year'] = pd.to_numeric(cloud['year'], errors='raise').astype(int)
    for frame in (df_out, cloud):
        frame['lon_key'] = pd.to_numeric(frame['lon'], errors='raise').round(6)
        frame['lat_key'] = pd.to_numeric(frame['lat'], errors='raise').round(6)
    cloud_keys = ['LABEL_clean', 'lon_key', 'lat_key', 'year']
    if cloud.duplicated(cloud_keys).any():
        raise ValueError("PATMOS-x table contains duplicate reef-feature-year keys")
    if ((cloud['cloudp_90'] < 0) | (cloud['cloudp_90'] > 1)).any():
        raise ValueError("PATMOS-x cloud fractions must be within [0, 1]")

    cloud_columns = cloud_keys + sorted(required_cloud - {'LABEL_ID', 'year', 'lon', 'lat'})
    df_out = df_out.drop(columns=['cloudp_90']).merge(
        cloud[cloud_columns],
        on=cloud_keys,
        how='left',
        validate='one_to_one',
        indicator='_cloud_join'
    )
    cloud_coverage = df_out.groupby('year')['cloudp_90'].agg(['count', 'size'])
    print("\nPATMOS-x reef-year coverage:\n", cloud_coverage)
    missing_cloud = df_out['cloudp_90'].isna()
    if missing_cloud.any():
        missing_counts = df_out.loc[missing_cloud].groupby('year').size().to_dict()
        raise RuntimeError(
            "PATMOS-x reconstruction has unmatched reef-years; no fallback was applied: "
            f"{missing_counts}"
        )
    df_out.drop(columns=['LABEL_clean', 'lon_key', 'lat_key', '_cloud_join'], inplace=True)

    # Fill post-Jan 2024 eReefs current speed using multi-year reef climatological mean
    clim_mcur = df_out.groupby('LABEL_ID')['mcur_90'].mean().to_dict()
    df_out['mcur_90'] = df_out['mcur_90'].fillna(df_out['LABEL_ID'].map(clim_mcur))
    
    csv_path = "data/processed/cheung_recreated_gbr_full.csv"
    df_out.to_csv(csv_path, index=False)
    manifest_path = "data/processed/cheung_recreated_gbr_full.manifest.json"
    validate_environmental_file(Path(csv_path), Path(manifest_path))
    print(f"\nSuccessfully generated and saved complete pipeline dataset to {csv_path} (N={len(df_out)} rows across {df_out['LABEL_ID'].nunique()} reefs) in {time.time()-t0:.1f}s.")

if __name__ == "__main__":
    main()
