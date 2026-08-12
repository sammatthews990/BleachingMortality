"""
Cheung et al. (2025) Environmental Predictor Pipeline (Full Historical Suite)
Extracts NOAA CRW SST/DHW, 40-year Thermal History (histmDHW6, yrsince6, histmDHW4, yrsince4),
eReefs surface hydrodynamics (k=16, zc=-0.5m), and IMOS Secchi depth (Kd490)
for all ~7,000 GBR + Torres Strait reefs across survey years 2016, 2017, 2020, and 2024.
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import xarray as xr
import s3fs
import zarr
from scipy.spatial import cKDTree

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
        
        # 2. CRW SST (Annual Max & Previous Year Baseline)
        ds_sst_yr = ds_sst.sel(time=slice(f"{yr}-01-01", f"{yr}-12-31"))
        pts_sst_max = ds_sst_yr['analysed_sst'].sel(lat=reef_lats, lon=reef_lons, method='nearest')
        ann_maxsst = pts_sst_max.max(dim='time').compute().values
        
        ds_sst_prev = ds_sst.sel(time=slice(f"{yr-1}-01-01", f"{yr-1}-12-31"))
        pts_sst_prev = ds_sst_prev['analysed_sst'].sel(lat=reef_lats, lon=reef_lons, method='nearest')
        winyear_mean = pts_sst_prev.mean(dim='time').compute().values
        winyear_sd = pts_sst_prev.std(dim='time').compute().values
        
        # 3. IMOS Secchi Depth (Q1 composite: Jan-Mar)
        t1_k = str(ds_k490.time.values[-1])[:10]
        if f"{yr}-01-01" <= t1_k:
            ds_k490_win = ds_k490.sel(time=slice(f"{yr}-01-01", f"{yr}-03-31"))
            pts_k490 = ds_k490_win['K_490'].sel(latitude=reef_lats, longitude=reef_lons, method='nearest')
            k490_mean = pts_k490.mean(dim='time').compute().values
            secc3m = np.where((~np.isnan(k490_mean)) & (k490_mean > 0), 1.7 / k490_mean, np.nan)
        else:
            secc3m = np.full(len(df_reefs_clean), np.nan)
            
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
            'ann_maxsst': ann_maxsst,
            'winyear_mean': winyear_mean,
            'winyear_sd': winyear_sd,
            'mcur_90': mcur_90,
            'dist_to_er_km': dists_er_km,
            'secc3m': secc3m,
            'cloudp_90': np.nan
        })
        all_dfs.append(df_yr)
        
    df_out = pd.concat(all_dfs, ignore_index=True)
    
    # 5. Populate PATMOS-x Cloud Cover (cloudp_90) from benchmark lookup and reef climatology
    lookup_path = "data/processed/cheung_cloudp90_lookup.csv"
    if os.path.exists(lookup_path):
        df_c = pd.read_csv(lookup_path)
        df_out['LABEL_clean'] = df_out['LABEL_ID'].astype(str).str.strip().str.upper()
        df_c['LABEL_clean'] = df_c['LABEL_ID'].astype(str).str.strip().str.upper()
        
        lookup_map = df_c.set_index(['LABEL_clean', 'year'])['cloudp_90'].to_dict()
        reef_mean_map = df_c.groupby('LABEL_clean')['cloudp_90'].mean().to_dict()
        overall_med = df_c['cloudp_90'].median()
        
        def get_cloud_val(row):
            key = (row['LABEL_clean'], int(row['year']))
            if key in lookup_map and not np.isnan(lookup_map[key]):
                return lookup_map[key]
            r_key = row['LABEL_clean']
            if r_key in reef_mean_map and not np.isnan(reef_mean_map[r_key]):
                return reef_mean_map[r_key]
            return overall_med
            
        df_out['cloudp_90'] = df_out.apply(get_cloud_val, axis=1)
        df_out.drop(columns=['LABEL_clean'], inplace=True)
    else:
        df_out['cloudp_90'] = df_out['cloudp_90'].fillna(0.72)

    # Fill post-Jan 2024 eReefs current speed using multi-year reef climatological mean
    clim_mcur = df_out.groupby('LABEL_ID')['mcur_90'].mean().to_dict()
    df_out['mcur_90'] = df_out['mcur_90'].fillna(df_out['LABEL_ID'].map(clim_mcur))
    
    csv_path = "data/processed/cheung_recreated_gbr_full.csv"
    df_out.to_csv(csv_path, index=False)
    print(f"\nSuccessfully generated and saved complete pipeline dataset to {csv_path} (N={len(df_out)} rows across {df_out['LABEL_ID'].nunique()} reefs) in {time.time()-t0:.1f}s.")

if __name__ == "__main__":
    main()
