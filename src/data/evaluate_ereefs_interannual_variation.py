"""
Evaluate interannual variability of eReefs current speeds across GBR reefs
during peak thermal stress / bleaching windows (Jan-Mar).
"""
import numpy as np
import pandas as pd
import s3fs
import zarr

def evaluate_ereefs_variability():
    print("Opening eReefs daily hydrodynamic Zarr store...")
    fs = s3fs.S3FileSystem(anon=True)
    store = s3fs.S3Map('gbr-dms-data-public/aims-ereefs-agg-hydrodynamic-4km-daily/data.zarr', s3=fs)
    zg = zarr.open_group(store, mode='r')
    
    # Load coordinates
    t_raw = zg['time'][:]
    times = pd.to_datetime(t_raw, unit='D', origin='1990-01-01')
    lats = zg['latitude'][:]
    lons = zg['longitude'][:]
    
    years = [2016, 2017, 2020, 2021, 2022, 2023]
    window_means = {}
    
    lat_indices = np.linspace(50, len(lats) - 50, 15, dtype=int)
    lon_indices = np.linspace(50, len(lons) - 50, 15, dtype=int)
    
    print(f"Sample grid points: {len(lat_indices)} x {len(lon_indices)} = {len(lat_indices)*len(lon_indices)}")
    
    for yr in years:
        start_date = f"{yr}-01-01"
        end_date = f"{yr}-03-31"
        mask = (times >= start_date) & (times <= end_date)
        t_indices = np.where(mask)[0]
        
        if len(t_indices) == 0:
            continue
            
        print(f"Extracting {yr} (N={len(t_indices)} days)...")
        t_min, t_max = t_indices[0], t_indices[-1] + 1
        
        # Read block for selected lat/lon sample
        vals = []
        for i_lat in lat_indices:
            for i_lon in lon_indices:
                pt_series = zg['mean_cur'][t_min:t_max, 0, i_lat, i_lon]
                vals.append(np.nanmean(pt_series))
                
        window_means[yr] = np.array(vals)
        
    df_win = pd.DataFrame(window_means).dropna()
    print("\n=== Interannual Correlation Matrix (Jan-Mar Mean Current Speed) ===")
    corr_matrix = df_win.corr()
    print(corr_matrix.round(3))
    
    print("\n=== Summary Statistics Across Sample Points ===")
    overall_mean = df_win.values.mean()
    inter_year_std = df_win.std(axis=1).mean()
    print(f"Mean current speed across GBR sample: {overall_mean:.4f} m/s")
    print(f"Mean inter-annual std per location: {inter_year_std:.4f} m/s")
    print(f"Relative variability (std/mean): {(inter_year_std/overall_mean)*100:.2f}%")

if __name__ == "__main__":
    evaluate_ereefs_variability()

