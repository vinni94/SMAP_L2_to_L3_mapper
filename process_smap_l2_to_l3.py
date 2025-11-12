# Imports
from __future__ import division
import argparse
from pathlib import Path
import numpy as np
import xarray as xr
import h5py
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from datetime import datetime
import pandas as pd
import gc
from tqdm import tqdm
import re
from numpy import cos, sin, tan, mod, sqrt, all

# ---------------------------
# ARGUMENT PARSER
# ---------------------------
parser = argparse.ArgumentParser(description="SMAP Processing Pipeline (3km/9km)")

parser.add_argument("--res", type=int, choices=[3, 9], required=True,
                    help="Spatial resolution: 3 for 3km L2, 9 for 9km L3")
parser.add_argument("--year", type=int, required=True, help="Year to process")
parser.add_argument("--lat_start", type=float, required=True, help="Northern latitude (top)")
parser.add_argument("--lat_end", type=float, required=True, help="Southern latitude (bottom)")
parser.add_argument("--lon_start", type=float, required=True, help="Western longitude (left)")
parser.add_argument("--lon_end", type=float, required=True, help="Eastern longitude (right)")
parser.add_argument("--output_root", type=str, required=True, help="Root folder for output NetCDFs")
parser.add_argument("--input_dir", type=str, required=True, help="Input directory with SMAP HDF5 files")
parser.add_argument("--region", type=str, required=True, help="Give a shortname of region of interest")
parser.add_argument("--easegrid_9km_file", type=str, required=True, help="Give path of easegrid global lat lon file")
parser.add_argument("--easegrid_3km_file", type=str, required=True, help="Give path of easegrid global lat lon file")

args = parser.parse_args()

# ---------------------------
# CONFIGURATION
# ---------------------------
res = args.res
year = args.year
lat_start, lat_end = args.lat_start, args.lat_end
lon_start, lon_end = args.lon_start, args.lon_end
fdir = Path(args.input_dir)
output_root = Path(args.output_root)
output_root.mkdir(exist_ok=True, parents=True)
region = args.region

#############################
# Read Easegrid Global lat lon file
#############################
if res == 3:
    easegrid_latlon_3km = xr.open_dataset(args.easegrid_3km_file)
elif res == 9:
    easegrid_latlon_9km = xr.open_dataset(args.easegrid_9km_file)
    
##############################
# The EASE2 helper functions
##############################

def easeconv_normalize_degrees(dlon):
    #
    # Return dlon to within the range -180 <= dlon <= 180
    #  can handle array inputs
    #
    out = dlon
    if np.size(out) > 1:
        while (out < -180.0).sum() > 0:
            out[out < -180.0] = out[out > 180.0] + 360.0
        while (out > 180.0).sum() > 0:
            out[out > 180.0] = out[out > 180.0] - 360.0
    else:
        while out < -180.0:
            out = out + 360.0
        while out > 180.0:
            out = out - 360.0
    return out


def ease2_map_info(iopt, isc, ind):
    """ internally used routine

 (map_equatorial_radius_m, map_eccentricity, \
  e2, map_reference_latitude, map_reference_longitude, \
  map_second_reference_latitude, sin_phi1, cos_phi1, kz, \
  map_scale, bcols, brows, r0, s0, epsilon) = ease2_map_info(iopt, isc, nd)

       defines EASE2 grid information

	inputs:
         iopt: projection type 8=EASE2 N, 9=EASE2 S, 10=EASE2 T/M
         isc:  scale factor 0..5 grid size is (basesize(ind))/2^isc
         ind:  base grid size index   (map units per cell in m

        NSIDC .grd file for isc=0
           project type    ind=0     ind=1         ind=2       ind=3
	      N         EASE2_N25km EASE2_N30km EASE2_N36km  EASE2_N24km
              S         EASE2_S25km EASE2_S30km EASE2_S36km  EASE2_S24km
              T/M       EASE2_T25km EASE2_M25km EASE2_M36km  EASE2_M24km

          cell size (m) for isc=0 (scale is reduced by 2^isc)
           project type  ind=0     ind=1            ind=2        ind=3
	      N        25000.0     30000.0         36000.0      24000.0
              S        25000.0     30000.0         36000.0      24000.0
              T/M   T25025.26 M25025.26000 M36032.220840584 M24021.480560389347

	  for a given base cell size (e.g., ind=0) isc is related to
          NSIDC CETB .grd file names according to
	     isc        N .grd name   S .grd name   T .grd name
	      0	      EASE2_N25km     EASE2_S25km     EASE2_T25km
	      1	      EASE2_N12.5km   EASE2_S12.5km   EASE2_T12.5km
	      2	      EASE2_N6.25km   EASE2_S6.25km   EASE2_T6.25km
	      3	      EASE2_N3.125km  EASE2_S3.125km  EASE2_T3.125km
	      4	      EASE2_N1.5625km EASE2_S1.5625km EASE2_T1.5625km

  outputs
    map_equatorial_radius_m  EASE2 Earth equitorial radius (km) [WGS84]
    map_eccentricity         EASE2 Earth eccentricity [WGS84]
    map_reference_latitude   Reference latitude (deg)
    map_reference_longitude  Reference longitude (deg)
    map_second_reference_latitude Secondary reference longitude* (deg)
    sin_phi1, cos_phi1 kz    EASE2 Cylin parameters*
    map_scale                EASE2 map projection pixel size (km)
    bcols, brows,            EASE2 grid size in pixels
    r0, s0                   EASE2 base projection size in pixels
    epsilon                  EASE2 near-polar test factor
    """

    DTR = 0.01745329241994

    m = 2 ** np.floor(isc)  # compute power-law scale factor

    map_equatorial_radius_m = 6378137.0  # WGS84
    map_eccentricity = 0.081819190843  # WGS84
    e2 = map_eccentricity * map_eccentricity
    map_reference_longitude = 0.0
    epsilon = 1.e-6

    #  map-specific parameters
    if iopt == 8:  # EASE2 grid north
        map_reference_latitude = 90.0
        if ind == 1:  # EASE2_N30km.gpd
            base = 30000.0
            nx = 600
            ny = 600
        elif ind == 2:  # EASE2_N36km.gpd
            base = 36000.0
            nx = 500
            ny = 500
        elif ind == 3:  # EASE2_N24km.gpd
            base = 24000.0
            nx = 750
            ny = 750
        else:  # EASE2_N25km.gpd
            base = 25000.0
            nx = 720
            ny = 720
        map_second_reference_latitude = 0.0
        sin_phi1 = 0.0
        cos_phi1 = 1.0
        kz = cos_phi1

    elif iopt == 9:  # EASE2 grid south
        map_reference_latitude = -90.0
        if ind == 1:  # EASE2_S30km.gpd
            base = 30000.0
            nx = 600
            ny = 600
        elif ind == 2:  # EASE2_S36km.gpd
            base = 36000.0
            nx = 500
            ny = 500
        elif ind == 3:  # EASE2_S24km.gpd
            base = 24000.0
            nx = 750
            ny = 750
        else:  # EASE2_S25km.gpd
            base = 25000.0
            nx = 720
            ny = 720
        map_second_reference_latitude = 0.0
        sin_phi1 = 0.0
        cos_phi1 = 1.0
        kz = cos_phi1

    elif iopt == 10:  # EASE2 cylindrical
        map_reference_latitude = 0.0
        map_second_reference_latitude = 30.0
        sin_phi1 = np.sin(DTR * map_second_reference_latitude)
        cos_phi1 = np.cos(DTR * map_second_reference_latitude)
        kz = cos_phi1 / np.sqrt(1.0 - e2 * sin_phi1 * sin_phi1)
        if ind == 1:  # EASE2_M25km.gpd
            base = 25025.26000
            nx = 1388
            ny = 584
        elif ind == 2:  # EASE2_M36km.gpd
            base = 36032.220840584
            nx = 964
            ny = 406
        elif ind == 3:  # EASE2_M24km.gpd
            base = 24021.480560389347
            nx = 1446
            ny = 609
        elif ind == 4:  # EASE2_M09km.gpd
            base = 9008.055210146
            nx = 3856
            ny = 1624
        elif ind == 5:  # EASE2_M03km.gpd
            print('Entered bscale 5 in ease2_map_info')
            base = 3002.6850700487
            nx = 11568
            ny = 4872
        else:  # EASE2_T25km.gpd
            base = 25025.26000
            nx = 1388
            ny = 540

    else:
        print('*** invalid EASE2 projection code ***')

    # grid info
    if isc >= 0:
        map_scale = base / m
        bcols = np.ceil(nx * m)
        brows = np.ceil(ny * m)
        r0 = (nx * m - 1) / 2
        s0 = (ny * m - 1) / 2
    else:
        map_scale = base * m
        bcols = np.ceil(nx / np.float(m))
        brows = np.ceil(ny / np.float(m))
        r0 = (nx / np.float(m) - 1) / 2
        s0 = (ny / np.float(m) - 1) / 2

    return (map_equatorial_radius_m, map_eccentricity, \
            e2, map_reference_latitude, map_reference_longitude, \
            map_second_reference_latitude, sin_phi1, cos_phi1, kz, \
            map_scale, bcols, brows, r0, s0, epsilon)


def easegrid(iopt, alat, alon, ascale):
    """EASE grid transformation
    (thelon thelat)=easegrid(iopt,lat,lon,ascale)

    computes the forward "ease" grid transform

    given a lat,lon (alat,alon) and the scale (ascale) the image
    transformation coordinates (thelon,thelat) are comuted
    using the "ease grid" (version 1.0) transformation given in fortran
    source code supplied by nsidc.

    the radius of the earth used in this projection is imbedded into
    ascale while the pixel dimension in km is imbedded in bscale
    the base values are: radius earth= 6371.228 km
                 pixel dimen =25.067525 km
    then, bscale = base_pixel_dimen
          ascale = radius_earth/base_pixel_dimen

    iopt is ease type: iopt=11=north, iopt=12=south, iopt=13=cylindrical
    """
    # ported from easegrid.m by JPB 21 Sept 2011
    pi2 = np.pi / 2.0
    dtr = pi2 / 90.0

    if iopt == 11:  # ease grid north
        thelon = ascale * sin(alon * dtr) * sin(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * cos(alon * dtr) * sin(dtr * (45.0 - 0.5 * alat))
    elif iopt == 12:  # ease grid south
        thelon = ascale * sin(alon * dtr) * cos(dtr * (45.0 - 0.5 * alat))
        thelat = ascale * cos(alon * dtr) * cos(dtr * (45.0 - 0.5 * alat))
    elif iopt == 13:  # ease cylindrical
        thelon = ascale * pi2 * alon * cos(30.0 * dtr) / 90.0
        thelat = ascale * sin(alat * dtr) / cos(30.0 * dtr)

    return thelon, thelat


def ease2grid(iopt, alat, alon, ascale, bscale):
    """EASE2 grid transformation

    (thelon thelat)=ease2grid(iopt,lat,lon,ascale,bscale)

	given a lat,lon (alat,alon) and the scale (ascale) the image
        transformation coordinates (thelon,thelat) are comuted
	using the "ease2 grid" (version 2.0) transformation given in IDL
	source code supplied by MJ Brodzik

	RADIUS EARTH=6378.137 KM (WGS 84)
	MAP ECCENTRICITY=0.081819190843 (WGS84)

	inputs:
	  iopt: projection type 8=EASE2 N, 9-EASE2 S, 10=EASE2 T/M
	  alon, alat: lon, lat (deg) to convert (can be outside of image)
          ascale and bscale should be integer valued)
	  ascale: grid scale factor (0..5)  pixel size is (bscale/2^ascale)
	  bscale: base grid scale index (ind=int(bscale))

          see ease2helper.py for definitions of isc and ind

	outputs:
	  thelon: X coordinate in pixels (can be outside of image)
	  thelat: Y coordinate in pixels (can be outside of image)
    """

    DTR = 0.01745329241994
    ind = round(bscale)
    isc = round(ascale)
    dlon = alon
    phi = DTR * alat
    lam = dlon

    # get base EASE2 map projection parameters
    (map_equatorial_radius_m, map_eccentricity, \
     e2, map_reference_latitude, map_reference_longitude, \
     map_second_reference_latitude, sin_phi1, cos_phi1, kz, \
     map_scale, bcols, brows, r0, s0, epsilon) = ease2_map_info(iopt, isc, ind)

    dlon = dlon - map_reference_longitude
    dlon = easeconv_normalize_degrees(dlon)
    lam = DTR * dlon

    sin_phi = np.sin(phi)

    q = (1.0 - e2) * ((sin_phi / (1.0 - e2 * sin_phi * sin_phi)) \
                      - (1.0 / (2.0 * map_eccentricity)) \
                      * np.log((1.0 - map_eccentricity * sin_phi) \
                               / (1.0 + map_eccentricity * sin_phi)))

    if iopt == 8:  # EASE2 grid north
        qp = 1.0 - ((1.0 - e2) / (2.0 * map_eccentricity) \
                    * np.log((1.0 - map_eccentricity) \
                             / (1.0 + map_eccentricity)))
        rho = map_equatorial_radius_m * np.sqrt(qp - q)
        if np.size(rho) > 1:
            rho[np.abs(qp - q) < epsilon] = 0.0
        else:
            if np.abs(qp - q) < epsilon:
                rho = 0

        x = rho * np.sin(lam)
        y = -rho * np.cos(lam)

    elif iopt == 9:  # EASE2 grid south
        qp = 1.0 - ((1.0 - e2) / (2.0 * map_eccentricity) \
                    * np.log((1.0 - map_eccentricity) \
                             / (1.0 + map_eccentricity)))
        rho = map_equatorial_radius_m * np.sqrt(qp + q)
        if np.size(rho) > 1:
            rho[np.abs(qp - q) < epsilon] = 0.0
        else:
            if np.abs(qp - q) < epsilon:
                rho = 0

        x = rho * np.sin(lam)
        y = rho * np.cos(lam)

    elif iopt == 10:  # EASE2 cylindrical
        x = map_equatorial_radius_m * kz * lam
        y = (map_equatorial_radius_m * q) / (2.0 * kz)

    else:
        print('*** invalid EASE2 projection specificaion in ease2grid')

    idx_col = round(r0 + (x / map_scale))
    idx_row = round(s0 - (y / map_scale))

    #Point cannot exceed map boundaries
    if iopt == 10 and bscale == 2:
        rows = 406
        cols = 964
    elif iopt == 10 and bscale == 5:
        rows = 4872
        cols = 11568
    elif iopt == 10 and bscale == 4:
        rows = 1624
        cols = 3856
    else:
        print('rows and cols not found')
    col = max(0, min(cols, idx_col))
    row = max(0, min(rows, idx_row))
    return (int(row),int(col))

# ==========================================================
# SMAP L2 (3 km)
# ==========================================================
if res == 3:
    print("🛰️ Running SMAP L2 (3 km) Processing Pipeline")
    
    iopt, ascale, bscale = 10, 0, 5
    n_row, n_col = 4872, 11568

    # fdir = Path("/dodrio/scratch/projects/2022_200/project_input/rsda/l_data/obs_satellite/SMAP/SPL2SMAP_E.003/SMAP_downloads/")
    out_dir = output_root / f"SMAP_{res}km_processed"
    out_dir.mkdir(exist_ok=True, parents=True)

    # ----------------------------------------------------------
    # Helper: Parse scene center
    # ----------------------------------------------------------
    def parse_scene_center(fname):
        pattern = re.compile(r'(\d+)([EW])(\d+)([NS])')
        match = pattern.search(fname)
        if match:
            lon_val, lon_dir, lat_val, lat_dir = match.groups()
            lon = int(lon_val) * (1 if lon_dir == 'E' else -1)
            lat = int(lat_val) * (1 if lat_dir == 'N' else -1)
            return lon, lat
        return None, None

    start = ease2grid(iopt, lat_start, lon_start, ascale, bscale)
    end = ease2grid(iopt, lat_end, lon_end, ascale, bscale)
    rID_start, cID_start = start
    rID_end, cID_end = end
    d_rows = np.arange(rID_start, rID_end + 1)
    d_cols = np.arange(cID_start, cID_end + 1)
    lat_subset = easegrid_latlon_3km.latitude.values[np.ix_(d_rows, d_cols)]
    lon_subset = easegrid_latlon_3km.longitude.values[np.ix_(d_rows, d_cols)]

    files = sorted(fdir.glob(f"{year}*/*.h5"))
    filtered_files = []
    lat_min, lat_max = lat_end, lat_start
    lon_min, lon_max = lon_start, lon_end

    for f in files:
        lon, lat = parse_scene_center(f.name)
        if lon is None or lat is None:
            continue
        if lon_min <= lon <= lon_max and lat_min <= lat <= lat_max:
            filtered_files.append(f)

    print(f"Selected {len(filtered_files)} granules out of {len(files)} total")

    # ----------------------------------------------------------
    # Read single granule
    # ----------------------------------------------------------
    def read_smap_granule(f):
        with h5py.File(f, "r") as tmp:
            row = tmp["Soil_Moisture_Retrieval_Data_3km"]["EASE_row_index_3km"][:]
            col = tmp["Soil_Moisture_Retrieval_Data_3km"]["EASE_column_index_3km"][:]
            sm = tmp["Soil_Moisture_Retrieval_Data_3km"]["soil_moisture_3km"][:]
            qf = tmp["Soil_Moisture_Retrieval_Data_3km"]["retrieval_qual_flag_3km"][:]

        valid_mask = ((qf & 1) == 0) & ~np.isnan(sm)
        row_valid, col_valid, sm_valid = row[valid_mask], col[valid_mask], sm[valid_mask].astype(np.float32)

        mask = (row_valid >= rID_start) & (row_valid <= rID_end) & \
               (col_valid >= cID_start) & (col_valid <= cID_end)
        row_subset = row_valid[mask] - rID_start
        col_subset = col_valid[mask] - cID_start
        sm_subset = sm_valid[mask]

        if len(sm_subset) == 0:
            return None

        ssm_step = np.full((len(d_rows), len(d_cols)), np.nan, dtype=np.float32)
        ssm_step[row_subset, col_subset] = sm_subset
        return ssm_step

    # ----------------------------------------------------------
    # Merge and write
    # ----------------------------------------------------------
    tiles_by_time = defaultdict(list)
    for f in filtered_files:
        smap_time = pd.to_datetime(f.name.split("_")[5], format="%Y%m%dT%H%M%S")
        tiles_by_time[smap_time].append(f)

    month_groups = defaultdict(list)
    for t, files_for_time in tiles_by_time.items():
        month_groups[t.month].append((t, files_for_time))

    def merge_tiles(files_for_time):
        ssm_accum = np.zeros((len(d_rows), len(d_cols)), dtype=np.float32)
        ssm_count = np.zeros((len(d_rows), len(d_cols)), dtype=np.int32)
        for f in files_for_time:
            ssm_step = read_smap_granule(f)
            if ssm_step is None:
                continue
            valid_mask = ~np.isnan(ssm_step)
            ssm_accum[valid_mask] += ssm_step[valid_mask]
            ssm_count[valid_mask] += 1
        with np.errstate(invalid="ignore", divide="ignore"):
            ssm_merged = np.where(ssm_count > 0, ssm_accum / ssm_count, np.nan)
        return ssm_merged

    for month, time_file_list in sorted(month_groups.items()):
        print(f"\n📦 Processing {year}-{month:02d}, {len(time_file_list)} SMAP times")
        results = []
        with ProcessPoolExecutor(max_workers=12) as executor:
            futures = {executor.submit(merge_tiles, files_for_time): t for t, files_for_time in time_file_list}
            for fut in tqdm(as_completed(futures), total=len(futures), desc=f"Month {month:02d}"):
                t = futures[fut]
                ssm_merged = fut.result()
                if np.all(np.isnan(ssm_merged)) or np.all(ssm_merged == 0):
                    continue
                results.append((t, ssm_merged))
        if not results:
            print(f"⚠️ No valid data for month {month}, skipping...")
            continue

        results.sort(key=lambda x: x[0])
        time_list, ssm_list = zip(*results)
        ssm_all = np.stack(ssm_list, axis=0)
        time_index = pd.DatetimeIndex(time_list)

        nc_file = out_dir / f"SMAP_{res}km_processed_{region}_{year}_{month:02d}.nc"
        ds = xr.Dataset(
            {"ssm": (("time", "y", "x"), ssm_all)},
            coords={"time": time_index, "lat": (("y", "x"), lat_subset), "lon": (("y", "x"), lon_subset)},
        )
        ds.to_netcdf(nc_file, format="NETCDF4", encoding={"ssm": {"zlib": True, "complevel": 4}}, unlimited_dims=["time"])
        print(f"✅ Saved monthly NetCDF: {nc_file}")
        del ssm_all, ds
        gc.collect()

    print("\n✅ Completed all months for 3 km processing.")

# ==========================================================
# SMAP L3 (9 km)
# ==========================================================
elif res == 9:
    print("🛰️ Running SMAP L3 (9 km) Processing Pipeline with parallel merging")
    
    iopt, ascale, bscale = 10, 0, 4
    n_row, n_col = 1624, 3856
    
    start = ease2grid(iopt, lat_start, lon_start, ascale, bscale)
    end = ease2grid(iopt, lat_end, lon_end, ascale, bscale)
    rID_start, cID_start = start
    rID_end, cID_end = end
    d_rows = np.arange(rID_start, rID_end + 1)
    d_cols = np.arange(cID_start, cID_end + 1)
    lat_subset = easegrid_latlon_9km.latitude.values[np.ix_(d_rows, d_cols)]
    lon_subset = easegrid_latlon_9km.longitude.values[np.ix_(d_rows, d_cols)]

    out_dir = output_root / f"SMAP_{res}km_processed"
    out_dir.mkdir(exist_ok=True, parents=True)
    output_nc = out_dir / f"SMAP_{res}km_processed_{region}_{year}.nc"

    # ----------------------------------------------------------
    # Helper to parse datetime from filename
    # ----------------------------------------------------------
    def parse_smap_datetime(filename):
        match = re.search(r"SMAP_L2_SM_P_E_\d+_([AD])_(\d{8}T\d{6})_", filename)
        if not match:
            return None, None
        pass_type = match.group(1)
        timestamp = datetime.strptime(match.group(2), "%Y%m%dT%H%M%S")
        return timestamp, pass_type

    # ----------------------------------------------------------
    # Read single granule and subset to region
    # ----------------------------------------------------------
    def read_smap_granule_9km(fp):
        smap_time, pass_type = parse_smap_datetime(fp.name)
        if smap_time is None:
            return None

        nominal_time = smap_time.replace(hour=6 if pass_type == "D" else 18, minute=0, second=0)
        
        with h5py.File(fp, "r") as f:
            row = f['Soil_Moisture_Retrieval_Data']['EASE_row_index'][:]
            col = f['Soil_Moisture_Retrieval_Data']['EASE_column_index'][:]
            sm = f["Soil_Moisture_Retrieval_Data"]["soil_moisture"][:]
            flag = f["Soil_Moisture_Retrieval_Data"]["retrieval_qual_flag"][:]
            valid = ((flag & 1) == 0) & (~np.isnan(sm))

        row_valid, col_valid, sm_valid = row[valid], col[valid], sm[valid].astype(np.float32)
        mask = (row_valid >= rID_start) & (row_valid <= rID_end) & (col_valid >= cID_start) & (col_valid <= cID_end)
        row_subset = row_valid[mask] - rID_start
        col_subset = col_valid[mask] - cID_start
        sm_subset  = sm_valid[mask]
        if len(sm_subset) == 0:
            return None

        sm_grid = np.full((len(d_rows), len(d_cols)), np.nan, dtype=np.float32)
        sm_grid[row_subset, col_subset] = sm_subset
        return nominal_time, sm_grid

    # ----------------------------------------------------------
    # List all files
    # ----------------------------------------------------------
    files = sorted(fdir.glob(f'{year}*/*.h5'))
    print(f"Found {len(files)} SMAP files for {year}")

    # ----------------------------------------------------------
    # Parallel read & merge
    # ----------------------------------------------------------
    results = []
    with ProcessPoolExecutor(max_workers=12) as executor:
        futures = {executor.submit(read_smap_granule_9km, fp): fp for fp in files}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Reading SMAP granules"):
            res_item = fut.result()
            if res_item is None:
                continue
            results.append(res_item)

    if not results:
        print("⚠️ No valid data found!")
        exit()

    # Sort by time
    results.sort(key=lambda x: x[0])
    time_list, sm_list = zip(*results)

    # Stack into array
    soil_moisture_arr = np.stack(sm_list, axis=0)
    time_arr = np.array(time_list)

    # Create xarray Dataset
    ds = xr.Dataset(
        {"ssm": (("time", "y", "x"), soil_moisture_arr)},
        coords={
            "time": time_arr,
            "lat": (("y", "x"), lat_subset),
            "lon": (("y", "x"), lon_subset),
        },
    )

    ds["ssm"].attrs.update({
        "long_name": "Surface Soil Moisture (masked where retrieval_qual_flag != 0)",
        "units": "cm3/cm3",
        "source": "SMAP L3 Enhanced 9 km",
    })
    ds.attrs.update({
        "title": "SMAP L3 Enhanced Soil Moisture 9 km (AM/PM, masked for valid retrievals)",
        "history": f"Created on {datetime.now().isoformat()}",
    })

    ds.to_netcdf(output_nc, format="NETCDF4", encoding={"ssm": {"zlib": True, "complevel": 4, "dtype": "float32"}})
    print(f"✅ Saved masked dataset to {output_nc}")

print("\n All processing complete.")

