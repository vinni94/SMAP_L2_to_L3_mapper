# 🛰️ SMAP L2 to L3 Soil Moisture Mapper

This repository provides a **Python-based processing pipeline** for transforming NASA SMAP (Soil Moisture Active Passive) **Level 2 (3 km)** and **Level 3 (9 km)** products into regionally subsetted, time-sorted, and quality-masked **NetCDF datasets**.

The script supports large-scale processing with **parallel computation** using `ProcessPoolExecutor`, automatically filters granules by region and date, and outputs monthly or yearly processed NetCDF files ready for analysis.

---

## 🚀 Features

- ✅ Supports both **L2 (3 km)** and **L3 (9 km)** SMAP soil moisture products  
- ✅ Region-based spatial subsetting using **EASE-Grid 2.0 projection**  
- ✅ Quality flag masking to exclude unreliable retrievals  
- ✅ **Parallelized tile merging** for efficient computation  
- ✅ Outputs compressed NetCDF4 files (`.nc`)  
- ✅ Easily configurable via command-line arguments  

---

## 📦 Installation

Clone this repository and install required dependencies:

```bash
git clone https://github.com/vinni94/SMAP_L2_to_L3_mapper.git
cd SMAP_L2_to_L3_mapper
pip install -r requirements.txt
```

## ⚙️ Dependencies

The script uses the following Python libraries:

```text
numpy
xarray
h5py
pandas
tqdm
argparse
```
## 🚀 Usage

This script processes SMAP L2 (3 km) and L3 (9 km) soil moisture data for a specified region and year.  

### 1. Prepare the Run Script

Update the `run_smap_pipeline.sh` file with:

- `INPUT_DIR`: Path to the folder containing SMAP HDF5 files.  
- `OUTPUT_ROOT`: Path where processed NetCDFs will be saved.  
- `PYTHON_SCRIPT`: Path to the main Python processing script (`process_smap_l2_to_l3.py`).  
- `LAT_START`, `LAT_END`, `LON_START`, `LON_END`: Latitude/longitude bounds of your region.  
- `REGION`: Shortname for your region (e.g., `EU`).  
- `RESOLUTIONS`: List of resolutions to process (`3` for L2, `9` for L3).  
- `YEARS`: List of years to process.  
- `EASEGRID_9KM` / `EASEGRID_3KM`: Paths to global EASEgrid latitude/longitude files.  

### 2. Run the Pipeline

Make the script executable and run it:

```bash
chmod +x run_smap_pipeline.sh
./run_smap_pipeline.sh
```

## 🗂️ Output Structure

For each region and year, the script generates processed NetCDF files:

output_root/
└── SMAP_9km_processed/
├── SMAP_9km_processed_<REGION><YEAR>.nc # 9 km dataset
└── SMAP_3km_processed/
├── SMAP_3km_processed<REGION><YEAR>01.nc # 3 km monthly datasets
├── SMAP_3km_processed<REGION><YEAR>_02.nc
└── ...


Each NetCDF file contains the following variables:

- `ssm`: Surface soil moisture (masked where retrieval flags indicate poor quality)  
- `lat`: Latitude coordinates of the subsetted region  
- `lon`: Longitude coordinates of the subsetted region  
- `time`: Observation timestamps corresponding to the measurements


## 🪪 Citation

If you use this script for research or publication, please cite:  

> NASA SMAP Level 2 & 3 Soil Moisture Products  
> Jet Propulsion Laboratory (JPL) / National Snow and Ice Data Center (NSIDC)  

Reference the official SMAP documentation for details on the datasets:  
[https://nsidc.org/data/smap](https://nsidc.org/data/smap)

