#!/bin/bash

# ---------------------------
# CONFIGURATION
# ---------------------------
INPUT_DIR="/staging/leuven/stg_00024/l_data/obs_satellite/SMAP/SPL2SMP_E.006"
OUTPUT_ROOT="/scratch/leuven/378/vsc37824/SMAP_Processed"
PYTHON_SCRIPT="process_smap_l2_to_l3.py"

# ---------------------------
# Region extent (EU)
# ---------------------------
LAT_START=60   # northernmost latitude
LAT_END=35     # southernmost latitude
LON_START=-10  # westernmost longitude
LON_END=40     # easternmost longitude
REGION="EU"

# ---------------------------
# Resolutions and years to process
# ---------------------------
RESOLUTIONS=(9)  # Can add 3 for L2 processing
YEARS=(2016)

# ---------------------------
# Global EASEGRID lat/lon files
# ---------------------------
EASEGRID_9KM="/scratch/leuven/378/vsc37824/SMAP_3km_processed/NSIDC0772_LatLon_EASE2_M09km_v1.1.nc"
EASEGRID_3KM="/scratch/leuven/378/vsc37824/SMAP_3km_processed/NSIDC0772_LatLon_EASE2_M03km_v1.1.nc"

# ---------------------------
# Loop over resolutions and years
# ---------------------------
for RES in "${RESOLUTIONS[@]}"; do
    for YEAR in "${YEARS[@]}"; do
        echo "Running SMAP processing: RES=${RES}, YEAR=${YEAR}, REGION=${REGION}"

        python "$PYTHON_SCRIPT" \
            --res "$RES" \
            --year "$YEAR" \
            --lat_start "$LAT_START" \
            --lat_end "$LAT_END" \
            --lon_start "$LON_START" \
            --lon_end "$LON_END" \
            --input_dir "$INPUT_DIR" \
            --output_root "$OUTPUT_ROOT" \
            --region "$REGION" \
            --easegrid_9km_file "$EASEGRID_9KM" \
            --easegrid_3km_file "$EASEGRID_3KM"

        echo "Completed processing for RES=${RES}, YEAR=${YEAR}, REGION=${REGION}"
    done
done

echo "All processing complete!"

