# PREPROCESSING SCRIPTS README

This folder contains preprocessing scripts for refining key datasets used in the burn severity analysis. The scripts clean and subset data related to wildfire perimeters, burn severity records, forest disturbances, and administrative boundaries.

## Contents
1) 0_export_admin_boundary_shps.py — Extracts and exports Malheur National Forest boundary for spatial analysis. **[Must be run before other scripts in folder]**
2) preprocess_ads_data.py — Processes Aerial Detection Survey (ADS) outbreak layers.
3) preprocess_mnf_fire_perimeters.py — Subsets historical fire perimeters and extracts relevant burn severity records.

## Libraries required
* geopandas
* pandas
* fiona
* zipfile
* contextily
* matplotlib
* os