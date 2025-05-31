# SCRIPTS

## Contents
(1) The "preprocessing_scripts" folder which contains that were run locally due to file storage contraints and whose outputs are uploaded to the "subsetted_data" folder (CHAR/subsetted_data) 
(2) A set of four scripts comprising the main workflow of the outbreak and burn severity analysis. This folder contains scripts that process wildfire burn severity data, correlate fire severity with elevation, and analyze disturbance interactions. The workflow follows a sequential order, relying on outputs from the preprocessing scripts.

## Execution Order
1) 01_extract_outbreak_polygons.py
* Uses preprocessed ADS data to extract WSB outbreak polygons
* Filters years 1965–1990 for analysis
* Dissolves outbreak boundaries for yearly aggregation
* **Outputs**: wsb_outbreaks_1965_1990.geojson, time series of WSB outbreak progression plots
* **Dependencies**: Relies on output of preprocessing script "preprocess_mnf_fire_perimeters.ipynb"

2) 02_attribute_mnf_fires.py
* Cleans MNF fire perimeter dataset
* Joins Monitoring Trends in Burn Severity (MTBS) records with fire perimeters
* Removes duplicate incidents using overlap analysis
* Filters fires within relevant years (1984–2002) for burn severity study
* **Outputs**: mnf_fire_perimeters_cleaned.geojson, target_fires.geojson

3) 03_calculate_NBR.py
* Searches pre- and post-fire Landsat 5 imagery
* Computes NBR and dNBR using Near Infrared (NIR) and Shortwave Infrared (SWIR2) bands
* Applies Key & Benson (2006) burn severity classification thresholds
* Generates multi-year post-fire recovery analysis (1, 2, 5 years post-fire)
* Compares results to MTBS classifications
* **Outputs**: fire_analysis_results.geojson, maps & plots showing NBR/dNBR changes over time, statistical summaries for burn severity categories

4) 04_elevation_analysis.py
* Correlates fire severity (dNBR) with elevation
* Extracts elevation data using OpenTopoData API
* Performs statistical correlation analysis
* Visualizes burn severity distributions by elevation bins
* Generates six-panel plots summarizing findings
* **Outputs**: fire_elevation_analysis.geojson (fire severity v. elevation results), scatter plots, stacked burn severity histograms, and correlation summaries

**NOTE: Execution order is critical. Running notebooks out of order will yield invalid results.**

## Libraries required
* geopandas
* pandas
* rasterio
* contextily
* numpy
* matplotlib
* pystac_client
* planetary_computer
* rioxarray
* ipywidgets
* seaborn
* scipy
* requests
