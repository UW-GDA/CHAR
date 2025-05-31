"""
COMPLETE WORKFLOW DOCUMENTATION:

Fire Analysis Elevation and Burn Severity Module

Analyzes wildfire burn severity using Landsat 5 NBR and correlates with elevation data.
Clean, readable code for graduate students without unnecessary complexity.

Functions:
- clean_year(): Convert messy year data to clean integers
- find_fire(): Find a fire by name in the dataset
- search_landsat5_imagery(): Get Landsat 5 images around a target date
- calculate_nbr(): Calculate NBR from Landsat bands
- get_burn_severity_stats(): Calculate burn severity statistics
- get_elevation_points(): Get elevation data using OpenTopoData API
- analyze_elevation_vs_severity(): Correlate fire severity with elevation
- plot_results(): Make comprehensive 6-panel plots
- analyze_fire(): Main analysis function with customizable post-fire delay
- run_fire_elevation_analysis_by_date(): Complete workflow function

Usage:
    results = run_fire_elevation_analysis_by_date('Monument Rock', '1989-07-28')
    results = run_fire_elevation_analysis_by_date('Canyon Creek', '1989-08-04', post_fire_years=2)

The function will:
  1. Find your fire by name
  2. Get Landsat imagery before and after the fire
  3. Calculate burn severity (dNBR)
  4. Download elevation data
  5. Analyze how elevation affects fire severity
  6. Make comprehensive plots with all severity categories
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import rioxarray as rxr
import pystac_client
import planetary_computer as pc
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import requests
import time
import warnings
warnings.filterwarnings('ignore')

def clean_year(year_val):
    """Convert messy year data to clean integers"""
    if pd.isna(year_val):
        return None
    year_str = str(year_val)
    if '-' in year_str:
        return int(year_str.split('-')[0])
    return int(float(year_str))

def find_fire(fire_name, fire_file="../subsetted_data/mnf_fires_all.geojson"):
    """Find a fire by name in the dataset"""
    print(f"Looking for fire: {fire_name}")
    
    # Load fires
    fires = gpd.read_file(fire_file)
    fires['year'] = fires['FIRE_YEAR'].apply(clean_year)
    fires = fires.dropna(subset=['year']).to_crs('EPSG:4326')
    
    # Search
    matches = fires[fires['INCIDENT'].str.contains(fire_name, case=False, na=False)]
    
    if len(matches) == 0:
        print(f"No fires found with '{fire_name}' in the name")
        return None
    
    if len(matches) > 1:
        print(f"Found {len(matches)} matches:")
        for i, fire in matches.iterrows():
            print(f"  {fire['INCIDENT']} ({int(fire['year'])})")
    
    fire = matches.iloc[0]
    print(f"Using: {fire['INCIDENT']} ({int(fire['year'])})")
    return fire

def search_landsat5_imagery(geometry, target_date, days_buffer=30):
    """
    Search for Landsat 5 imagery around a specific target date
    
    Parameters:
    -----------
    geometry : shapely geometry
        Fire boundary polygon
    target_date : datetime
        The target date to search around
    days_buffer : int
        Number of days to search before/after the target date (default: 30)
    """
    
    start_date = (target_date - timedelta(days=days_buffer)).strftime('%Y-%m-%d')
    end_date = (target_date + timedelta(days=days_buffer)).strftime('%Y-%m-%d')
    
    print(f"Searching Landsat 5 imagery from {start_date} to {end_date}")
    print(f"  Target date: {target_date.strftime('%Y-%m-%d')}")
    
    # Initialize STAC catalog
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1/",
        modifier=pc.sign_inplace
    )
    
    # Search for Landsat 5 data
    search = catalog.search(
        collections=["landsat-c2-l2"],
        intersects=geometry.__geo_interface__,
        datetime=f"{start_date}/{end_date}",
        query={
            "platform": {"in": ["landsat-5"]},
            "eo:cloud_cover": {"lt": 50}  # Relaxed cloud cover for 1980s
        }
    )
    
    items = list(search.items())
    
    if len(items) == 0:
        print("No images found with cloud filter, trying without cloud filter...")
        search = catalog.search(
            collections=["landsat-c2-l2"],
            intersects=geometry.__geo_interface__,
            datetime=f"{start_date}/{end_date}",
            query={"platform": {"in": ["landsat-5"]}}
        )
        items = list(search.items())
    
    if len(items) == 0:
        # Try expanding search window
        expanded_days = days_buffer * 2
        start_date = (target_date - timedelta(days=expanded_days)).strftime('%Y-%m-%d')
        end_date = (target_date + timedelta(days=expanded_days)).strftime('%Y-%m-%d')
        
        print(f"Expanding search window to {start_date} to {end_date}")
        
        search = catalog.search(
            collections=["landsat-c2-l2"],
            intersects=geometry.__geo_interface__,
            datetime=f"{start_date}/{end_date}",
            query={"platform": {"in": ["landsat-5"]}}
        )
        items = list(search.items())
    
    if len(items) == 0:
        print(f"No Landsat 5 images found around {target_date.strftime('%Y-%m-%d')}")
        return []
    
    # Sort by date proximity to target, then cloud cover
    items_with_info = []
    target_date_clean = target_date
    if hasattr(target_date, 'tz') and target_date.tz is not None:
        target_date_clean = target_date.tz_localize(None)
    
    for item in items:
        cloud_cover = item.properties.get('eo:cloud_cover', 100)
        item_date = pd.to_datetime(item.properties['datetime'])
        if hasattr(item_date, 'tz') and item_date.tz is not None:
            item_date = item_date.tz_localize(None)
        
        date_diff = abs((item_date - target_date_clean).days)
        date_str = item_date.strftime('%Y-%m-%d')
        items_with_info.append((item, cloud_cover, date_diff, date_str))
    
    # Sort by date proximity first (closer to target = better), then cloud cover
    items_with_info.sort(key=lambda x: (x[2], x[1]))
    sorted_items = [item for item, cloud, date_diff, date_str in items_with_info]
    
    print(f"Found {len(sorted_items)} images:")
    for i, (item, cloud, date_diff, date_str) in enumerate(items_with_info[:3]):
        print(f"  Image {i+1}: {date_str}, {cloud:.1f}% clouds, {date_diff} days from target")
    
    return sorted_items

def calculate_nbr(landsat_item, fire_polygon):
    """Calculate NBR from Landsat bands"""
    signed_item = pc.sign(landsat_item)
    
    date = landsat_item.properties.get('datetime', 'Unknown')
    clouds = landsat_item.properties.get('eo:cloud_cover', 'Unknown')
    print(f"Processing: {date}, {clouds}% clouds")
    
    # Find the right band names (they vary)
    available_bands = list(signed_item.assets.keys())
    print(f"Available bands: {available_bands}")
    
    # Look for NIR and SWIR bands
    nir_band = None
    swir_band = None
    
    for band in available_bands:
        if band in ['B04', 'nir08']:
            nir_band = band
        elif band in ['B07', 'swir22']:
            swir_band = band
    
    if nir_band is None or swir_band is None:
        print(f"ERROR: Can't find NIR ({nir_band}) or SWIR ({swir_band}) bands")
        return None
    
    print(f"Using {nir_band} (NIR) and {swir_band} (SWIR)")
    
    # Load NIR and SWIR bands
    nir = rxr.open_rasterio(signed_item.assets[nir_band].href)
    swir = rxr.open_rasterio(signed_item.assets[swir_band].href)
    
    # Convert fire to raster CRS and buffer it
    fire_gdf = gpd.GeoDataFrame([1], geometry=[fire_polygon], crs='EPSG:4326')
    fire_proj = fire_gdf.to_crs(nir.rio.crs).geometry.iloc[0]
    fire_buffered = fire_proj.buffer(2000)  # 2km buffer
    
    # Crop to fire area
    nir_crop = nir.rio.clip([fire_buffered], drop=True)
    swir_crop = swir.rio.clip([fire_buffered], drop=True)
    
    # Scale to reflectance (Landsat Collection 2)
    nir_refl = nir_crop * 0.0000275 - 0.2
    swir_refl = swir_crop * 0.0000275 - 0.2
    
    # Calculate NBR
    nbr = (nir_refl - swir_refl) / (nir_refl + swir_refl)
    
    # Clip to exact fire boundary
    nbr_fire = nbr.rio.clip([fire_proj], drop=True)
    
    print(f"NBR calculated, shape: {nbr_fire.shape}")
    return nbr_fire

def get_burn_severity_stats(dnbr_values):
    """Calculate burn severity statistics"""
    high = np.sum(dnbr_values > 0.66)
    mod_high = np.sum((dnbr_values >= 0.44) & (dnbr_values <= 0.66))
    mod_low = np.sum((dnbr_values >= 0.25) & (dnbr_values < 0.44))
    low = np.sum((dnbr_values >= 0.1) & (dnbr_values < 0.25))
    unburned = np.sum(dnbr_values < 0.1)
    total = len(dnbr_values)
    
    return {
        'high': (high, 100 * high / total),
        'mod_high': (mod_high, 100 * mod_high / total),
        'mod_low': (mod_low, 100 * mod_low / total),
        'low': (low, 100 * low / total),
        'unburned': (unburned, 100 * unburned / total),
        'total': total
    }

def get_elevation_points(fire_polygon, n_points=400):
    """Get elevation data using OpenTopoData API"""
    print("Getting elevation data...")
    
    # Get fire bounds with buffer
    fire_gdf = gpd.GeoDataFrame([1], geometry=[fire_polygon], crs='EPSG:4326')
    fire_utm = fire_gdf.to_crs('EPSG:3857')
    buffered = fire_utm.buffer(3000).to_crs('EPSG:4326')  # 3km buffer
    bounds = buffered.bounds.iloc[0]
    
    # Create grid of points
    n_side = int(np.sqrt(n_points))
    lats = np.linspace(bounds[1], bounds[3], n_side)
    lons = np.linspace(bounds[0], bounds[2], n_side)
    
    elevations = []
    coords = []
    
    for lat in lats:
        for lon in lons:
            try:
                url = f"https://api.opentopodata.org/v1/srtm30m?locations={lat},{lon}"
                response = requests.get(url, timeout=5)
                
                if response.status_code == 200:
                    data = response.json()
                    if data['status'] == 'OK' and len(data['results']) > 0:
                        elev = data['results'][0]['elevation']
                        if elev is not None:
                            elevations.append(elev)
                            coords.append((lon, lat))
                
                time.sleep(0.05)  # Be nice to the API
                
            except:
                continue
    
    print(f"Got {len(elevations)} elevation points")
    
    if len(elevations) < 50:
        print("Not enough elevation data")
        return None
    
    # Simple interpolation to create a grid
    from scipy.interpolate import griddata
    
    target_lats = np.linspace(bounds[1], bounds[3], 100)
    target_lons = np.linspace(bounds[0], bounds[2], 100)
    lon_grid, lat_grid = np.meshgrid(target_lons, target_lats)
    
    elev_grid = griddata(coords, elevations, (lon_grid, lat_grid), method='cubic')
    
    # Fill missing values
    if np.isnan(elev_grid).any():
        elev_grid_nn = griddata(coords, elevations, (lon_grid, lat_grid), method='nearest')
        elev_grid = np.where(np.isnan(elev_grid), elev_grid_nn, elev_grid)
    
    # Convert to xarray
    import xarray as xr
    elevation = xr.DataArray(
        elev_grid,
        coords={'y': target_lats, 'x': target_lons},
        dims=['y', 'x']
    ).rio.write_crs('EPSG:4326')
    
    print(f"Elevation range: {float(elevation.min()):.0f} - {float(elevation.max()):.0f} m")
    return elevation

def analyze_elevation_vs_severity(dnbr, elevation, fire_polygon):
    """Correlate fire severity with elevation"""
    print("Analyzing elevation vs severity...")
    
    # Make sure both datasets match
    if dnbr.rio.crs != elevation.rio.crs:
        elevation = elevation.rio.reproject_match(dnbr)
    
    # Clip both to fire boundary
    fire_gdf = gpd.GeoDataFrame([1], geometry=[fire_polygon], crs='EPSG:4326')
    fire_proj = fire_gdf.to_crs(dnbr.rio.crs)
    
    dnbr_fire = dnbr.rio.clip(fire_proj.geometry, drop=True)
    elev_fire = elevation.rio.clip(fire_proj.geometry, drop=True)
    
    # Match shapes
    if dnbr_fire.shape != elev_fire.shape:
        elev_fire = elev_fire.rio.reproject_match(dnbr_fire)
    
    # Get clean data
    dnbr_flat = dnbr_fire.values.flatten()
    elev_flat = elev_fire.values.flatten()
    
    # Remove NaN values
    valid = ~(np.isnan(dnbr_flat) | np.isnan(elev_flat))
    dnbr_clean = dnbr_flat[valid]
    elev_clean = elev_flat[valid]
    
    print(f"Analyzing {len(dnbr_clean)} pixels")
    
    # Calculate correlation
    from scipy.stats import pearsonr
    correlation, p_value = pearsonr(dnbr_clean, elev_clean)
    
    # Bin by elevation
    n_bins = 6
    elev_bins = np.linspace(elev_clean.min(), elev_clean.max(), n_bins + 1)
    bin_centers = (elev_bins[:-1] + elev_bins[1:]) / 2
    
    severity_by_elev = []
    for i in range(n_bins):
        in_bin = (elev_clean >= elev_bins[i]) & (elev_clean < elev_bins[i + 1])
        if i == n_bins - 1:  # Include max value in last bin
            in_bin = (elev_clean >= elev_bins[i]) & (elev_clean <= elev_bins[i + 1])
        
        if np.sum(in_bin) > 0:
            bin_dnbr = dnbr_clean[in_bin]
            
            # Calculate all severity categories
            high_sev = 100 * np.sum(bin_dnbr > 0.66) / len(bin_dnbr)
            mod_high_sev = 100 * np.sum((bin_dnbr >= 0.44) & (bin_dnbr <= 0.66)) / len(bin_dnbr)
            mod_low_sev = 100 * np.sum((bin_dnbr >= 0.25) & (bin_dnbr < 0.44)) / len(bin_dnbr)
            low_sev = 100 * np.sum((bin_dnbr >= 0.1) & (bin_dnbr < 0.25)) / len(bin_dnbr)
            unburned = 100 * np.sum(bin_dnbr < 0.1) / len(bin_dnbr)
            
            severity_by_elev.append({
                'elevation': bin_centers[i],
                'mean_dnbr': np.mean(bin_dnbr),
                'high_severity': high_sev,
                'moderate_high': mod_high_sev,
                'moderate_low': mod_low_sev,
                'low_severity': low_sev,
                'unburned': unburned,
                'pixels': len(bin_dnbr)
            })
    
    return {
        'dnbr': dnbr_clean,
        'elevation': elev_clean,
        'correlation': correlation,
        'p_value': p_value,
        'elevation_bins': pd.DataFrame(severity_by_elev)
    }

def plot_results(results, fire_name):
    """Make comprehensive plots like the original"""
    fig = plt.figure(figsize=(16, 10))
    
    # Scatter plot with trend line
    ax1 = plt.subplot(2, 3, 1)
    plt.scatter(results['elevation'], results['dnbr'], alpha=0.5, s=2, c='red')
    
    z = np.polyfit(results['elevation'], results['dnbr'], 1)
    p = np.poly1d(z)
    x_trend = np.linspace(results['elevation'].min(), results['elevation'].max(), 100)
    plt.plot(x_trend, p(x_trend), 'blue', linewidth=2, alpha=0.8)
    
    plt.xlabel('Elevation (m)')
    plt.ylabel('dNBR')
    plt.title(f'Fire Severity vs Elevation\nr = {results["correlation"]:.3f}, p = {results["p_value"]:.4f}')
    plt.grid(True, alpha=0.3)
    
    # Mean dNBR by elevation
    ax2 = plt.subplot(2, 3, 2)
    bins = results['elevation_bins']
    plt.plot(bins['elevation'], bins['mean_dnbr'], 'o-', linewidth=2, markersize=8, color='darkred')
    plt.xlabel('Elevation (m)')
    plt.ylabel('Mean dNBR')
    plt.title('Mean Fire Severity by Elevation')
    plt.grid(True, alpha=0.3)
    
    # Stacked bar chart of ALL severity categories by elevation
    ax3 = plt.subplot(2, 3, 3)
    x_pos = np.arange(len(bins))
    
    # Get all severity percentages
    high_pct = bins['high_severity']
    mod_high_pct = bins['moderate_high']
    mod_low_pct = bins['moderate_low']
    low_pct = bins['low_severity']
    unburned_pct = bins['unburned']
    
    # Create stacked bars with all categories
    plt.bar(x_pos, high_pct, label='High (>0.66)', color='darkred')
    plt.bar(x_pos, mod_high_pct, bottom=high_pct, label='Mod-High (0.44-0.66)', color='red')
    plt.bar(x_pos, mod_low_pct, bottom=high_pct + mod_high_pct, label='Mod-Low (0.25-0.44)', color='orange')
    plt.bar(x_pos, low_pct, bottom=high_pct + mod_high_pct + mod_low_pct, label='Low (0.1-0.25)', color='yellow')
    plt.bar(x_pos, unburned_pct, bottom=high_pct + mod_high_pct + mod_low_pct + low_pct, label='Unburned (<0.1)', color='green')
    
    plt.xlabel('Elevation Bin')
    plt.ylabel('Percentage')
    plt.title('Burn Severity by Elevation')
    plt.xticks(x_pos, [f"{elev:.0f}m" for elev in bins['elevation']], rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    # Elevation histogram
    ax4 = plt.subplot(2, 3, 4)
    plt.hist(results['elevation'], bins=30, alpha=0.7, color='brown', edgecolor='black')
    plt.xlabel('Elevation (m)')
    plt.ylabel('Pixels')
    plt.title('Elevation Distribution')
    plt.grid(True, alpha=0.3)
    
    # dNBR histogram
    ax5 = plt.subplot(2, 3, 5)
    plt.hist(results['dnbr'], bins=30, alpha=0.7, color='red', edgecolor='black')
    plt.xlabel('dNBR')
    plt.ylabel('Pixels')
    plt.title('dNBR Distribution')
    plt.grid(True, alpha=0.3)
    
    # Summary statistics box
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    corr_abs = abs(results['correlation'])
    strength = "Strong" if corr_abs > 0.7 else "Moderate" if corr_abs > 0.3 else "Weak"
    direction = "positive" if results['correlation'] > 0 else "negative"
    significance = "significant" if results['p_value'] < 0.05 else "not significant"
    
    summary_text = f"""ELEVATION ANALYSIS SUMMARY

Fire: {fire_name}
Pixels: {len(results['dnbr']):,}

Elevation:
• Range: {results['elevation'].min():.0f} - {results['elevation'].max():.0f} m
• Mean: {results['elevation'].mean():.0f} m

dNBR:
• Range: {results['dnbr'].min():.3f} - {results['dnbr'].max():.3f}
• Mean: {results['dnbr'].mean():.3f}

Correlation:
• r = {results['correlation']:.4f}
• p = {results['p_value']:.6f}
• {strength} {direction} correlation
• {significance} (α=0.05)"""
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.suptitle(f'{fire_name} - Fire Severity vs Elevation Analysis', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.show()

def analyze_fire(fire_name, fire_date, post_fire_years=1, fire_file="../subsetted_data/mnf_fires_all.geojson"):
    """
    Main function - analyze a fire by name and date with customizable post-fire delay
    
    Parameters:
    -----------
    fire_name : str
        Name of the fire
    fire_date : str
        Fire date in 'YYYY-MM-DD' format
    post_fire_years : float
        Years after fire to search for post-fire imagery (default: 1)
    fire_file : str
        Path to fire boundaries file
    
    Examples:
    ---------
    analyze_fire("Canyon Creek", "1989-08-04")  # 1 year delay
    analyze_fire("Canyon Creek", "1989-08-04", post_fire_years=2)  # 2 year delay
    analyze_fire("Canyon Creek", "1989-08-04", post_fire_years=0.5)  # 6 month delay
    """
    print(f"FIRE ANALYSIS: {fire_name}")
    
    # Find the fire
    fire = find_fire(fire_name, fire_file)
    if fire is None:
        return None
    
    # Set up dates
    fire_date_parsed = pd.to_datetime(fire_date)
    pre_fire_target = fire_date_parsed - timedelta(days=60)  # 2 months before fire
    post_fire_target = fire_date_parsed + timedelta(days=int(365 * post_fire_years))
    
    print(f"Fire Information:")
    print(f"  Name: {fire['INCIDENT']}")
    print(f"  Fire Date: {fire_date_parsed.strftime('%Y-%m-%d')}")
    print(f"  Pre-fire target: {pre_fire_target.strftime('%Y-%m-%d')}")
    print(f"  Post-fire target: {post_fire_target.strftime('%Y-%m-%d')} ({post_fire_years} years after)")
    
    # Get PRE-FIRE satellite images
    print(f"\n=== SEARCHING FOR PRE-FIRE IMAGERY ===")
    pre_images = search_landsat5_imagery(fire.geometry, pre_fire_target)
    if len(pre_images) == 0:
        print("No pre-fire imagery found")
        return None
    
    # Get POST-FIRE satellite images
    print(f"\n=== SEARCHING FOR POST-FIRE IMAGERY ===")
    post_images = search_landsat5_imagery(fire.geometry, post_fire_target)
    if len(post_images) == 0:
        print("No post-fire imagery found")
        return None
    
    # Calculate NBR
    print(f"\n=== CALCULATING PRE-FIRE NBR ===")
    pre_nbr = calculate_nbr(pre_images[0], fire.geometry)
    if pre_nbr is None:
        return None
    
    print(f"\n=== CALCULATING POST-FIRE NBR ===")
    post_nbr = calculate_nbr(post_images[0], fire.geometry)
    if post_nbr is None:
        return None
    
    # Calculate dNBR
    print(f"\n=== CALCULATING dNBR ===")
    if pre_nbr.rio.crs != post_nbr.rio.crs:
        post_nbr = post_nbr.rio.reproject_match(pre_nbr)
    
    dnbr = pre_nbr - post_nbr
    dnbr_values = dnbr.values.flatten()
    dnbr_clean = dnbr_values[~np.isnan(dnbr_values)]
    
    print(f"dNBR stats: mean={np.mean(dnbr_clean):.3f}, range={np.min(dnbr_clean):.3f} to {np.max(dnbr_clean):.3f}")
    
    # Burn severity stats
    severity = get_burn_severity_stats(dnbr_clean)
    print(f"Burn Severity ({severity['total']} pixels):")
    print(f"  High (>0.66):       {severity['high'][0]:4d} ({severity['high'][1]:4.1f}%)")
    print(f"  Moderate-High:      {severity['mod_high'][0]:4d} ({severity['mod_high'][1]:4.1f}%)")
    print(f"  Moderate-Low:       {severity['mod_low'][0]:4d} ({severity['mod_low'][1]:4.1f}%)")
    print(f"  Low (0.1-0.25):     {severity['low'][0]:4d} ({severity['low'][1]:4.1f}%)")
    print(f"  Unburned (<0.1):    {severity['unburned'][0]:4d} ({severity['unburned'][1]:4.1f}%)")
    
    # Get elevation data
    print(f"\n=== GETTING ELEVATION DATA ===")
    elevation = get_elevation_points(fire.geometry)
    if elevation is None:
        print("Skipping elevation analysis")
        return {'fire': fire, 'dnbr': dnbr, 'severity': severity}
    
    # Elevation analysis
    print(f"\n=== ANALYZING ELEVATION VS SEVERITY ===")
    elev_results = analyze_elevation_vs_severity(dnbr, elevation, fire.geometry)
    
    # Print detailed results
    print("======================================================================")
    print("FIRE SEVERITY + ELEVATION ANALYSIS RESULTS")
    print("======================================================================")
    
    print(f"Fire Information:")
    print(f"Name: {fire['INCIDENT']}")
    print(f"ID: {fire.get('UNQE_FIRE_', 'Unknown')}")
    print(f"Year: {int(fire['year'])}")
    print(f"Fire Date: {fire_date_parsed.strftime('%Y-%m-%d')}")
    print(f"Post-fire delay: {post_fire_years} years")
    print(f"Satellite: Landsat 5 (30m resolution)")
    
    print(f"Spatial Analysis:")
    print(f"Total analyzed pixels: {len(elev_results['dnbr']):,}")
    print(f"Fire area: {fire.geometry.area:.8f} square degrees")
    
    print(f"Elevation Summary:")
    print(f"Range: {elev_results['elevation'].min():.0f} - {elev_results['elevation'].max():.0f} m")
    print(f"Mean: {elev_results['elevation'].mean():.0f} m")
    print(f"Relief: {elev_results['elevation'].max() - elev_results['elevation'].min():.0f} m")
    
    print(f"Fire Severity Summary:")
    print(f"dNBR Range: {elev_results['dnbr'].min():.3f} - {elev_results['dnbr'].max():.3f}")
    print(f"Mean dNBR: {elev_results['dnbr'].mean():.3f}")
    
    print(f"Elevation-Fire Severity Correlation:")
    print(f"Pearson correlation: {elev_results['correlation']:.4f}")
    print(f"P-value: {elev_results['p_value']:.6f}")
    print(f"Statistical significance: {'Yes' if elev_results['p_value'] < 0.05 else 'No'} (α=0.05)")
    
    # Interpretation
    corr_abs = abs(elev_results['correlation'])
    if corr_abs > 0.3:
        strength = "strong" if corr_abs > 0.7 else "moderate"
    else:
        strength = "weak"
    
    direction = "positive" if elev_results['correlation'] > 0 else "negative"
    
    print(f"Interpretation:")
    print(f"There is a {strength} {direction} correlation between elevation and fire severity.")
    
    if elev_results['correlation'] > 0.1:
        print(f"Higher elevations experienced more severe burning on average.")
    elif elev_results['correlation'] < -0.1:
        print(f"Lower elevations experienced more severe burning on average.")
    else:
        print(f"Elevation does not strongly influence fire severity in this fire.")
    
    # Elevation-based severity table
    bins_df = elev_results['elevation_bins']
    print(f"Detailed Severity by Elevation Zones:")
    print(f"{'Elevation (m)':<12} {'Mean dNBR':<10} {'High Sev %':<11} {'Mod-High %':<11} {'Pixels':<8}")
    print("------------------------------------------------------------------")
    for _, row in bins_df.iterrows():
        print(f"{row['elevation']:<12.0f} {row['mean_dnbr']:<10.3f} "
              f"{row['high_severity']:<11.1f} {row['moderate_high']:<11.1f} {row['pixels']:<8.0f}")
    
    # Make plots
    print(f"\n=== CREATING PLOTS ===")
    plot_results(elev_results, f"{fire['INCIDENT']} ({fire_date_parsed.strftime('%Y')})")
    
    print(f"Analysis complete!")
    
    return {
        'fire': fire,
        'dnbr': dnbr,
        'elevation': elevation,
        'severity': severity,
        'elevation_analysis': elev_results
    }

def run_fire_elevation_analysis_by_date(fire_name, fire_date, post_fire_years=1, fire_file="../subsetted_data/mnf_fires_all.geojson"):
    """
    Complete workflow function with customizable post-fire delay
    
    Parameters:
    -----------
    fire_name : str
        Name of the fire to analyze
    fire_date : str
        Fire date in 'YYYY-MM-DD' format  
    post_fire_years : float
        Years after fire to search for post-fire imagery (default: 1.0)
    fire_file : str
        Path to fire boundaries file
        
    Returns:
    --------
    dict : Analysis results including fire info, imagery, severity stats, and plots
    
    Examples:
    ---------
    # Standard 1-year delay
    results = run_fire_elevation_analysis_by_date('Monument Rock', '1989-07-28')
    
    # Custom 2-year delay for more vegetation recovery
    results = run_fire_elevation_analysis_by_date('Canyon Creek', '1989-08-04', post_fire_years=2.0)
    
    # 6-month delay for immediate post-fire analysis
    results = run_fire_elevation_analysis_by_date('Some Fire', '1990-06-15', post_fire_years=0.5)
    """
    
    print(f"=== FIRE ELEVATION ANALYSIS WORKFLOW ===")
    print(f"Fire: {fire_name}")
    print(f"Fire Date: {fire_date}")
    print(f"Post-fire delay: {post_fire_years} years")
    print("="*50)
    
    return analyze_fire(fire_name, fire_date, post_fire_years, fire_file)

if __name__ == "__main__":
    print("Fire Analysis Elevation and Burn Severity Module - COMPLETE VERSION")
    print("===================================================================")
    print("Functions:")
    print("analyze_fire() - Main analysis function with customizable post-fire delay")
    print("run_fire_elevation_analysis_by_date() - Complete workflow function")
    print()
    print("Quick Examples:")
    print("results = analyze_fire('Canyon Creek', '1989-08-04')  # 1 year delay")
    print("results = analyze_fire('Canyon Creek', '1989-08-04', post_fire_years=2)  # 2 year delay")
    print("results = run_fire_elevation_analysis_by_date('Monument Rock', '1989-07-28', post_fire_years=1.5)")
    print()
    print("Key Features:")
    print("- Customizable post-fire delay (0.5, 1, 2, 5 years, etc.)")
    print("- Simplified search logic - searches around specific target dates")
    print("- Clean separation of pre-fire and post-fire imagery")
    print("- Comprehensive elevation vs fire severity analysis")
    print("- Six-panel visualization plots")
    print("- Statistical correlation analysis")