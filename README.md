# **🔥 CHAR: Climate, Herbivory, and Risk 🌲🐛**
## **Investigating the Role of Spruce Budworm Outbreaks in Fire Vulnerability**

### **Team Members**
- Rachel Potter  
- Iris Mire  

---

## **Project Summary**
Spruce budworm (**Choristoneura occidentalis**) is a forest pest that weakens trees, increasing their susceptibility to wildfires. In this project, we examine the spatial dynamics of **Western spruce budworm outbreaks** and their role in **fire vulnerability** in the Pacific Northwest. Using geospatial analysis, climate data, and fire behavior modeling, we aim to identify key predictors of fire risk in infested forests.

🌍 **Range Map of the Western Spruce Budworm:**  
![Spruce Budworm Range](range_map.png)

---

## **Background**
Western spruce budworms are defoliators that consume conifer needles, particularly in **spruce and fir forests**. Over several years of infestation (3–10 years), trees experience stress due to reduced foliage, weakening their defenses and making them more flammable. Understanding how these outbreaks interact with climate variability and fire regimes is crucial for forest management and wildfire mitigation.

---

## **Questions & Objectives**
🔍 **Key Research Questions:**
- How do spruce budworm outbreaks influence fire vulnerability in Pacific Northwest forests?
- What spatial and temporal patterns exist in budworm outbreaks, and how do they correlate with climate variability? 
- Can fire behavior models predict increased fire risk in budworm-affected forests? 

Our objective is to analyze historical outbreak patterns, integrate climate and fire data, and assess fire vulnerability using geospatial approaches.

---

## **Datasets 📊**
We will use the following datasets:
- **🌲 Forest Health Monitoring Program (USDA Forest Service)** – [Link](https://www.fs.usda.gov/foresthealth/) (Budworm outbreak data)
- **🏞️ State Department of Natural Resources** – Regional forest health reports
- **🛰️ NASA MODIS (Terra/Aqua)** – [Link](https://modis.gsfc.nasa.gov/) (Satellite imagery for vegetation and fire detection)
- **🌡️ NOAA Climate Data** – [Link](https://www.ncdc.noaa.gov/) (Temperature, precipitation, and drought indices)
- **🔥 National Interagency Fire Center (NIFC)** – [Link](https://www.nifc.gov/) (Fire occurrence and severity data)
- **☁️ ECMWF Climate Reanalysis** – [Link](https://www.ecmwf.int/en/forecasts/datasets) (Long-term climate variability)

---

## **Tools & Packages 🛠️** 
We will conduct all analyses in **Python 🐍**, using the following libraries:
- **Geospatial Analysis**: `geopandas`, `rasterio`, `shapely`
- **Data Processing & Visualization**: `pandas`, `numpy`, `matplotlib`, `seaborn`
- **Machine Learning & Modeling**: `scikit-learn`, `xgboost`
- **Remote Sensing**: `Google Earth Engine (gee)`, `satpy`
- **Fire Behavior Modeling**: `FARSITE`, `FlamMap`

---

## **Methodology 📈 **
1. **Data Collection & Cleaning**: Aggregate outbreak, fire, and climate data, ensuring temporal and spatial alignment.
2. **Geospatial Analysis**: Use remote sensing to map budworm infestations over time and integrate fire occurrence data.
3. **Climate Correlation**: Analyze climate trends and their relationship with outbreak severity.
4. **Fire Vulnerability Modeling**: Use fire behavior models and statistical techniques (e.g., Random Forest, regression) to assess fire risk in affected areas.
5. **Visualization & Interpretation**: Create maps, graphs, and statistical summaries to communicate findings.

---

## **Expected Outcomes**
🎯 **Goals we hope to achieve:**
- A comprehensive geospatial analysis of **budworm outbreaks** and their role in fire susceptibility.
- Identification of key climatic and ecological factors that drive **fire risk** in infested forests.
- Predictive models that could help **inform forest management and wildfire prevention strategies**.
- Open-source dataset and visualization tools for further research.

---

## **References**
📚 **Sources & Additional Reading:**
- USDA Forest Service: [https://www.fs.usda.gov/foresthealth/](https://www.fs.usda.gov/foresthealth/)
- MODIS Fire and Vegetation Data: [https://modis.gsfc.nasa.gov/](https://modis.gsfc.nasa.gov/)
- NOAA Climate Data: [https://www.ncdc.noaa.gov/](https://www.ncdc.noaa.gov/)
- National Interagency Fire Center: [https://www.nifc.gov/](https://www.nifc.gov/)
- ECMWF Climate Reanalysis: [https://www.ecmwf.int/en/forecasts/datasets](https://www.ecmwf.int/en/forecasts/datasets)

📖 **Relevant Academic Articles:**
- *Fettig, C. J., et al. (2007). "Effects of climate change on insects and pathogens in coniferous forests of the western United States and Canada." Environmental Reviews, 15(1), 1-17.* [DOI](https://doi.org/10.1139/a06-016)
- *Hummel, S., & Agee, J. K. (2003). "Western spruce budworm defoliation effects on fuel dynamics and potential fire behavior in mixed-conifer forests." Forest Ecology and Management, 186(1-3), 13-28.* [DOI](https://doi.org/10.1016/S0378-1127(03)00228-2)
- *Raffa, K. F., et al. (2008). "Cross-scale drivers of natural disturbances prone to anthropogenic amplification: The dynamics of bark beetle eruptions." BioScience, 58(6), 501-517.* [DOI](https://doi.org/10.1641/B580607)

