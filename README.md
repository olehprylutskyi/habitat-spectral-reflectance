# Spectral reflectance of habitat types: case study on Southern Bug valley, Ukraine

## Description
Scripts for Google Eart Engine (javascript) and R (as RMarkdown) for obtaining and visualizing satellite derived spectral reflectance data for predefined polygons (habitat types in this case).

GEE is used for obtaining cloud-masked, pre-processed, 10-meters in resolution, multi-band surface reflectance satellite imagery for user-defined area, medianized for user-defined timespan from the Sentinel-2 Level 2A imagery collection, available through Google Earth Engine platform.

RMarkdown is used for obtain pixel reflectance data for user-defined ground truth (as polygon ESRI shapefiles), and produce a range of visualisation for the purpose of exploratory analysis of spectral reflectance of studied types of land cover (EUNIS habitat types in our case). Of course, the scripts can be used for other purposes, boiled down to quick visual assessment of spectral reflectance of known spatial polygons.


## Requirements
- Google Earh Engine web service
- R 4.0 or higher
- RStudio IDE

# Outputs
- 10-meter in resolution medianized satellite imagery
- R-objects with the raw pixel data for furter analysis
- Spectral reflectance curves for user-defined set of habitat types
- Reflectance per bands for user-defined set of habitat types
- Density of values (as violin plots) for selected habitat types per band
- NDVI as pixel data for selected period

# Knowing issues
1. Timespan selected manually, though better using loops.
