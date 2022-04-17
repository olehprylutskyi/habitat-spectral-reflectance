# Spectral reflectance of habitat types: case study on Southern Bug valley, Ukraine

This bunch of scripts are deprecated. Use instead our integrated R package [spectralR](https://github.com/olehprylutskyi/spectralR), which does the same, but avoid downloading large files to the local machine, as well as using multiple tools.

## Description
Scripts for Google Eart Engine (javascript) and R (as RMarkdown) for obtaining and visualizing satellite derived spectral reflectance data for predefined polygons (habitat types in this case).

GEE is used for obtaining cloud-masked, pre-processed, 10-meters in resolution, multi-band surface reflectance satellite imagery for user-defined area, medianized for user-defined timespan from the Sentinel-2 Level 2A imagery collection, available through Google Earth Engine platform.

RMarkdown is used for obtain pixel reflectance data for user-defined ground truth (as polygon ESRI shapefiles), and produce a range of visualisation for the purpose of exploratory analysis of spectral reflectance of studied types of land cover (EUNIS habitat types in our case). Of course, the scripts can be used for other purposes, boiled down to quick visual assessment of spectral reflectance of known spatial polygons.


## Requirements
- Google Earh Engine web service
- R 4.0 or higher
- RStudio IDE

## Outputs
- 10-meter in resolution medianized satellite imagery
- R-objects with the raw pixel data for furter analysis
- Spectral reflectance curves for user-defined set of habitat types
- Reflectance per bands for user-defined set of habitat types
- Density of values (as violin plots) for selected habitat types per band
- NDVI as pixel data for selected period

## Output examples
![Spectral reflectance curves for user-defined set of habitat types](https://github.com/olehprylutskyi/habitat-spectral-reflectance/blob/main/SRC_within_bands_S35_S36_R1B_X18_T19_2019_15May-15Jun.png)
![Reflectance per bands for user-defined set of habitat types](https://github.com/olehprylutskyi/habitat-spectral-reflectance/blob/main/Reflectance_within_bands_S35_S36_R1B_X18_T19_2019_15May-15Jun.png)
![Density of values (as violin plots) for selected habitat types per band](https://github.com/olehprylutskyi/habitat-spectral-reflectance/blob/main/Violins_R_2019_15May-15Jun.png)

## Known issues
1. Timespan selected manually, though using loops would be better.

## References
Shyriaieva, D., Prylutskyi, O. (2021). Exploratory analysis of the spectral reflectance curves of habitat types: a case study on Southern Bug River valley, Ukraine. In: 63rd IAVS Annual Symposium: Book of Abstracts, p. 153.
