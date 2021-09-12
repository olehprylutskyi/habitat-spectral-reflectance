/*
Script written to be executed in Google Earth Engine (https://code.earthengine.google.com/)
Scrpt extractes Sentinel-2 Level 2A satellite imagery collection for user-defined
date range and year, masks clouds by user-defined treshold, specifies required bands, generates median
values per bands, and save the result as a single tiff to Google Drive.
*/

Map.setOptions("HYBRID");

//  SET Area Of Interest (boundary for imagery)

// change path to your vector bounds of Areo Of Interest, saved as GEE assest

var AOI = ee.FeatureCollection('users/olegpril12/Southern_Buh/AOI_restricted');

var name_area = 'SBuh';

Map.centerObject(AOI, 12);


// SET DATES

var year = 2019;

var STARTDAY=196; // julian number of start day 
var ENDDAY=258; // julian number of end day

// Julian number of day obtained through the web-service
// https://people.biology.ucsd.edu/patrick/julian_cal.html

// 122-152 - (1-30 May)
// 153-182 - (1-30 June)
//
// 
// Vegetative periods

// 105-135 - (15 Apr - 15 May) 
// 135-166 - (15 May - 15 Jun)
// 166-196 - (15 Jun - 15 Jul)
// 196-227 - (15 Jul - 15 Aug)
// 227-258 - (15 Aug - 15 Sep)
// 258-288 - (15 Sep - 15 Oct)

// leaf-less preriods

// 288-319 - (15 Oct - 15 Nov)
// 074-105 - (15 Mar - 15 Apr)

var cloud_treshold=10;

function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 2 << 10;
  var cirrusBitMask = 2 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask);
//  return image.updateMask(mask).divide(10000); // divide to 10000 to make resulting raster in unit fraction.
  
}

// NDVI
var addNDVI = function(image) {
  return image.addBands(image.expression('int16((b("B8") - b("B4")) / (b("B8") + b("B4"))*1000)').rename('NDVI'));
};

// Map the function over one year of data and take the median.
// Load Sentinel-2 atmospherically corrected surface reflectance data Level-2A .
// https://explorer.earthengine.google.com/#detail/COPERNICUS%2FS2_SR
// https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR


var dataset_2a = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(AOI)
                  .filter(ee.Filter.calendarRange(year, year, 'year'))
                  .filter(ee.Filter.dayOfYear(STARTDAY, ENDDAY))
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_treshold))
                  .map(maskS2clouds)
                  .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']);
    print(year);
    print(dataset_2a);
    print('Number of Scenes', dataset_2a.size());
    
var med_dataset_2a = addNDVI(dataset_2a.median()).uint16();


var rgbVis_natural_color = {
  min: 0,
  max: 1500,
  bands: ['B4', 'B3', 'B2'],
};

// add to the map

Map.addLayer(med_dataset_2a.clip(AOI), rgbVis_natural_color, 'natural_color');

// Add AOI as vector boundary
var empty = ee.Image().byte();
var outline = empty.paint({
  featureCollection: AOI,
  color: 1,
  width: 2
});

Map.addLayer(outline, {palette: 'FF0000'}, 'AOI');

// Transform to 3-bands RGB image

var natural_color2rgb = med_dataset_2a.select(['B4', 'B3', 'B2'])
                 .divide(5.8).uint8();

// Export imagery to Google Drive

  Export.image.toDrive({
  image: med_dataset_2a.clip(AOI),
  description: name_area+'_S2_full_'+year+'_'+STARTDAY+'_'+ENDDAY,
  folder: 'GEE data',
  scale: 10,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
  });
