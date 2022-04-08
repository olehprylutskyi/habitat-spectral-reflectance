get.pixel.data <- function(sf_data, startday, endday, cloud_threshold, scale_value){
  ee_df <-  sf_as_ee(sf_data) # convert sf to ee featureCollection object

  # create an envelope region of interest to filter image collection
  region <- ee_df$geometry()$bounds()

  # cloud mask function for Sentinel-2
  # from https://github.com/ricds/DL_RS_GEE/blob/main/rgee_data_acquisition.R
  maskS2clouds <-  function(image) {
    qa = image$select('QA60');

    # Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = bitwShiftL(1,10)
    cirrusBitMask = bitwShiftL(1, 11)

    # Both flags should be set to zero, indicating clear conditions.
    mask_data = qa$bitwiseAnd(cloudBitMask)$eq(0)$And(qa$bitwiseAnd(cirrusBitMask)$eq(0));

    return(image$updateMask(mask_data)$divide(10000))
  }

  # Make median multi-band image from Sentinel L2A image collection for given
  # date range and region

  sentinel2A <-  ee$ImageCollection("COPERNICUS/S2_SR")$
    filterDate(startday, endday)$
    filterBounds(region)$
    filter(ee$Filter$lt('CLOUDY_PIXEL_PERCENTAGE', cloud_threshold))$
    map(maskS2clouds)$
    select(c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"))$
    median()

  # This property of the table stores the land cover labels.
  label <- "class"

  # Overlay the polygons (or any other vector features) on the imagery to get training.
  training <- sentinel2A$sampleRegions(
    collection = ee_df,
    properties = list(label),
    scale = scale_value
  )

  # rgee uses three different approach to upload and download data from and to the server. For small
  # dataset () default value "getInfo" is recommended, while for large vector objects / outputs
  # using intermediate container (Google Drive or Google Cloud Storage) is required. We tested
  # Google Drive approach and it showed good performance with downloading of ca. 90K pixel values.
  # function "get.pixel.data" estimates the size of input data / GEE object to be downloaded, and then
  # pick more appropriate method on their own, based on both total area of vector polygons and user defined scale
  # value (pixel size).

  if(as.numeric(sum(st_area(sf_data) / scale_value^2)) < 15000){
  # Convert training to the sf object directly
  values <-  ee_as_sf(training,
                      maxFeatures = 10000000000)
  } else {
  # Initialize Google Earth Engine API for using Google Drive as a container
  ee_Initialize(user = 'ndef', drive = TRUE)
  # Convert training to the sf object (with saving via google drive)
  values <- ee_as_sf(training,
  overwrite = TRUE,
  via = "drive",
  container = "rgee_backup",
  crs = NULL,
  maxFeatures = 10000000000,
  selectors = NULL,
  lazy = FALSE,
  public = TRUE,
  add_metadata = TRUE,
  timePrefix = TRUE,
  quiet = FALSE
)
  }

  # make a list of label values (types of surface) and its numerical IDs
  classes_cheatsheet <- as.data.frame(levels(factor(sf_df$label)))
  classes_cheatsheet$class <- rownames(as.data.frame(levels(factor(sf_df$label))))
  colnames(classes_cheatsheet) <- c("label", "class")
  classes_cheatsheet <-  classes_cheatsheet %>%
    mutate(across(label, as.factor)) %>%
    mutate(across(class, as.numeric))

  # Get final dataframe with class labels and pixel values
  reflectance <-  values %>%
    left_join(classes_cheatsheet, by="class") %>%
    st_drop_geometry() %>%
    select(-class)
}
