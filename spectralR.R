#### Description ####

# This code is aimed to obtain, process, and visualize spectral reflectance data for
# the user-defined land(water) surface classes, for visual exploring in which
# wavelength the classes differ. Input should be a shapefile with polygons
# of surface classes (it might be different habitat types, crops, any other things). The 
# single source of spectral data are Sentinel2 L2A satellite mission optical bands
# pixel data so far, obtained through Google Earth Engine service.

# The workflow depends on rgee R package, which provides a bridge between R and Python 
# API for Google Earth Engine. All the operations with satellite imageries run in a cloud,
# and afterwards obtained pixel data visualize locally. Therefore, despite of extent of 
# input data, the most resource hungry operations do not overload your local machine. But that means that
# you need a stable Internet connection for using API.

# For using rgee you should have a Google Earth Engine account. 
# If you don't, first register her using your Google account: https://earthengine.google.com/new_signup/

# Depends on operating system you use and your current Python configuration, it may require some 
# additional R and Python packages for running rgee. See the following links for instructions.
# Quick Start User's Guide for rgee: https://www.rdocumentation.org/packages/rgee/versions/1.0.7
# Official documentation for rgee: https://r-spatial.github.io/rgee/index.html 
# rgee source code: https://github.com/r-spatial/rgee

# We strongly encourage you to follow official rgee installation guide and messages
# arrived during installation process.

# The overall workflow is following:
# 1. Load user's ESRI shapefile containing polygons for user-defined surface classes,
# as well as the text or numerical field with classes names (labels).
# 2. Apply rgee functionality to retrieve multi-band pixel data for classes polygons from 
# Google Earth Engine service.
# 3. Visualize retrieved pixel data locally, mainly using ggplot2 approach.

# Essential requirements:
# - stable Internet connection (for using API)
# - Installed and correctly pre-configured Python environment (v. 3.5 or above)
# - valid Google Earth Engine account

# Install and set up rgee
remotes::install_github("r-spatial/rgee")

# Load the library
library(rgee)

## It is necessary just once to complete installation necessary dependencies
ee_install()
# If something went wrong in this step, see https://r-spatial.github.io/rgee/index.html#installation

# Check non-R dependencies
ee_check() 

# rgee developers recommend installing the version of the Earth Engine Python API
# which rgge was tested with, using the following command
ee_install_upgrade()

# Initialize Google Earth Engine API
ee_Initialize()

# If everything is OK on this step and you see a message of successful initiation of 
# API with your GEE username in console, - congratulations, you managed to install 
# and configure rgee!

#### Use case 1. Basic habitat types of 'Homilsha Forests' National Park neighborhood. ####

#### Environment preparation ####

# Reset R's brain before new analysis session started
rm(list = ls())

# Load required packages
library(tidyverse)
library(rgee)
library(sf)
library(geojsonio)
library(reshape2)

# Initialize Google Earth Engine API for current session
ee_Initialize()

#### Upload and process vector data ####

### Step-by-step ###

# Upload a shapefile with polygons of knowing surface classes.
# The shapefile must contain a text field with classes labels
nc <-  st_read("test_shapefile.shp", quiet = TRUE)

# Specify which field contains classes labels
label <-  "veget_type"

# Create a 'class' variable for integer class IDs
names(nc)[names(nc) == label] <- 'label'

# make a list of label values (types of surface) and its numerical IDs
classes_cheatsheet <- as.data.frame(levels(factor(nc$label)))
classes_cheatsheet$class <- rownames(as.data.frame(levels(factor(nc$label))))
colnames(classes_cheatsheet) <- c("label", "class")
classes_cheatsheet <-  classes_cheatsheet %>%
  mutate(across(label, as.factor)) %>%
  mutate(across(class, as.numeric))

# Add class IDs
nc <- left_join(nc, classes_cheatsheet, by = "label")
nc$class <- as.numeric(nc$class)
# delete objects with NA in the target variable
nc <- nc[!is.na(nc$label) ,]

rm(label) # drop spare variables to lighten your R Environment

### Using function ###

# Function takes shapefile with polygons of different classes of surface (habitats,
# crops, vegetation, etc.), and retrieves ready-for-sampling sf object. More about sf
# (), one of the basic R spatial data format, here: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Specify:
# - shapefile name (if should lie within working directory, using absolute 
# paths were not tested)
# - name of the field which contains class labels

prepare.vector.data <- function(shapefile_name, label_field){
  # Upload a shapefile with polygons of knowing surface classes.
  # The shapefile must contain a text field with classes labels
  nc <-  st_read(shapefile_name, quiet = TRUE)
  
  # Specify which field contains classes labels
  # label <-  label_field
  
  # Create a 'class' variable for integer class IDs
  names(nc)[names(nc) == label_field] <- 'label'
  
  # make a list of label values (types of surface) and its numerical IDs
  classes_cheatsheet <- as.data.frame(levels(factor(nc$label)))
  classes_cheatsheet$class <- rownames(as.data.frame(levels(factor(nc$label))))
  colnames(classes_cheatsheet) <- c("label", "class")
  classes_cheatsheet <-  classes_cheatsheet %>%
    mutate(across(label, as.factor)) %>%
    mutate(across(class, as.numeric))
  
  # Add class IDs
  nc <- left_join(nc, classes_cheatsheet, by = "label")
  nc$class <- as.numeric(nc$class)
  # delete objects with NA in the target variable
  nc <- nc[!is.na(nc$label) ,]
}


sf_df <- prepare.vector.data("test_shapefile.shp", "veget_type")

#### Obtain pixel values from Sentinel 2A image collection ####

### Step-by-step ###

ee_df <-  sf_as_ee(sf_df) # convert sf to ee featureCollection object

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

# Set cloud_threshold applied for individual imageries
cloud_treshold = 10

# Make median multi-band image from Sentinel L2A image collection for given
# date range and region

sentinel2A <-  ee$ImageCollection("COPERNICUS/S2_SR")$
	filterDate("2019-05-15", "2019-06-30")$
	filterBounds(region)$
  filter(ee$Filter$lt('CLOUDY_PIXEL_PERCENTAGE', cloud_treshold))$
  map(maskS2clouds)$
	select(c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"))$
	median()


# This property of the table stores the land cover labels.
label <- "class"

# Overlay the polygons (or any other vector features) on the imagery to get training.
training <- sentinel2A$sampleRegions(
  collection = ee_df,
  properties = list(label),
  scale = 100
)

# Convert training to the sf object
values <-  ee_as_sf(training,
                  maxFeatures = 10000000000)

# make a list of label values (types of surface) and its numerical IDs
classes_cheatsheet <- as.data.frame(levels(factor(df$label)))
classes_cheatsheet$class <- rownames(as.data.frame(levels(factor(df$label))))
colnames(classes_cheatsheet) <- c("label", "class")
classes_cheatsheet <-  classes_cheatsheet %>%
  mutate(across(label, as.factor)) %>%
  mutate(across(class, as.numeric))

# Get final dataframe with class labels and pixel values
reflectance <-  values %>%
  left_join(classes_cheatsheet, by="class") %>%
  st_drop_geometry() %>%
  select(-class)

save(reflectance, file = "reflectance_data")

#load(file = "./reflectance_data")



### Using function ###

# Function takes sf polygon object, obtained in the previous step, and retrieves data 
# frame with brightness values for each pixel intersected with polygons, for each 
# optical band of Sentinel-2 sensor, marked according to the label of surface class from
# the polygons.

# Specify:
# - polygons of surface classes as a sf object, prepared on previous step;
# - starting day for Sentinel image collection, as "YYYY-MM-DD". See Note 1 below;
# - final day for Sentinel image collection, as "YYYY-MM-DD";
# - cloud threshold (maximum per cent of cloud-covered pixels per image) by which individual 
# satellite imageries will be filtered;
# - scale of resulting satellite images in meters (pixel size). See Note 2 below;

# Resulting pixel data will be saved within working directory and can be loaded during
# next sessions.

# Note 1.
# Particular satellite imagery is typically not ready for instant analysis - it contains clouds, 
# cloud shadows, aerosols, as well as may cover not all the territory you of 
# your interest. Another issue is that each particular pixel slightly differs in reflectance 
# between images taken in different days due to difference in atmospheric conditions and angle
# of sunlight at the moments images were taken. Google Earth Engine has its own build-in
# algorithms for image pre-processing, atmospheric corrections and mosaicing, which allows to 
# obtain a ready-to-use, rectified image. Approach used in this script is that to find a median
# value for each pixel between several images within each of 10 optical band, and thereby make a 
# composite image. To define a set of imageries between which we are going to calculate median, 
# we should set a timeframe definining starting and final days. Sentinel-2 apparatus takes picture 
# once a 5 days, so if you set up month-long timesnap, you can expect that each pixel value will
# be calculated based on 5 to 6 values.

# Note 2.
# Finest resolution for Sentinel data - 10 m, while using larger scale values decreases 
# required computational resources and size of resulting dataframe. Although sampling 
# satellite data performs in a cloud, there are some limitations for geocalculatons placed 
# by GEE itself. If you are about to sample large areas, consider setting higher 'scale' 
# value (100, 1000).
# More about GEE best practices: https://developers.google.com/earth-engine/guides/best_practices


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
  
  # Convert training to the sf object
  values <-  ee_as_sf(training,
                      maxFeatures = 10000000000)
  
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


reflectance = get.pixel.data(sf_df, "2019-05-15", "2019-06-30", 10, 100)

save(reflectance, file = "reflectance_data") # save pixel data for further sessions


#### Making spectral reflectance curves ####

# Quantitative overview of pixel data
load(file = "./reflectance_data") # restore prevoiusly saved pixel data

summary(factor(reflectance$label)) # how many pixels in each class? 


# Number of spectral values for different classes

ggplot(reflectance, aes(x=label))+
  geom_bar()+
  labs(x = 'Classes of surface', y = 'Number of pixels',
       title = "Total pixel data for different classes of surface",
       caption = 'Data: Sentinel-2 Level-2A')+
  theme_minimal()


## Prepare pixel data to the ready-to-visualize mode

# Function takes reflectance data as an input and retrieve ggplot object, which can be 
# visualized directly or customized using ggplot functionality. Default aesthetic is 
# smoother curve (geom_smooth), method - gam. See https://ggplot2.tidyverse.org/reference/geom_smooth.html 
# for more details.

# Specify:
# - reflectance data as dataframe with pixel values for Sentinel optical bands
# B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12;


spectral.curves.plot <- function(data){
  # Create "dummy" wavelength object, containing mean wavelengths (nm) for Sentinel 2A 
  # (https://en.wikipedia.org/wiki/Sentinel-2), for bands 2-12
  
  dummy_wavelength <-  c(492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4)
  bands <-  c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12")
  waves <-  cbind(bands, dummy_wavelength)
  colnames(waves)[1] <-  "variable"

  # Reshape the dataframe to make it appropriate to ggplot2 syntax
  p <- tibble::as_tibble(data) %>%
    reshape2::melt(id = "label") %>%
    left_join(as.data.frame(waves)) %>%
    mutate(across(label, as.factor)) %>%
    mutate(across(dummy_wavelength, as.numeric)) %>%
    mutate(across(variable, as.factor)) %>%
    mutate(across(value, as.numeric)) %>%
    mutate(variable = factor(variable, 
                             levels = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"))) %>%
    na.omit() %>%
    ggplot(aes(x=dummy_wavelength, y= value, colour = label))+
    geom_smooth(aes(fill = label))
}

p1 <- spectral.curves.plot(reflectance)

p1

p1+
  labs(x = 'Wavelength, nm', y = 'Reflectance',
       colour = "Surface classes",
       fill = "Surface classes",
       title = "Spectral reflectance curves for different classes of surface",
       caption = 'Data: Sentinel-2 Level-2A')+
  theme_minimal()


# Transform the data to acquire statistical summaries (mean, mean-standard deviation, 
# mean+standard deviation) for each group of variables

# Function takes reflectance data as an input and retrieve ggplot object, which can be 
# visualized directly or customized using ggplot functionality. Default aesthetics are 
# line with statistical summary for each satellite band (geom_line + geom_pointrange),
# See https://ggplot2.tidyverse.org/reference/geom_linerange.html 
# https://ggplot2.tidyverse.org/reference/geom_path.html for more details.

# Specify:
# - reflectance data as dataframe with pixel values for Sentinel optical bands
# B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12;


stat.summary.plot <- function(data){
  # Create "dummy" wavelength object, containing mean wavelengths (nm) for Sentinel 2A 
  # (https://en.wikipedia.org/wiki/Sentinel-2), for bands 2-12
  
  dummy_wavelength <-  c(492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4)
  bands <-  c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12")
  waves <-  cbind(bands, dummy_wavelength)
  colnames(waves)[1] <-  "variable"
  
  # Reshape the dataframe to make it appropriate to ggplot2 syntax
  p <- tibble::as_tibble(data) %>%
    reshape2::melt(id = "label") %>%
    left_join(as.data.frame(waves)) %>%
    mutate(across(label, as.factor)) %>%
    mutate(across(dummy_wavelength, as.numeric)) %>%
    mutate(across(variable, as.factor)) %>%
    mutate(across(value, as.numeric)) %>%
    mutate(variable = factor(variable, 
                             levels = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"))) %>%
    na.omit() %>%
    drop_na() %>%
    mutate(across(variable, as.factor)) %>%
    mutate(across(value, as.numeric)) %>%
    group_by(variable, label) %>% 
    summarise(mean_refl = mean(value), min_refl = mean(value)-sd(value), max_refl = mean(value)+sd(value)) %>%
    left_join(as.data.frame(waves)) %>%
    mutate(across(dummy_wavelength, as.numeric)) %>%
    rename(band = variable, wavelength = dummy_wavelength) %>%
    mutate(band = factor(band, 
                         levels = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"))) %>%
    ggplot(aes(x=band, y=mean_refl, colour = label))+
    geom_line(aes(group = label))+
    geom_pointrange(aes(ymin = min_refl, ymax = max_refl), width = 0.2)
}


p2 <- stat.summary.plot(reflectance)

p2

p2 + 
  labs(x = 'Sentinel-2 bands', y = 'Reflectance',
       colour = "Surface classes",
       title = "Reflectance for different surface classes",
       caption='Data: Sentinel-2 Level-2A\nmean ± standard deviation')+
  theme_minimal()


# Reflectance per bands on high classification level

# Function takes reflectance data as an input and retrieve ggplot object, which can be 
# visualized directly or customized using ggplot functionality. Default aesthetics is 
# violin plot for each satellite band (geom_violin),
# See https://ggplot2.tidyverse.org/reference/geom_violin.html for more details.

# Specify:
# - reflectance data as dataframe with pixel values for Sentinel optical bands
# B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12;


violin.plot <- function(data){
  # Create "dummy" wavelength object, containing mean wavelengths (nm) for Sentinel 2A 
  # (https://en.wikipedia.org/wiki/Sentinel-2), for bands 2-12
  
  dummy_wavelength <-  c(492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4)
  bands <-  c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12")
  waves <-  cbind(bands, dummy_wavelength)
  colnames(waves)[1] <-  "variable"
  
  # Reshape the dataframe to make it appropriate to ggplot2 syntax
  p <- tibble::as_tibble(data) %>%
    reshape2::melt(id = "label") %>%
    left_join(as.data.frame(waves)) %>%
    mutate(across(label, as.factor)) %>%
    mutate(across(dummy_wavelength, as.numeric)) %>%
    mutate(across(variable, as.factor)) %>%
    mutate(across(value, as.numeric)) %>%
    mutate(variable = factor(variable, 
                             levels = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"))) %>%
    na.omit() %>%
    ggplot(aes(x=label, y= value, fill = label))+
    geom_violin(trim = FALSE)+
    facet_wrap( ~ variable, ncol = 2)
}

p3 <- violin.plot(reflectance)

p3

p3 + 
  labs(x='Surface class',y='Reflectance',
       fill="Surface classes",
       title = "Reflectance for different surface classes",
       caption='Data: Sentinel-2 Level-2A')+
  theme_minimal()


#### Use case 2. Habitat types of Southern Buh valley ####

#### DO NOT RUN THIS! UNDER DEVELOPMENT #

# For large data (vector file larger than 1 MB, more than 1500 pixel values to get.

# Environment preparation
rm(list = ls()) # Reset R's brain before new analysis session started

# Load required packages
library(tidyverse)
library(rgee)
library(sf)
library(geojsonio)
library(reshape2)

# You may define all the custom function mentioned above before proceeding.

# Upload and process vector data

# sf_df <- prepare.vector.data("training_polygons_shapefile.shp", "eunis_2020")
# EEException: Invalid GeoJSON geometry.

sf_df <- prepare.vector.data("half_file.shp", "eunis_2020")

# st_is_valid(sf_df) # check polygon's validity

# st_make_valid(sf_df) # fix invalid polygons if any

# Initialize Google Earth Engine API for current session
ee_Initialize(user = 'ndef', drive = TRUE)

### Step-by-step ###

ee_df <-  sf_as_ee(sf_df) # convert sf to ee featureCollection object

# Upload small geometries to EE asset
#assetId <- sprintf("%s/%s", ee_get_assethome(), 'sf_df')
#ee_df <-  sf_as_ee(sf_df,
#                   overwrite = TRUE,
#                   assetId = assetId,
#                   via = "getInfo_to_asset")


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

# Set cloud_threshold applied for individual imageries
cloud_treshold = 10

# Make median multi-band image from Sentinel L2A image collection for given
# date range and region

sentinel2A <-  ee$ImageCollection("COPERNICUS/S2_SR")$
  filterDate("2019-05-15", "2019-06-30")$
  filterBounds(region)$
  filter(ee$Filter$lt('CLOUDY_PIXEL_PERCENTAGE', cloud_treshold))$
  map(maskS2clouds)$
  select(c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"))$
  median()


# This property of the table stores the land cover labels.
label <- "class"

# Overlay the polygons (or any other vector features) on the imagery to get training.
training <- sentinel2A$sampleRegions(
  collection = ee_df,
  properties = list(label),
  scale = 10
)


# calculate total area of polygons
sf_df$area <- st_area(sf_df)
rowsum(sf_df$area)


library(raster)
area(sf_df)

# find expected number of pixels to get sampled (rows in resulting dataframe)




# Convert training to the sf object (with saving via google drive)
values <- ee_as_sf(
  training,
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





# Obtain pixel values from Sentinel 2A image collection
reflectance = get.pixel.data(sf_df, "2019-05-15", "2019-06-30", 10, 20)

# Spectral reflectance curves
p1 <- spectral.curves.plot(reflectance)

p1

p1+
  labs(x = 'Wavelength, nm', y = 'Reflectance',
       colour = "Surface classes",
       fill = "Surface classes",
       title = "Spectral reflectance curves for different classes of surface",
       caption = 'Data: Sentinel-2 Level-2A')+
  theme_minimal()

# Statistical summary plot
p2 <- stat.summary.plot(reflectance)

p2

p2 + 
  labs(x = 'Sentinel-2 bands', y = 'Reflectance',
       colour = "Surface classes",
       title = "Reflectance for different surface classes",
       caption='Data: Sentinel-2 Level-2A\nmean ± standard deviation')+
  theme_minimal()

# Violin plots
p3 <- violin.plot(reflectance)

p3

p3 + 
  labs(x='Surface class',y='Reflectance',
       fill="Surface classes",
       title = "Reflectance for different surface classes",
       caption='Data: Sentinel-2 Level-2A')+
  theme_minimal()

#### The end ####
