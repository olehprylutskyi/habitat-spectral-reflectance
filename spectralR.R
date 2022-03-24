#### Installation of required components and startup configuration ####

# This code is aimed to obtain, process, and visualize spectral reflectance data for
# the user-defined land(water) surface classes, for visual exploring in which
# wavelength the classes differ. Input should be a shapefile with polygons
# of surface classes (it might be different habitat types, crops, any other units). The 
# single source of spectral data so far is Sentinel2 L2A satellite mission optical bands
# pixel data, obtained through Google Earth Engine service.

# The workflow depends on rgee R package, which provides a bridge between R and Python 
# API for Google Earth Engine. All the operations with satellite imageries run in a cloud,
# and afterwards obtained pixel data visualize locally. Therefore, despite of extent of 
# input data, the most resource hungry operations do not load your local machine.

# For using rgee you should have a Google Earth Engine account. 
# If not, first register her using your Google account: https://earthengine.google.com/new_signup/

# Essential requirements:
# - stable Internet connection (for using API)
# - Installed and correctly configured Python environment (v. 3.5 or above)
# - active Google Earth Engine account

# Depends on operating system you use and your current Python configuration, it may require some 
# additional packages. See the following links for instructions.
# Quick Start User's Guide for rgee: https://www.rdocumentation.org/packages/rgee/versions/1.0.7
# Official documentation for rgee: https://r-spatial.github.io/rgee/index.html 
# rgee source code: https://github.com/r-spatial/rgee

# We strongly encourage you to follow official rgee installation guide and messages
# arrived during installation process


# Install rgee
remotes::install_github("r-spatial/rgee")

# Load the library
library(rgee)

## It is necessary just once to complete installation necessary dependencies
ee_install()
# If something went wrong in this step, see https://r-spatial.github.io/rgee/index.html#installation

# Check non-R dependencies
ee_check() 

# Initialize Google Earth Engine API
ee_Initialize()

#### Environment preparation ####

# Reset R's brain before new analysis session started
rm(list = ls())

# Config environment
library(tidyverse)
library(rgee)
library(sf)
library(geojsonio)
library(reshape2)

# Initialize Google Earth Engine API for current session
ee_Initialize()

#### Upload and process vector data ###

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


#### Obtain pixel values from Sentinel 2A image collection ####

ee_nc <-  sf_as_ee(nc) # convert sf to ee featureCollection object

# create an envelope region of interest to filter image collection
region <- ee_nc$geometry()$bounds()


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
	filterDate("2019-05-15", "2019-06-30")$
	filterBounds(region)$
  map(maskS2clouds)$
	select(c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"))$
	median()


# This property of the table stores the land cover labels.
label <- "class"

# Overlay the polygons (or any other vector features) on the imagery to get training.
training <- sentinel2A$sampleRegions(
  collection = ee_nc,
  properties = list(label),
  scale = 100
)

# Convert training to the sf object
values <-  ee_as_sf(training,
                  maxFeatures = 100000)


# Get final dataframe with class labels and pixel values
reflectance <-  values %>%
  left_join(classes_cheatsheet, by="class") %>%
  st_drop_geometry() %>%
  select(-class)

save(reflectance, file = "reflectance_data")

load(file = "./reflectance_data")


#### Making spectral reflectance curves ####

# Quantitative overview of pixel data
summary(factor(reflectance$label))


# Number of spectral values for different classes

ggplot(reflectance, aes(x=label))+
  geom_bar()+
  labs(x = 'Classes of surface', y = 'Number of pixels',
       title = "Total pixel data for different classes of surface",
       caption = 'Data: Sentinel-2 Level-2A')+
  theme_minimal()


## Prepare pixel data to the ready-to-visualize mode

# Create "dummy" wavelength object, containing mean wavelengths (nm) for Sentinel 2A 
# (https://en.wikipedia.org/wiki/Sentinel-2), for bands 2-12

dummy_wavelength <-  c(492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4)
bands <-  c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12")
waves <-  cbind(bands, dummy_wavelength)
colnames(waves)[1] <-  "variable"


# Reshape the dataframe to make it appropriate to ggplot2 syntax
df1 <- tibble::as_tibble(reflectance) %>%
  reshape2::melt(id = "label") %>%
  left_join(as.data.frame(waves)) %>%
  mutate(across(label, as.factor)) %>%
  mutate(across(dummy_wavelength, as.numeric)) %>%
  mutate(across(variable, as.factor)) %>%
  mutate(across(value, as.numeric)) %>%
  mutate(variable = factor(variable, 
                           levels = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"))) %>%
  na.omit()


# Spectral reflectance curves for all habitat types
#png("SRC_all_classes.png", width = 200,
#    height = 170, units = "mm", res = 150)
ggplot(df1, aes(x=dummy_wavelength, y= value, colour = label))+
  geom_smooth(aes(fill = label))+
  labs(x = 'Wavelength, nm', y = 'Reflectance',
       colour = "Surface classes",
       fill = "Surface classes",
       title = "Spectral reflectance curves for different classes of surface",
       caption = 'Data: Sentinel-2 Level-2A')+
  theme_minimal()
#dev.off()

# Transform the data to acquire statistical summaries (mean, mean-standard deviation, 
# mean+standard deviation) for each group of variables
# Summarize
df2 <-  df1 %>%
  drop_na() %>%
  mutate(across(variable, as.factor)) %>%
  mutate(across(value, as.numeric)) %>%
  group_by(variable, label) %>% 
  summarise(mean_refl = mean(value), min_refl = mean(value)-sd(value), max_refl = mean(value)+sd(value)) %>%
  left_join(as.data.frame(waves)) %>%
  mutate(across(dummy_wavelength, as.numeric)) %>%
  rename(band = variable, wavelength = dummy_wavelength) %>%
  mutate(band = factor(band, 
                       levels = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")))



# Reflectance ber bands on high classification level
# png("Reflectance_per_bands.png", width = 200,
#    height = 150, units = "mm", res = 150)
ggplot(df2, aes(x=band, y=mean_refl, colour = label))+
  geom_line(aes(group = label))+
  geom_pointrange(aes(ymin = min_refl, ymax = max_refl), width = 0.2)+
  labs(x = 'Sentinel-2 bands', y = 'Reflectance',
       colour = "Surface classes",
       title = "Reflectance for different surface classes",
       caption='Data: Sentinel-2 Level-2A\nmean ± standard deviation')+
  theme_minimal()
# dev.off()


# Differences between habitat types per channel
#png("Violin_all_variables_EUNIS.png", width = 250,
#    height = 210, units = "mm", res = 150)
ggplot(df1, aes(x=label, y= value, fill = label))+
  geom_violin(trim = FALSE)+
  facet_wrap( ~ variable, ncol = 2)+
  labs(x='Surface class',y='Reflectance',
       fill="Surface classes",
       title = "Reflectance for different surface classes",
       caption='Data: Sentinel-2 Level-2A')+
  theme_minimal()
#dev.off()

#### The end ####
