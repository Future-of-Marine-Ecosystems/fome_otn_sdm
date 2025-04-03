# Load in environmental data from copernicus
# Reid Steele and Esteban Salazar
# March 20 2025

# Load required libraries
library(raster)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(lubridate)
library(biooracler)


#############################################################################
##################### 2. Make Environmental Data Rasters ####################
#############################################################################

# Bio-ORACLE

# Code adapted from https://github.com/bio-oracle/nc-multiple-layers-tutorial/blob/main/downloadMultipleDatasets.Rmd

# Download directory
dir = './Data/bio-oracle/'

# Create the directory where we will store our data if it does not already exist
if(dir.exists(dir) == F){dir.create(dir)}

# Lets look at available Bio-ORACLE Layers
layers = biooracler::list_layers()
View(layers)

# example of selecting variables from four datasets
datasets <- list(
  
  # Download temperature
  list(dataset_id = "tas_baseline_2000_2020_depthsurf", # Bio-ORACLE dataset name
       variables = c("tas_mean"), # Variable name(s)
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'), # Must have two time coordinates, can be identical
                          latitude = ylim, longitude = xlim)), # Space Coordinates
  
  # Download salinity
  list(dataset_id = "so_baseline_2000_2019_depthsurf", # Bio-ORACLE dataset name
       variables = "so_mean", # Variable name(s)
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'), # Must have two time coordinates, can be identical
                          latitude = ylim, longitude = xlim)), # Space Coordinates
  
  # Download chlorophyll
  list(dataset_id = "chl_baseline_2000_2018_depthsurf", # Bio-ORACLE dataset name
       variables = c("chl_mean"), # Variable name(s)
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'), # Must have two time coordinates, can be identical
                          latitude = ylim, longitude = xlim)), # Space Coordinates
  
  # Download bathymetry
  list(dataset_id = "terrain_characteristics", # Bio-ORACLE dataset name
       variables = c("bathymetry_mean"), # Variable name(s)
       constraints = list(time = c('1970-01-01T00:00:00Z', '1970-01-01T00:00:00Z'), # Must have two time coordinates, can be identical
                          latitude = ylim, longitude = xlim)) # Space Coordinates
  
  ) # End datasets list

# Download rasters
for (dataset in datasets) {
  
  # Input dataset information
  dataset_id <- dataset$dataset_id
  variables <- dataset$variables
  constraints <- dataset$constraints
  
  # Download layers into our directory
  download_layers(dataset_id, variables = variables, constraints = constraints, directory = dir)
}



# Now we need to combine the rasters

# Load in and stack rasters
bio_oracle <- raster::stack(rast(paste0(dir, list.files(dir))))

# Grab number of layers
nlayers = nlayers(bio_oracle)
nlayers # Should be equal to our number of variables

# Translate rasters to terra for use with biomod2
bio_oracle = terra::rast(bio_oracle)

# Set all rasters to the same time (Bathymetry is 1970)
time(bio_oracle) <- c(rep(as.Date("2010-01-01 UTC"), nlayers))

# Plot environmental variables
plot(bio_oracle)

########## Return to presentation ##########


