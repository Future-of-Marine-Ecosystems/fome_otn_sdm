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

# example of selecting variables from three datasets
datasets <- list(
  
  # Download temperature
  list(dataset_id = "tas_baseline_2000_2020_depthsurf",
       variables = c("tas_mean"),
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'), latitude = ylim, longitude = xlim)),
  
  # Download salinity
  list(dataset_id = "so_baseline_2000_2019_depthsurf",
       variables = "so_mean",
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'), latitude = ylim, longitude = xlim)),
  
  # Download chlorophyll
  list(dataset_id = "chl_baseline_2000_2018_depthsurf",
       variables = c("chl_mean"),
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'), latitude = ylim, longitude = xlim)),
  
  # Download bathymetry
  list(dataset_id = "terrain_characteristics",
       variables = c("bathymetry_mean"),
       constraints = list(time = c('1970-01-01T00:00:00Z', '1970-01-01T00:00:00Z'), latitude = ylim, longitude = xlim))
  
  ) # End datasets list

# Download rasters
for (dataset in datasets) {
  
  # Input dataset information
  dataset_id <- dataset$dataset_id
  variables <- dataset$variables
  constraints <- dataset$constraints
  
  download_layers(dataset_id, variables = variables, constraints = constraints, directory= dir)
}


# Load in and stack rasters
bio_oracle <- raster::stack(rast(paste0(dir, list.files(dir))))

# Grab number of layers
nlayers = nlayers(bio_oracle)

# Change bathymetry date to match others
bio_oracle$layers[[4]]$z = '2010-01-01'

# Translate rasters to terra
bio_oracle = terra::rast(bio_oracle)

# Set all to same time
time(bio_oracle) <- c(rep(as.Date("2010-01-01 UTC"), nlayers))

# Plot environmental variables
plot(bio_oracle)




