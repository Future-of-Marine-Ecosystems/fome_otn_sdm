# Load in environmental data from copernicus
# Reid Steele and Esteban Salazar
# March 20 2025

# Load required libraries
library(CopernicusMarine)
library(raster)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(lubridate)

# Set file download location and which block to download (max = 13)
# USERS NEED TO CHANGE THESE
file_download_location = "./Data/copernicus/"

# Create directories to download files and save extracted data to
# dir.create(file.path(file_download_location))
# dir.create(file.path(file_save_location))

# Select the data set
tibble_files_list <- cms_list_stac_files("GLOBAL_MULTIYEAR_PHY_001_030", "cmems_mod_glo_phy_my_0.083deg_P1M-m")

# Filter out unwanted dates
filt = which((as.numeric(substr(tibble_files_list$current_path, 114, 117)) >= 2014) & 
               (as.numeric(substr(tibble_files_list$current_path, 118, 119)) >= 5))
tibble_files_list = tibble_files_list[filt,]

print(nrow(tibble_files_list))

# set layers_to_extract to etag
layers_to_extract = rownames(tibble_files_list)
# Zero pad
layers_to_extract = ifelse(nchar(layers_to_extract) < 4, paste0(0, layers_to_extract), layers_to_extract)
layers_to_extract = ifelse(nchar(layers_to_extract) < 4, paste0(0, layers_to_extract), layers_to_extract)
layers_to_extract = ifelse(nchar(layers_to_extract) < 4, paste0(0, layers_to_extract), layers_to_extract)
layers_to_extract = ifelse(nchar(layers_to_extract) < 4, paste0(0, layers_to_extract), layers_to_extract)
layers_to_extract = ifelse(nchar(layers_to_extract) < 4, paste0(0, layers_to_extract), layers_to_extract)

# Set days to extract
# Silky dates correspond to  4,384 to 10,500... It only downloads until 10,408 which corresponds to 2021-06-30

# Download rasters
for (ii in 1:nrow(tibble_files_list))
{
  # Download data and convert to a raster
  cms_download_stac(tibble_files_list[ii,, drop = FALSE], destination= file_download_location, show_progress = F, overwrite = T)
  print(ii)
}


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

