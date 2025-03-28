# "R for movement ecologists: Species Distribution Models
# Reid Steele & Esteban Salazar
# March 13, 2025

# Libraries
library(tidyverse)
library(data.table)
library(glatos)
library(terra)
library(tidyterra)
library(geodata)
library(raster)
library(terra)
library(biomod2)

# source functions
source('functions.R')


#############################################################################
################### 1.Format and Visualize Acoustic Data  ###################
#############################################################################


# Data directory
data_dir = './Data/'

# Load in data
detect_data = otn_read(data_dir)

# Map

# Generate x and y limits
xlim = c(floor(min(detect_data$deploy_long))-1, ceiling(max(detect_data$deploy_long))+1)
ylim = c(floor(min(detect_data$deploy_lat))-1, ceiling(max(detect_data$deploy_lat))+1)

# Set world object
NorthAmerica <- gadm(country = country_codes("North America")$ISO3,
                     level = 0, resolution = 2,
                     path = getwd())

# map method
detmap_y = ggplot() +
  geom_spatvector(data = NorthAmerica) +
  coord_sf(xlim = xlim, ylim = ylim) +
  geom_point(data = detect_data, aes(x = deploy_long, y = deploy_lat, color = detectedby))
detmap_y

# Remove false detections
detect_filt = false_detections(detect_data, 3600, show_plot = T)

# Filter out false detections
detect_filt = filter(detect_filt, passed_filter == T)


# Map again

# Generate x and y limits
xlim = c(floor(min(detect_filt$deploy_long))-1, ceiling(max(detect_filt$deploy_long))+1)
ylim = c(floor(min(detect_filt$deploy_lat))-1, ceiling(max(detect_filt$deploy_lat))+1)

# Set world object
NorthAmerica <- gadm(country = country_codes("North America")$ISO3,
                     level = 0, resolution = 2,
                     path = getwd())

# map method
detmap_f = ggplot() +
  geom_spatvector(data = NorthAmerica) +
  coord_sf(xlim = xlim, ylim = ylim) +
  geom_point(data = detect_filt, aes(x = deploy_long, y = deploy_lat, color = detectedby))
detmap_f

# Data explorations

# How many animals
length(unique(detect_filt$animal_id))

# Number of detections by receiver group
dets_rg = group_by(detect_filt, detectedby) %>% summarize(n = n())
dets_rg

# Seasonal pattern
ggplot(detect_filt, aes(x = julianday, y = deploy_lat, color = detectedby)) + geom_point()

# Calculate detection events
events = detection_events(detect_filt, time_sep = 86400)

# Species name
spp_to_model = unique(detect_data$scientificname)
spp_to_model

