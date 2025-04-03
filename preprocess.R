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

# Either open the R Project in the GitHub repository or set working directory here:
setwd(getwd())


#############################################################################
################### 1.Format and Visualize Acoustic Data  ###################
#############################################################################


# Data directory
data_dir = './Data/'

# Load in data
detect_data = otn_read(data_dir)

# Create a map of detections

# Generate x and y limits
xlim = c(floor(min(detect_data$deploy_long))-1, ceiling(max(detect_data$deploy_long))+1)
ylim = c(floor(min(detect_data$deploy_lat))-1, ceiling(max(detect_data$deploy_lat))+1)

# Set world object
NorthAmerica <- gadm(country = country_codes("North America")$ISO3,
                     level = 0, resolution = 2,
                     path = getwd())

# map method
detmap_y = ggplot() +
  geom_spatvector(data = NorthAmerica) + # Map baselayer
  coord_sf(xlim = xlim, ylim = ylim) + # lat and lon limits
  geom_point(data = detect_data, aes(x = deploy_long, y = deploy_lat, color = detectedby))
detmap_y

# Some detections appear to be spatial outliers

# Remove false detections
detect_filt = false_detections(detect_data, 3600, show_plot = T)
# 3600 is arbitrary and possibly a bit high, could reduce

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
# Looks better!


# Data explorations

# How many animals in the dataset?
length(unique(detect_filt$animal_id))

# Where are the animals primarily detected?
dets_rg = group_by(detect_filt, detectedby) %>% summarize(n = n())
dets_rg

# Is there a seasonal migration pattern?
ggplot(detect_filt, aes(x = julianday, y = deploy_lat, color = detectedby)) + geom_point()

# Any issues with sampling effort?
ggplot(detect_filt, aes(x = julianday, y = deploy_lat, color = detectedby)) + 
  geom_point() + facet_wrap(~yearcollected)

# Calculate detection events
# This filters detections down into "events" based on the period of time between detections
# Detectons are collated into one event until the animal is not detected for time_sep seconds
# In this case, we are arbitrarily using one day, larger values = fewer events and vice versa
events = detection_events(detect_filt, time_sep = 86400)

# Species name, which we will use for modelling
spp_to_model = unique(detect_data$scientificname)
spp_to_model

########## Return to presentation ##########

