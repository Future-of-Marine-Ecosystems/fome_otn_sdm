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
library(plotly)
library(biomod2)

# source functions
source('functions.R')

# Data directory
data_dir = './Data/'

# Load in data
detect_data = otn_read(data_dir)

# How many animals
length(unique(detect_data$animal_id))

# Number of detections by receiver group
dets_rg = group_by(detect_data, detectedby) %>% summarize(n = n())
dets_rg

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

# Calculate detection events
events = detection_events(detect_filt, time_sep = 86400)

# test = otn_events(detect_filt)

# Species name
spp_to_model = unique(detect_data$scientificname)
spp_to_model

############################################################################################
############################################################################################
############################################################################################

# Download environmental data

# Set up data for species distribution modelling - detections
det_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model, resp.var = rep(1, nrow(detect_filt)), expl.var = envcrop_selvar,
  resp.xy = detect_filt[, c("decimalLongitude", "decimalLatitude")],
  PA.nb.rep = 1, PA.nb.absences = 1000, PA.strategy = "random", filter.raster = TRUE)

# Set up data for species distribution modelling - events
evt_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model, resp.var = rep(1, nrow(events)), expl.var = envcrop_selvar,
  resp.xy = events[, c("deploy_long", "deploy_lat")],
  PA.nb.rep = 1, PA.nb.absences = 1000, PA.strategy = "random", filter.raster = TRUE)








