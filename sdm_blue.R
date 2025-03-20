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
#############################################################################
######################## Create Pseudo-absences (PA) ########################
#############################################################################

library(dismo)

range(detect_filt$deploy_long)
range(detect_filt$deploy_lat)
date_sequence<- range(detect_filt$detection_timestamp_utc)

ext <- extent(-70.05312, -58.81101,42.83728, 48.85370)
n= 1000
p = spp_env_df %>% dplyr::select("deploy_long", "deploy_lat") 
mask <- raster('/Users/estebansal/Downloads/cmems_mod_glo_phy_anfc_0.083deg-sst-anomaly_P1D-m_1742477208324.nc')
pa <- randomPoints(mask, n, p, ext=NULL, extf=1.1, excludep=TRUE)
pa <- as.data.frame(pa)

#Create random dates for PA
random_dates <- sample(date_sequence, size = 1000, replace = TRUE)
pa <- cbind(pa, random_dates)
pa$Response <- 0
names(pa) <- c("lon", "lat", "date", "response")
pa$id <- "PA"

pa <- pa%>% dplyr::select("id", "date", "lon", "lat", "response")


# Format True Prescences
spp_env_df <- detect_filt
spp_env_df <- spp_env_df %>% dplyr::select(2, 19, 21,22, 37,38)
names(spp_env_df) <- c('id', 'date', 'lon', 'lat', 'sst', 'response')
spp_env_df <- as.data.frame(spp_env_df)
spp_env_df$Response <- 1


spp_PA <- rbind(spp_env_df, pa)


#############################################################################
########################### Extract enviormental data #######################
#############################################################################
  
library(rerddapXtracto)
spp_PA # Our data frame
# Set spatio-temporal limits
xpos <- spp_PA$lon
ypos <- spp_PA$lat
tpos <- spp_PA$date



# Sea Surface Temperature (SST) 
sst_dataInfo <- rerddap::info('jplMURSST41') #daily, 0.025 degrees, gaps 

ssta <- rxtracto(sst_dataInfo, parameter = 'analysed_sst', xcoord = xpos, ycoord = ypos, tcoord =tpos, xlen = .1, ylen = .1,  progress_bar = TRUE)

# Extract SST mean... 
sst_mean <- ssta$`mean analysed_sst`
summary(sst_mean)



#CHL 
chl_dataInfo <- rerddap::info('erdMH1chla1day_R2022SQ')  




#############################################################################
######################## Format data to BIOMOD2 #############################
#############################################################################

length(sst_mean)
spp_PA$sst <- sst_mean

# Prepare data for BIOMOD2 -------------------------------------------------
Var_extracted <- spp_PA %>% dplyr::select("sst")#Extracted enviormental data
resp.xy <- spp_PA %>% dplyr::select("lon", "lat") # Locations (Longitude, Latitude)
MyRespVar <- spp_PA %>% dplyr::select("response") # Prescence/Absence data (binary)

#Rename object to keep it cleaner
myRespXY <- resp.xy
myrespVar <- MyRespVar
expl.var <- Var_extracted


# Set up data for species distribution modelling - detections
det_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model,
  expl.var = expl.var,
  resp.var = MyRespVar,
  resp.xy = myRespXY)

# Set up data for species distribution modelling - events
evt_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model, 
  resp.var = rep(1, nrow(events)), expl.var = envcrop_selvar,
  resp.xy = events[, c("deploy_long", "deploy_lat")],
  PA.nb.rep = 1, PA.nb.absences = 1000, PA.strategy = "random", filter.raster = TRUE)

# Modelling in Biomod 

myBiomodModelOut <- BIOMOD_Modeling(bm.format = det_sdm_data,
                                    modeling.id = 'Example',
                                    models = c('RF', 'GLM', 'GAM'),
                                    CV.strategy = 'kfold',
                                    #CV.nb.rep = 5,
                                    #CV.perc = 0.8,
                                    CV.k = 5, 
                                    OPT.strategy = 'default',
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 3,
                                    seed.val = 42, 
                                    do.progress = TRUE)


# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodModelOut)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'mean')


# Ensemble model

myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean', 'EMcv', 'EMwmean'),
                                      metric.select = c('ROC'),
                                      metric.select.thresh = c(0.6),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional')
myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

