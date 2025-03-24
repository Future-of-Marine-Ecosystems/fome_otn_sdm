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

# Map...
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

# Species name
spp_to_model = unique(detect_data$scientificname)
spp_to_model

#############################################################################
########################## Climatology of env data  #########################
#############################################################################

# Read location where you have your raster files 
env_dir <- '/Users/estebansal/Documents/GitHub/fome_otn_sdm/Data/copernicus'
env_files <- list.files(env_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE) #binary files... make sure you're readin in the correct format(.nc .tif .grd)

# Load the raster files with raster package and then transform to terra package 
env_stack <- raster::stack(env_files) #raster
env_stack <- terra::rast(env_files) #terra

# Subset enviormental variables for superficial layers only !!!
names(env_stack)
#"zos" ... Sea surface height
#"thetao_depth=0.49402499"... sea surface temperature
#"so_depth=0.49402499" ... salinity 

zos_grep <- grep("zos", names(env_stack))
thetao_grep <- grep("thetao_depth=0.49402499", names(env_stack))
so_grep <- grep("so_depth=0.49402499", names(env_stack))

# stack superficial layers
zos_stack <- subset(env_stack, c(zos_grep))
thetao_stack <- subset(env_stack,c(thetao_grep))
so_stack <- subset(env_stack,c(so_grep))

# Create a global climatology (average across months per env layer)
global_zos <- mean(zos_stack)
global_thetao <- mean(thetao_stack)
global_so <- mean(so_stack)

global_climatology <- c(global_zos, global_thetao, global_so)
names(global_climatology) <- c("ssh", "sst", "sal") # change names


###Create a regional climatology (average across months per env layer)

# Select spatial range .... 
ext <- c(-72, -57, 41, 50) # spatial extent (xmin, xmax, ymin, ymax)
zos_stack <- crop(zos_stack, ext)
thetao_stack <- crop(thetao_stack, ext)
so_stack <- crop(so_stack, ext)

# Create regional climatology
zos_mean <- mean(zos_stack)
thetao_mean <- mean(thetao_stack)
so_mean <- mean(so_stack)

regional_climatology <- c(zos_mean, thetao_mean, so_mean)
names(regional_climatology) <- c("ssh", "sst", "sal")

regional_climatology # raster stack to be used to format biomod

#############################################################################
########################## Extract environmental data ######################
#############################################################################
coords <- spp_env_df[, 3:4] # select coordinates
env_extract <- terra::extract(climatology, coords) #extract
env_extract <- env_extract[ ,-1] # delete first column
spp_env_df <- cbind(spp_env_df, env_extract) # bind them with data frame 

#############################################################################
######################## Format data to BIOMOD2 #############################
#############################################################################

# Prepare data for BIOMOD2 -------------------------------------------------
Var_extracted <- spp_env_df %>% dplyr::select("ssh", "sst", "sal")#Extracted enviormental data
resp.xy <- spp_env_df %>% dplyr::select("lon", "lat") # Locations (Longitude, Latitude)
MyRespVar <- spp_env_df %>% dplyr::select("response") # Prescence/Absence data (binary)

#Rename object to keep it cleaner
myRespXY <- resp.xy
myrespVar <- MyRespVar
expl.var <- Var_extracted


# Set up data for species distribution modelling - detections
det_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model,
  expl.var = regional_climatology,
  resp.var = MyRespVar,
  resp.xy = myRespXY, 
  PA.nb.rep = 1,
  PA.nb.absences = 1000, 
  PA.strategy = "random")


#############################################################################
########################### Modelling in Biomod2  ###########################
#############################################################################

# Modelling in Biomod 
myBiomodModelOut <- BIOMOD_Modeling(bm.format = det_sdm_data,
                                    modeling.id = 'Example',
                                    models = c('RF', 'GLM', 'GAM'),
                                    CV.strategy = 'kfold',
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


#############################################################################
############################# Ensemble modelling ############################
#############################################################################

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

#############################################################################
############################# Ensemble modelling ############################
#############################################################################


# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'Current_proj',
                                             new.env = regional_climatology,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')

myBiomodEMProj
plot(myBiomodEMProj[1])
