# "R for movement ecologists: Species Distribution Models
# Reid Steele & Esteban Salazar
# BONUS FILE: Create a global SDM using species occurrence databases
# This file DOES NOT require you to run any of the other files in the repository first
# April 9, 2025

# Libraries
library(tidyverse)
library(data.table)
library(glatos)
library(terra)
library(tidyterra)
library(geodata)
library(biomod2)
library(ggpubr)
library(robis)
library(biooracler)


#############################################################################
#################### Download occurrence data from OBIS #####################
#############################################################################

# Our species for this test:
spp_to_model = "Prionace glauca"

# Download blue shark occurrences from OBIS (global)
bshark = occurrence(scientificname = spp_to_model)

# Plot occurrences
world = world(path = getwd())

# Map
bshark_map = ggplot() +
  geom_spatvector(data = world) + # Map baselayer
  coord_sf() + # lat and lon limits
  geom_point(data = bshark, aes(x = decimalLongitude, y = decimalLatitude))
bshark_map



#############################################################################
############ Download global environmental data from Bio-ORACLE #############
#############################################################################



# Download directory
dir = './Data/bio-oracle_g/'

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
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z') # Must have two time coordinates, can be identical
                          )), # No Space Coordinates
  
  # Download salinity
  list(dataset_id = "so_baseline_2000_2019_depthsurf", # Bio-ORACLE dataset name
       variables = "so_mean", # Variable name(s)
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z') # Must have two time coordinates, can be identical
                          )), # No Space Coordinates
  
  # Download chlorophyll
  list(dataset_id = "chl_baseline_2000_2018_depthsurf", # Bio-ORACLE dataset name
       variables = c("chl_mean"), # Variable name(s)
       constraints = list(time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z') # Must have two time coordinates, can be identical
                          )), # No Space Coordinates
  
  # Download bathymetry
  list(dataset_id = "terrain_characteristics", # Bio-ORACLE dataset name
       variables = c("bathymetry_mean"), # Variable name(s)
       constraints = list(time = c('1970-01-01T00:00:00Z', '1970-01-01T00:00:00Z') # Must have two time coordinates, can be identical
                          )) # No Space Coordinates
  
) # End datasets list

# Download rasters
# Since the dataset is now global, everything is going to take quite a bit longer
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




#############################################################################
################# 3. Reformat data to biomod2 Data Format ###################
#############################################################################

# Data to model
spp_env_df = bshark

# Prepare data for BIOMOD2 -------------------------------------------------
lonlat <- spp_env_df %>% dplyr::select("decimalLongitude", "decimalLatitude") # Locations (Longitude, Latitude)
presences <- rep(1, nrow(spp_env_df))


# Set up data for species distribution modelling 
# Note that now that the data is global, this is going to take a while
det_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model,
  expl.var = bio_oracle,
  resp.var = presences,
  resp.xy = lonlat, 
  PA.nb.rep = 1,
  PA.nb.absences = 1000, # Number of Pseudoabsences
  # Since we have more presence points we could crank up the pseudoabsences, but it will take longer
  # Try playing around with the pseudoabsences value above to see what happens if you're interested
  # NOTE: More pseudoabsences means more run time
  PA.strategy = "random")


#############################################################################
################# 4. Running Species Distribution Models  ###################
#############################################################################

# Modelling in Biomod 
# Note that now that the data is global, this is going to take a while
myBiomodModelOut <- BIOMOD_Modeling(bm.format = det_sdm_data,
                                    modeling.id = 'Example',
                                    models = c('RF', 'GLM', 'GAM'),
                                    CV.strategy = 'kfold',
                                    CV.k = 5, 
                                    OPT.strategy = 'default',
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 3,
                                    seed.val = 42, 
                                    do.progress = TRUE,
                                    nb.cpu = 3)


# Plot projections
# Note that now that the data is global, this is going to take a while

# Random forest
bm_projection_rf <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, new.env = bio_oracle,
  models.chosen = get_built_models(myBiomodModelOut)[1],
  proj.name = "Current", metric.binary = "all",
  metric.filter = "all"
)
p1 = plot(bm_projection_rf)

# GLM
bm_projection_glm <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, new.env = bio_oracle,
  models.chosen = get_built_models(myBiomodModelOut)[2],
  proj.name = "Current", metric.binary = "all",
  metric.filter = "all"
)
p2 = plot(bm_projection_glm)

# GAM
bm_projection_gam <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, new.env = bio_oracle,
  models.chosen = get_built_models(myBiomodModelOut)[3],
  proj.name = "Current", metric.binary = "all",
  metric.filter = "all"
)
p3 = plot(bm_projection_gam)



#############################################################################
################ 5. Evaluating Species Distribution Models  #################
#############################################################################

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


########## Return to presentation ##########


#############################################################################
########################### 6. Ensemble modelling ###########################
#############################################################################

# Note that now that the data is global, this is going to take a while
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


# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'Current_proj',
                                             new.env = bio_oracle,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')

# Plot ensemble projection
p4 = plot(myBiomodEMProj, plot.output = 'list')[[1]] # Unweighted
p5 = plot(myBiomodEMProj, plot.output = 'list')[[3]] # Weighted

# Plot all projections
ggarrange(p4, p5,
          labels = c("Unweighted", 'Weighted'),
          ncol = 1, nrow = 2)

# Plot all projections
ggarrange(p1, p2, p3, p4,
          labels = c("RF", "GLM", "GAM", "Ensemble"),
          ncol = 2, nrow = 2)

########## Return to presentation ########## 
