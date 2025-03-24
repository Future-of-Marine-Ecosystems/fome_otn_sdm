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
library(biomod2)


#############################################################################
######################## Format data to BIOMOD2 #############################
#############################################################################

# Data to model
spp_env_df = detect_filt

# Prepare data for BIOMOD2 -------------------------------------------------
lonlat <- spp_env_df %>% dplyr::select("deploy_long", "deploy_lat") # Locations (Longitude, Latitude)
presences <- rep(1, nrow(spp_env_df))


# Set up data for species distribution modelling - detections
det_sdm_data <- BIOMOD_FormatingData(
  resp.name = spp_to_model,
  expl.var = bio_oracle,
  resp.var = presences,
  resp.xy = lonlat, 
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
