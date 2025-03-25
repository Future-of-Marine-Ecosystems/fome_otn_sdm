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
library(ggpubr)


#############################################################################
######################## Format data to BIOMOD2 #############################
#############################################################################

# Data to model
spp_env_df = events

# Prepare data for BIOMOD2 -------------------------------------------------
lonlat <- spp_env_df %>% dplyr::select("mean_longitude", "mean_latitude") # Locations (Longitude, Latitude)
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

# Projections

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
                                             new.env = bio_oracle,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')

# Plot ensemble projection
p4 = plot(myBiomodEMProj, plot.output = 'list')[[3]]

# Plot all projections
ggarrange(p1, p2, p3, p4,
          labels = c("RF", "GLM", "GAM", "Ensemble"),
          ncol = 2, nrow = 2)
