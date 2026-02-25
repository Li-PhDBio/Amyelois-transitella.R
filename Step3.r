library(biomod2)
library(terra)

species_data <- read_csv(file = "AT_cleaned.csv")
head(species_data)

# Select the name of the studied species
myRespName <- 'Amyelois transitella'

# Get corresponding presence/absence data
myResp <- species_data$`Amyelois transitella`

# Get corresponding XY coordinates
myRespXY <- species_data[, c('X_WGS84', 'Y_WGS84')]

climdata_current <- list.files("bioclim_back", pattern = ".tif", full.names = TRUE)
climate_data <- terra::rast(climdata_current)

climdata_world <- list.files("bioclim_world", pattern = ".tif", full.names = TRUE)
climate_world <- terra::rast(climdata_world)

# # Select multiple sets of pseudo-absences#
# # Transform true absences into potential pseudo-absences
myResp.PA <- ifelse(myResp == 1, 1, NA)
# myResp.PA.vect <- vect(cbind(myRespXY, myResp.PA), geom = c('X_WGS84','Y_WGS84')) 
# # Format Data with pseudo-absences : random method
myBiomodData.PA <- BIOMOD_FormatingData(resp.name = myRespName,
                                        resp.var = myResp.PA,
                                        expl.var = climate_data,
                                        resp.xy = myRespXY,
                                        data.type = "binary",
                                        PA.nb.rep = 5,
                                        PA.nb.absences = 1000,
                                        PA.strategy = 'sre',
                                        PA.sre.quant = 0.10,
                                        filter.raster = TRUE)
myBiomodData.PA
summary(myBiomodData.PA)
plot(myBiomodData.PA)


# # k-fold selection
# cv.k <- bm_CrossValidation(bm.format = myBiomodData,
#                            strategy = 'kfold',
#                            nb.rep = 2,
#                            k = 3)
# 
# # stratified selection (geographic)
cv.s <- bm_CrossValidation(bm.format = myBiomodData.PA,
                           strategy = 'strat',
                           k = 4,
                           balance = 'presences',
                           strat = 'y')
cv.s <- cbind('_allData_allRun' = TRUE , cv.s)

# cv.s$`_allData_allRun` <- TRUE
# head(cv.k)
head(cv.s)


data(ModelsTable)
ModelsTable

Models <- c('CTA', 'ANN', 'FDA', 'GAM', 'GBM', 'MARS', 'MAXNET', 'RF', 'SRE', 'XGBOOST')

# # bigboss parameters
# opt.b <- bm_ModelingOptions(data.type = 'binary',
#                           models = allModels,
#                           strategy = 'bigboss')

opt.bf <- bm_ModelingOptions(data.type = 'binary',
                             models = Models,
                             strategy = 'bigboss',
                             bm.format = myBiomodData.PA,
                             calib.lines = cv.s)
# # tuned parameters with formated data
# opt.t <- bm_ModelingOptions(data.type = 'binary',
#                             models = c('SRE', 'XGBOOST'),
#                             strategy = 'tuned',
#                             bm.format = myBiomodData)
# 
# opt.b
opt.bf

# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData.PA,
                                    modeling.id = 'Models',
                                    models = Models,
                                    CV.strategy = 'user.defined',
                                    CV.user.table = cv.s,
                                    OPT.strategy = 'user.defined',
                                    OPT.user = opt.bf,
                                    CV.do.full.models = TRUE,
                                    metric.eval = c('TSS', 'AUCroc', 'KAPPA'),
                                    var.import = 3)
# seed.val = 123)
# nb.cpu = 8)
myBiomodModelOut

# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
get_variables_importance(myBiomodModelOut)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'calibration', group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'calibration', group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1, 4, 8, 10, 13)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1, 4, 8, 10, 13)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[4],
                      fixed.var = 'median',
                      do.bivariate = TRUE)

# Explore models' outliers & residuals
bm_ModelAnalysis(bm.mod = myBiomodModelOut,
                 models.chosen = get_built_models(myBiomodModelOut)[c(1, 4, 8, 10, 13)])


# Model ensemble models

myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmedian', 'EMmean', 'EMwmean',
                                                  'EMca', 'EMci', 'EMcv'),
                                      metric.select = 'TSS',
                                      metric.select.thresh = c(0.60),
                                      metric.eval = c('TSS', 'AUCroc', 'KAPPA'),
                                      var.import = 3,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional')
myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodEM, dataset = 'calibration', group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, dataset = 'calibration', group.by = c('full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[7],
                      fixed.var = 'median',
                      do.bivariate = TRUE)
# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = climate_world,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE,
                                  keep.in.memory = FALSE)
myBiomodProj
plot(myBiomodProj)

# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')

myBiomodEMProj
plot(myBiomodEMProj)

# Load environmental variables extracted from BIOCLIM 
climdata_2140126 <- list.files("bioclim_126n", pattern = ".tif", full.names = TRUE)
climate_2140126 <- terra::rast(climdata_2140126)
climdata_4160126 <- list.files("bioclim_126f", pattern = ".tif", full.names = TRUE)
climate_4160126 <- terra::rast(climdata_4160126)
climdata_2140585 <- list.files("bioclim_585n", pattern = ".tif", full.names = TRUE)
climate_2140585 <- terra::rast(climdata_2140585)
climdata_4160585 <- list.files("bioclim_585f", pattern = ".tif", full.names = TRUE)
climate_4160585 <- terra::rast(climdata_4160585)
# Project onto future conditions
myBiomodProjFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        proj.name = 'Future2140126',
                                        new.env = climate_2140126,
                                        models.chosen = 'all',
                                        metric.binary = 'all',
                                        build.clamping.mask = TRUE)
myBiomodEMProjFuture <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                   bm.proj =myBiomodProjFuture,
                                                   models.chosen = 'all',
                                                   metric.binary = 'all',
                                                   metric.filter = 'all')

myBiomodProjFuture4160126 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        proj.name = 'Future4160126',
                                        new.env = climate_4160126,
                                        models.chosen = 'all',
                                        metric.binary = 'all',
                                        build.clamping.mask = TRUE)
myBiomodEMProjFuture4160126 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                   bm.proj =myBiomodProjFuture4160126,
                                                   models.chosen = 'all',
                                                   metric.binary = 'all',
                                                   metric.filter = 'all')

myBiomodProjFuture2140585 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        proj.name = 'Future2140585',
                                        new.env = climate_2140585,
                                        models.chosen = 'all',
                                        metric.binary = 'all',
                                        build.clamping.mask = TRUE)
myBiomodEMProjFuture2140585 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                   bm.proj =myBiomodProjFuture2140585,
                                                   models.chosen = 'all',
                                                   metric.binary = 'all',
                                                   metric.filter = 'all')

myBiomodProjFuture4160585 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                        proj.name = 'Future4160585',
                                        new.env = climate_4160585,
                                        models.chosen = 'all',
                                        metric.binary = 'all',
                                        build.clamping.mask = TRUE)
myBiomodEMProjFuture4160585 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                   bm.proj =myBiomodProjFuture4160585,
                                                   models.chosen = 'all',
                                                   metric.binary = 'all',
                                                   metric.filter = 'all')




# Get a summary report
BIOMOD_Report(bm.out = myBiomodEM, strategy = 'report')

# Get a pre-filled ODMAP
BIOMOD_Report(bm.out = myBiomodEM, strategy = 'ODMAP')

# Get a code report
BIOMOD_Report(bm.out = myBiomodEM, strategy = 'code')


