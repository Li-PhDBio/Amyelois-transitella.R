---
title: "Amyelois transitella"
output:
  html_document: 
    fig_width: 28
    fig_height: 12.8
    fig_caption: true
    theme: journal
    df_print: default
---

```{r Package}
library(sf)
library(readr)
library(sp)
library(terra)
library(dplyr)
library(caret)
library(dismo)
library(plyr)
library(biomod2)
library(tidyterra)
```

```{r Data.input}
#加载物种分布数据
species_data <- read_csv(
  file = "AT_biomod.csv")
head(species_data)

#Adds another column indicating these are presence points
species_data <- species_data %>%mutate(presence = 1)
head(species_data)


# 筛选数据
myRespName <- 'Amyelois.transitella' # 物种名称
myResp <- species_data$presence # 物种分布数据
myRespXY <- species_data[, c("x", "y")]

# 读取气候数据
climdata_current <- list.files("model_bio", pattern = ".tif", full.names = TRUE)
climate_data <- terra::rast(climdata_current)
```

```{r Absences.point}
# 生成多组假存在数据
myResp.PA <- ifelse(myResp == 1, 1, NA)

# 多组假存在数据
myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
                                     resp.var = myResp.PA,
                                     resp.xy = myRespXY,
                                     expl.var = climate_data,
                                     filter.raster = TRUE,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = c(275, 550),
                                     PA.strategy = 'random')
myBiomodData
#summary(myBiomodData)
plot(myBiomodData)
```

```{r Modle.option}
# k-fold selection
cv.k <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = 'kfold',
                           nb.rep = 2,
                           k = 3)
head(cv.k)
apply(cv.k, 2, table)
plot(cv.k)

bm_MakeFormula(resp.name = myRespName,
               expl.var = head(climate_data),
               type = 'quadratic',
               interaction.level = 0)

# default parameters
opt.d <- bm_ModelingOptions(data.type = 'binary',
                            models = c("ANN", "CTA","FDA","GAM","GBM","GLM","MARS","MAXNET","RF","RFd","SRE","XGBOOST"),
                            strategy = 'default')

opt.d
opt.d@models
opt.d@options$ANN.binary.nnet.nnet
names(opt.d@options$ANN.binary.nnet.nnet@args.values)
```

```{r Single.modle, warning=FALSE, paged.print=TRUE}
# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    models = c("ANN","CTA","FDA","GAM","GBM","GLM","MARS","MAXNET","RF","RFd","SRE","XGBOOST"),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    metric.eval = c('TSS','ROC','KAPPA'),
                                    var.import = 2,
                                    seed.val = 42)
myBiomodModelOut
```

```{r Modles.corce, warning=FALSE, paged.print=TRUE}
# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
get_variables_importance(myBiomodModelOut)

# Represent evaluation scores & variables importance
# 显示所有三个指标
bm_PlotEvalMean(bm.out = myBiomodModelOut, 
                dataset = 'calibration',
                metric.eval = c('TSS', 'ROC'))
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation',metric.eval = c('TSS', 'ROC'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'median') 
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'min')
```

```{r ensemblemodels, warning=FALSE, paged.print=TRUE}
# Model ensemble models
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by = 'all',
  em.algo = c('EMwmean', 'EMmean'),      # 两种集成方法
  metric.select = c('TSS'),              # TSS更稳定，避免过拟合
  metric.select.thresh = c(0.80),        # 保留验证TSS > 0.90的算法
  metric.eval = c('TSS', 'ROC', 'KAPPA'),
  var.import = 5,
  EMci.alpha = 0.05,
  EMwmean.decay = 'proportional'
)
myBiomodEM
```

```{r EMmodels.scores, warning=FALSE}
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'full.name',metric.eval = c('TSS', 'ROC'))
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('full.name', 'full.name'))
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
```


```{r Project.single.models, warning=FALSE, paged.print=TRUE}
# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = climate_data,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj)
```

```{r Project.ensemble.models, warning=FALSE, paged.print=TRUE}
# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'CurrentEM',
                                             new.env = climate_data,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
myBiomodEMProj
plot(myBiomodEMProj)
```

```{r 2140126, warning=FALSE, paged.print=TRUE}
# 2140126
# Load environmental variables extracted from BIOCLIM 2
ssp2140126 <- list.files("2140126", pattern = ".tif$", full.names = TRUE)
future2140126 <- stack(ssp2140126)
# Project onto future conditions
BiomodProj2140126 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                       proj.name = '2140126',
                                       new.env = future2140126,
                                       models.chosen = 'all',
                                       metric.binary = 'all',
                                       build.clamping.mask = TRUE)

# Load current and future binary projections
CurrentProj <- get_predictions(myBiomodProj, metric.binary = "TSS")
FutureProj2140126 <- get_predictions(BiomodProj2140126, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize2140126 <- BIOMOD_RangeSize(proj.current = CurrentProj, 
                                             proj.future = FutureProj2140126)

myBiomodRangeSize2140126$Compt.By.Models
plot(myBiomodRangeSize2140126$Diff.By.Pixel)

FutureProj2140126EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                  bm.proj = BiomodProj2140126,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')

# Project ensemble models (building single projections)
FutureProj2140126EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                  proj.name = '2140126EM',
                                                  new.env = future2140126,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')
FutureProj2140126EM
plot(FutureProj2140126EM)


# Represent main results 
gg_2140126 = bm_PlotRangeSize(bm.range = myBiomodRangeSize2140126,
                              do.count = TRUE,
                              do.perc = TRUE,
                              do.maps = TRUE,
                              do.mean = TRUE,
                              do.plot = TRUE)
```

```{r ssp5852140, fig.height=10, fig.width=12, warning=FALSE}
# 2140585
# Load environmental variables extracted from BIOCLIM 2
ssp2140585 <- list.files("2140585", pattern = ".tif$", full.names = TRUE)
future2140585 <- terra::rast(ssp2140585)
# Project onto future conditions
BiomodProj2140585 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                       proj.name = '2140585',
                                       new.env = future2140585,
                                       models.chosen = 'all',
                                       metric.binary = 'TSS',
                                       build.clamping.mask = TRUE)

# Load current and future binary projections
CurrentProj <- get_predictions(myBiomodProj, metric.binary = "TSS")
FutureProj2140585 <- get_predictions(BiomodProj2140585, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize2140585 <- BIOMOD_RangeSize(proj.current = CurrentProj, 
                                             proj.future = FutureProj2140585)

myBiomodRangeSize2140585$Compt.By.Models
plot(myBiomodRangeSize2140585$Diff.By.Pixel)



FutureProj2140585EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                  bm.proj = BiomodProj2140585,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')

# Project ensemble models (building single projections)
FutureProj2140585EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                  proj.name = '2140585EM',
                                                  new.env = future2140585,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')
FutureProj2140585EM
plot(FutureProj2140585EM)

# Represent main results 
gg_2140585 = bm_PlotRangeSize(bm.range = myBiomodRangeSize2140585,
                              do.count = TRUE,
                              do.perc = TRUE,
                              do.maps = TRUE,
                              do.mean = TRUE,
                              do.plot = TRUE)
```

```{r ssp1262140, fig.height=10, fig.width=12, warning=FALSE, paged.print=TRUE}
# 4160126
# Load environmental variables extracted from BIOCLIM 2
ssp4160126 <- list.files("4160126", pattern = ".tif$", full.names = TRUE)
future4160126 <- terra::rast(ssp4160126)
# Project onto future conditions
BiomodProj4160126 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                       proj.name = '4160126',
                                       new.env = future4160126,
                                       models.chosen = 'all',
                                       metric.binary = 'TSS',
                                       build.clamping.mask = TRUE)

# Load current and future binary projections
CurrentProj <- get_predictions(myBiomodProj, metric.binary = "TSS")
FutureProj4160126 <- get_predictions(BiomodProj4160126, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize4160126 <- BIOMOD_RangeSize(proj.current = CurrentProj, 
                                             proj.future = FutureProj4160126)

myBiomodRangeSize4160126$Compt.By.Models
plot(myBiomodRangeSize4160126$Diff.By.Pixel)



FutureProj4160126EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                  bm.proj = BiomodProj4160126,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')

# Project ensemble models (building single projections)
FutureProj4160126EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                  proj.name = '4160126EM',
                                                  new.env = future4160126,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')
FutureProj4160126EM
plot(FutureProj4160126EM)


# Represent main results 
gg_4160126 = bm_PlotRangeSize(bm.range = myBiomodRangeSize4160126,
                              do.count = TRUE,
                              do.perc = TRUE,
                              do.maps = TRUE,
                              do.mean = TRUE,
                              do.plot = TRUE)
```

```{r ssp5854160, fig.height=10, fig.width=12, warning=FALSE, paged.print=TRUE}
# 4160585
# Load environmental variables extracted from BIOCLIM 2
ssp4160585 <- list.files("4160585", pattern = ".tif$", full.names = TRUE)
future4160585 <- terra::rast(ssp4160585)
# Project onto future conditions
BiomodProj4160585 <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                       proj.name = '4160585',
                                       new.env = future4160585,
                                       models.chosen = 'all',
                                       metric.binary = 'TSS',
                                       build.clamping.mask = TRUE)

# Load current and future binary projections
CurrentProj <- get_predictions(myBiomodProj, metric.binary = "TSS")
FutureProj4160585 <- get_predictions(BiomodProj4160585, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize4160585 <- BIOMOD_RangeSize(proj.current = CurrentProj, 
                                             proj.future = FutureProj4160585)

myBiomodRangeSize4160585$Compt.By.Models
plot(myBiomodRangeSize4160585$Diff.By.Pixel)



FutureProj4160585EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                  bm.proj = BiomodProj4160585,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')

# Project ensemble models (building single projections)
FutureProj4160585EM <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                  proj.name = '4160585EM',
                                                  new.env = future4160585,
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all')
FutureProj4160585EM
plot(FutureProj4160585EM)


# Represent main results 
gg_4160585 = bm_PlotRangeSize(bm.range = myBiomodRangeSize4160585,
                              do.count = TRUE,
                              do.perc = TRUE,
                              do.maps = TRUE,
                              do.mean = TRUE,
                              do.plot = TRUE)
```

