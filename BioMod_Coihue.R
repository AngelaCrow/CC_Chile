
library("biomod2", quietly = TRUE)
library("sp", quietly = TRUE)
library("raster", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("rgdal", quietly = TRUE)



# Variables ambientales
# Hacer stack de los raster de las variables ambientales segun spp y extent
stfiles=list.files("F:\\COBERTURAS\\Worldclim\\WC_bios_asc", pattern = "*.asc$", full.name=T)
bioclima=stack(stfiles)

###
# M
shapePath <- 'C:/Proyectos/Carlos/CC_Chile'
shapeLayer <- "ecoregionsOI2"
M <- rgdal::readOGR(shapePath, shapeLayer)

# Mask rasters with M
VariablesCrop <- raster::crop(bioclima, M)
myExpl <- raster::mask(VariablesCrop, M) #Species variables delimited by M
myExpl <-stack(myExpl)
names(myExpl)
myExpl<-myExpl[[c(6,8,15,16,17)]]
plot(myExpl)

#DataSpecies_test <-read.csv("Sula_test.csv", h=T) #nombre de la megatabla segun sp y N
presencias <- read.csv("Coihue_pres2_5.csv")
ausencias <- read.csv("Coihue_aus2_5.csv")
DataSpecies<-rbind(presencias,ausencias)
colnames(DataSpecies)[4]<-"Coihue"

myRespName <- 'Coihue'
myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[c('X','Y')]

## Inicio, ajuste de los datos
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName)
myBiomodData
plot(myBiomodData)
### Parametrizacion
myBiomodOption <-BIOMOD_ModelingOptions() 

### Modelizacion
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS',
             'RF', 'MAXENT'),
  models.options = myBiomodOption, 
  NbRunEval=1,
  DataSplit=70,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"Coihue",sep=""))

myBiomodModelOut

myBiomodModelEval <- get_evaluations(myBiomodModelOut)
write.csv(myBiomodModelEval,"myBiomodModelEval.csv")

### Hacer predicciones sobre el raster
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'Coihue',
  selected.models = 'all',
  binary.meth = "TSS",
  compress = 'xz',
  clamping.mask = FALSE,
  output.format = '.grd')

#currentPred
currentPred <- stack("Coihue/proj_Coihue/proj_Coihue_Coihue.grd")
writeRaster(currentPred, filename="Coihue/.tif", overwrite=T, bylayer=TRUE, suffix='names')
plot(currentPred)

### ensemble_modeling#####
myBiomodEM <- BIOMOD_EnsembleModeling( 
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.8),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

# print summary
myBiomodEM

# get evaluation scores
myBiomodEMEval<-get_evaluations(myBiomodEM)
write.csv(myBiomodEMEval,"myBiomodEMEval.csv")

#projBiomodEnsemble
myBiomodEMcurrent <- BIOMOD_EnsembleForecasting( 
  EM.output = myBiomodEM,
  projection.output = myBiomodProj)

myBiomodEMcurrent
plot(myBiomodEMcurrent)

#myBiomodEFPred
myBiomodEMcurrentPred <- stack("Coihue/proj_Coihue/proj_Coihue_Coihue_ensemble.grd")
writeRaster(myBiomodEMcurrentPred, filename="Coihue/.tif", overwrite=T, bylayer=TRUE, suffix='names')
plot(currentPred)
