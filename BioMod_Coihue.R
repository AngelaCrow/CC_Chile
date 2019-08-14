
library("biomod2")
library("sp")
library("raster")
library("dplyr")
library("rgdal")
library("fuzzySim")

setwd("C:/Proyectos/Carlos/CC_Chile/")
outputFolder<-"C:/Proyectos/Carlos/CC_Chile/"

# Hacer stack de los raster de las variables ambientales segun spp y extent
bioclima<-stack(list.files(path="F:/COBERTURAS/CHELSA/Global/Bioclimas", pattern = "*.tif$", full.names=TRUE))

# Regionalization shapefile folder
polygonsOfInterest <- shapefile("C:/Proyectos/Carlos/CC_Chile/Presente_eco.shp")
#M<-raster("C:/Proyectos/Carlos/CC_Chile/mask/mask.asc") Export_Output

#Adelgace los puntos usando una malla de 0.041666
presencias <- read.csv("Coihue_pres2_5.csv")
ausencias <- read.csv("Coihue_aus2_5.csv")
occsData<-rbind(presencias,ausencias)
colnames(occsData)[4]<-"Coihue"

crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sp::coordinates(occsData) <- c("X","Y")
sp::proj4string(occsData) <- crs.wgs84

# Extract envorimental varibales with species occurrences
covarData <- raster::extract(bioclima, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>%
  dplyr::select(.dots=names(bioclima)) %>%
  complete.cases
covarData <- covarData[completeDataCases, ]

####Variables selection####
speciesCol <- match("Coihue", names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)

select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.txt"))

#Raster covariables selected for model calibration
selectedVariables <- bioclima[[select_var]]
# Variables ambientales
# Cortar bioclimas con mascara M

selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
myExpl <- raster::mask(selectedVariablesCrop,polygonsOfInterest) #Species variables delimited by M
myExpl<-stack(myExpl)
names(myExpl)
plot(myExpl)

dir.create("Presente")
writeRaster(myExpl,file.path("Presente/.tif"),bylayer = T, suffix='names',overwrite = TRUE)

###
# Regionalization shapefile folder
#shapePath <- 'C:/Proyectos/Carlos'
#shapeLayer <- "M_ecost_Clip"
#polygonsOfInterest <- rgdal::readOGR(shapePath, shapeLayer)

## BIOMOD####

#DataFormating####

DataSpecies<-rbind(presencias,ausencias)
colnames(DataSpecies)[4]<-"Coihue"

myRespName <- 'Coihue'
myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[c('X','Y')]
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
  models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS','RF'),
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

plot(myBiomodProj)

#Guardar raster por separado
#currentPred <- stack("Coihue/proj_Coihue/proj_Coihue_Coihue.grd")
#writeRaster(currentPred, filename="Coihue/.tif", overwrite=T, bylayer=TRUE, suffix='names')

### ensemble_modeling#####
myBiomodEM <- BIOMOD_EnsembleModeling( 
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.8),
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional')

# print summary
myBiomodEM

# get evaluation scores
myBiomodEMEval<-get_evaluations(myBiomodEM)
write.csv(myBiomodEMEval,"myBiomodEMEval.csv")

# Creating ensembles projections 
myBiomodEM_proj <-BIOMOD_EnsembleForecasting(EM.output  = myBiomodEM,
                                              projection.output = myBiomodProj,
                                              selected.models = 'all',
                                              proj.name = 'Coihue',
                                              binary.meth = "TSS")

currentPred <- stack("Coihue/proj_Coihue/proj_Coihue_Coihue_ensemble_TSSbin.grd")
writeRaster(currentPred, 
            file.path("Coihue/proj_Coihue/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)
#currentPred <-currentPred[7]
#writeRaster(currentPred, filename="Coihue/.tif", overwrite=T, bylayer=TRUE, suffix='names')

#--------------------------------- 2041-2060
CESM1CAM5_fc_rcp45<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2041-2060/bio/CESM1-CAM5_rcp45", pattern="*.tif$", full.names=TRUE))  
CESM1CAM5_fc_rcp85<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2041-2060/bio/CESM1-CAM5_rcp85", pattern="*.tif$", full.names=TRUE))
GFDLCM3_fc_rcp45<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2041-2060/bio/GFDL-CM3_rcp45", pattern="*.tif$", full.names=TRUE))  
GFDLCM3_fc_rcp85<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2041-2060/bio/GFDL-CM3_rcp85", pattern="*.tif$", full.names=TRUE))  

# Raster covariables selected for model transfer
#seleccion variables CC
env_CESM1CAM5_fc_rcp45<-CESM1CAM5_fc_rcp45[[select_var]]
env_CESM1CAM5_fc_rcp45 <- raster::crop(env_CESM1CAM5_fc_rcp45, polygonsOfInterest)
myExpl_CESM1CAM5_fc_rcp45 <- stack(raster::mask(env_CESM1CAM5_fc_rcp45,  polygonsOfInterest)) 
rm(CESM1CAM5_fc_rcp45,env_CESM1CAM5_fc_rcp45)

env_CESM1CAM5_fc_rcp85<-CESM1CAM5_fc_rcp85[[select_var]]
env_CESM1CAM5_fc_rcp85 <- raster::crop(env_CESM1CAM5_fc_rcp85, polygonsOfInterest)
myExpl_CESM1CAM5_fc_rcp85 <- stack(raster::mask(env_CESM1CAM5_fc_rcp85,  polygonsOfInterest)) 
rm(env_CESM1CAM5_fc_rcp85,CESM1CAM5_fc_rcp85)

env_GFDLCM3_fc_rcp45<-GFDLCM3_fc_rcp45[[select_var]]
env_GFDLCM3_fc_rcp45 <- raster::crop(env_GFDLCM3_fc_rcp45, polygonsOfInterest)
myExpl_GFDLCM3_fc_rcp45 <- stack(raster::mask(env_GFDLCM3_fc_rcp45,  polygonsOfInterest))
rm(env_GFDLCM3_fc_rcp45, GFDLCM3_fc_rcp45)

env_GFDLCM3_fc_rcp85<-GFDLCM3_fc_rcp85[[select_var]]
env_GFDLCM3_fc_rcp85 <- raster::crop(env_GFDLCM3_fc_rcp85, polygonsOfInterest)
myExpl_GFDLCM3_fc_rcp85 <- raster::mask(env_GFDLCM3_fc_rcp85,  polygonsOfInterest) 
rm(env_GFDLCM3_fc_rcp85,GFDLCM3_fc_rcp85)

###
#--------------------------------- 2061-2080
CESM1CAM5_fl_rcp45<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2061-2080/CESM1-CAM5_rcp45", pattern="*.tif$", full.names=TRUE))  
CESM1CAM5_fl_rcp85<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2061-2080/CESM1-CAM5_rcp85", pattern="*.tif$", full.names=TRUE))  
GFDLCM3_fl_rcp45<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2061-2080/GFDL-CM3_rcp45", pattern="*.tif$", full.names=TRUE))  
GFDLCM3_fl_rcp85<-stack(list.files(path = "F:/COBERTURAS/CHELSA/Global/cmip5/2061-2080/GFDL-CM3_rcp85", pattern="*.tif$", full.names=TRUE))  

env_CESM1CAM5_fl_rcp45 <-CESM1CAM5_fl_rcp45[[select_var]]
env_CESM1CAM5_fl_rcp45 <- raster::crop(env_CESM1CAM5_fl_rcp45, polygonsOfInterest)
myExpl_CESM1CAM5_fl_rcp45 <- stack(raster::mask(env_CESM1CAM5_fl_rcp45,  polygonsOfInterest)) 
rm(env_CESM1CAM5_fl_rcp45,CESM1CAM5_fl_rcp45)

env_CESM1CAM5_fl_rcp85<-CESM1CAM5_fl_rcp85[[select_var]]
env_CESM1CAM5_fl_rcp85 <- raster::crop(env_CESM1CAM5_fl_rcp85, polygonsOfInterest)
myExpl_CESM1CAM5_fl_rcp85 <- stack(raster::mask(env_CESM1CAM5_fl_rcp85,  polygonsOfInterest)) 
rm(env_CESM1CAM5_fl_rcp85,CESM1CAM5_fl_rcp85)

env_GFDLCM3_fl_rcp45<-GFDLCM3_fl_rcp45[[select_var]]
env_GFDLCM3_fl_rcp45 <- raster::crop(env_GFDLCM3_fl_rcp45, polygonsOfInterest)
myExpl_GFDLCM3_fl_rcp45 <- stack(raster::mask(env_GFDLCM3_fl_rcp45,  polygonsOfInterest)) 
rm(env_GFDLCM3_fl_rcp45,GFDLCM3_fl_rcp45)

env_GFDLCM3_fl_rcp85<-GFDLCM3_fl_rcp85[[select_var]]
env_GFDLCM3_fl_rcp85 <- raster::crop(env_GFDLCM3_fl_rcp85, polygonsOfInterest)
myExpl_GFDLCM3_fl_rcp85 <- stack(raster::mask(env_GFDLCM3_fl_rcp85,  polygonsOfInterest))
rm(env_GFDLCM3_fl_rcp85,GFDLCM3_fl_rcp85)

###################################################
### TRANSFERENCIA AL FUTURO
###################################################
###################################################

###CESM1CAM5_fc_rcp45####
### Hacer predicciones sobre el raster futuros
myBiomodEM_CESM1CAM5_fc_rcp45<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                           new.env = myExpl_CESM1CAM5_fc_rcp45,
                                                           selected.models = 'all',
                                                           proj.name = "CESM1CAM5_fc_rcp45",
                                                           binary.meth = "TSS"
)

myBiomodEM_CESM1CAM5_fc_rcp45
plot(myBiomodEM_CESM1CAM5_fc_rcp45)

FutureProj <- stack("Coihue/proj_CESM1CAM5_fc_rcp45/proj_CESM1CAM5_fc_rcp45_Coihue_ensemble_TSSbin.grd")

writeRaster(FutureProj, 
            file.path("Coihue/proj_CESM1CAM5_fc_rcp45/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="CESM1CAM5_fc_rcp45_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_CESM1CAM5_fc_rcp45/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_CESM1CAM5_fc_rcp45/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###CESM1CAM5_fc_rcp85####
### Hacer predicciones sobre el raster futuros
myBiomodEM_CESM1CAM5_fc_rcp85<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                           new.env = myExpl_CESM1CAM5_fc_rcp85,
                                                           selected.models = 'all',
                                                           proj.name = "CESM1CAM5_fc_rcp85",
                                                           binary.meth = "TSS"
)

myBiomodEM_CESM1CAM5_fc_rcp85
plot(myBiomodEM_CESM1CAM5_fc_rcp85)

FutureProj <- stack("Coihue/proj_CESM1CAM5_fc_rcp85/proj_CESM1CAM5_fc_rcp85_Coihue_ensemble_TSSbin.grd")
writeRaster(FutureProj, 
            file.path("Coihue/proj_CESM1CAM5_fc_rcp85/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="CESM1CAM5_fc_rcp85_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_CESM1CAM5_fc_rcp85/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_CESM1CAM5_fc_rcp85/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###CESM1CAM5_fl_rcp45####
### Hacer predicciones sobre el raster futuros
myBiomodEM_CESM1CAM5_fl_rcp45<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                           new.env = myExpl_CESM1CAM5_fl_rcp45,
                                                           selected.models = 'all',
                                                           proj.name = "CESM1CAM5_fl_rcp45",
                                                           binary.meth = "TSS"
)

myBiomodEM_CESM1CAM5_fl_rcp45
plot(myBiomodEM_CESM1CAM5_fl_rcp45)

#myBiomodEFPred
FutureProj <- stack("Coihue/proj_CESM1CAM5_fl_rcp45/proj_CESM1CAM5_fl_rcp45_Coihue_ensemble_TSSbin.grd")

writeRaster(FutureProj, 
            file.path("Coihue/proj_CESM1CAM5_fl_rcp45/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)


#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="CESM1CAM5_fl_rcp45_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_CESM1CAM5_fl_rcp45/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_CESM1CAM5_fl_rcp45/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###CESM1CAM5_fl_rcp85####
### Hacer predicciones sobre el raster futuros
myBiomodEM_CESM1CAM5_fl_rcp85<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                           new.env = myExpl_CESM1CAM5_fl_rcp85,
                                                           selected.models = 'all',
                                                           proj.name = "CESM1CAM5_fl_rcp85",
                                                           binary.meth = "TSS"
)

myBiomodEM_CESM1CAM5_fl_rcp85
plot(myBiomodEM_CESM1CAM5_fl_rcp85)

#myBiomodEFPred
FutureProj <- stack("Coihue/proj_CESM1CAM5_fl_rcp85/proj_CESM1CAM5_fl_rcp85_Coihue_ensemble_TSSbin.grd")
writeRaster(FutureProj, 
            file.path("Coihue/proj_CESM1CAM5_fl_rcp85/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="CESM1CAM5_fl_rcp85_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_CESM1CAM5_fl_rcp85/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_CESM1CAM5_fl_rcp85/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###GFDLCM3_fc_rcp45####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GFDLCM3_fc_rcp45<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                         new.env = myExpl_GFDLCM3_fc_rcp45,
                                                         selected.models = 'all',
                                                         proj.name = "GFDLCM3_fc_rcp45",
                                                         binary.meth = "TSS"
)

myBiomodEM_GFDLCM3_fc_rcp45
plot(myBiomodEM_GFDLCM3_fc_rcp45)

#myBiomodEFPred
FutureProj <- stack("Coihue/proj_GFDLCM3_fc_rcp45/proj_GFDLCM3_fc_rcp45_Coihue_ensemble_TSSbin.grd")
writeRaster(FutureProj, 
            file.path("Coihue/proj_GFDLCM3_fc_rcp45/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)


#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="GFDLCM3_fc_rcp45_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_GFDLCM3_fc_rcp45/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_GFDLCM3_fc_rcp45/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###GFDLCM3_fc_rcp85####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GFDLCM3_fc_rcp85<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                         new.env = myExpl_GFDLCM3_fc_rcp85,
                                                         selected.models = 'all',
                                                         proj.name = "GFDLCM3_fc_rcp85",
                                                         binary.meth = "TSS"
)

myBiomodEM_GFDLCM3_fc_rcp85
plot(myBiomodEM_GFDLCM3_fc_rcp85)

#myBiomodEFPred
FutureProj <- stack("Coihue/proj_GFDLCM3_fc_rcp85/proj_GFDLCM3_fc_rcp85_Coihue_ensemble_TSSbin.grd")
writeRaster(FutureProj, 
            file.path("Coihue/proj_GFDLCM3_fc_rcp85/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="GFDLCM3_fc_rcp85_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_GFDLCM3_fc_rcp85/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_GFDLCM3_fc_rcp85/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###GFDLCM3_fl_rcp45####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GFDLCM3_fl_rcp45<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                         new.env = myExpl_GFDLCM3_fl_rcp45,
                                                         selected.models = 'all',
                                                         proj.name = "GFDLCM3_fl_rcp45",
                                                         binary.meth = "TSS"
)

myBiomodEM_GFDLCM3_fl_rcp45
plot(myBiomodEM_GFDLCM3_fl_rcp45)

#myBiomodEFPred
FutureProj <- stack("Coihue/proj_GFDLCM3_fl_rcp45/proj_GFDLCM3_fl_rcp45_Coihue_ensemble_TSSbin.grd")
writeRaster(FutureProj, 
            file.path("Coihue/proj_GFDLCM3_fl_rcp45/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="GFDLCM3_fl_rcp45_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_GFDLCM3_fl_rcp45/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_GFDLCM3_fl_rcp45/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)

###GFDLCM3_fl_rcp85####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GFDLCM3_fl_rcp85<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                         new.env = myExpl_GFDLCM3_fl_rcp85,
                                                         selected.models = 'all',
                                                         proj.name = "GFDLCM3_fl_rcp85",
                                                         binary.meth = "TSS"
)

myBiomodEM_GFDLCM3_fl_rcp85
plot(myBiomodEM_GFDLCM3_fl_rcp85)

#myBiomodEFPred
FutureProj <- stack("Coihue/proj_GFDLCM3_fl_rcp85/proj_GFDLCM3_fl_rcp85_Coihue_ensemble_TSSbin.grd")
writeRaster(FutureProj, 
            file.path("Coihue/proj_GFDLCM3_fl_rcp85/.tif"), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

#RangeSize
myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="GFDLCM3_fl_rcp85_rs")

write.csv(myBiomodRangeSize$Compt.By.Models, "Coihue/proj_GFDLCM3_fl_rcp85/rangesize.csv")
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path("Coihue/proj_GFDLCM3_fl_rcp85/.tif"), 
            suffix='names',
            overwrite= TRUE)

rm(FutureProj, myBiomodRangeSize)
