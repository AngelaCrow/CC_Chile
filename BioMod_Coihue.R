
library("biomod2", quietly = TRUE)
library("sp", quietly = TRUE)
library("raster", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("rgdal", quietly = TRUE)



# Variables ambientales
# Hacer stack de los raster de las variables ambientales segun spp y extent
stfiles=list.files("F:/COBERTURAS/CHELSA", pattern = "*.tif$", full.name=T)
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
             'RF'),
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

currentPred_bin <- stack("Coihue/proj_Coihue/proj_Coihue_Coihue_TSSbin.grd")
plot(currentPred_bin)
sumaPres<-sum(currentPred_bin)
plot(sumaPres)
r <-function(x){
  x[x<=1]<-0; x[x>=2]<-1 
  return(x)}
Pres_bin<-calc(sumaPres, fun= r )
plot(Pres_bin)
writeRaster(Pres_bin, filename="Coihue/presente.tif", overwrite=T, suffix='names')
writeRaster(currentPred_bin, filename="Coihue/TSSbin.tif", overwrite=T, bylayer=TRUE, suffix='names')
plot(currentPred_bin)
