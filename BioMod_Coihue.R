
library(biomod2)
library(sp)
library(raster)

DataSpecies

#DataSpecies_test <-read.csv("Sula_test.csv", h=T) #nombre de la megatabla segun sp y N
myRespName <- 'Coihue'

# Variables ambientales
# Hacer stack de los raster de las variables ambientales segun spp y extent
stfiles=list.files("~/Google Drive/Proyectos/Invasoras_Conabio/ISLAS/Presente", pattern = "*.asc$", full.name=T)
myExpl=stack(stfiles)
myExpl
names(myExpl)
plot(myExpl)

myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[c('long','lat')]

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
currentPred <- stack("/proj_Coihue/proj_Coihue.grd")
writeRaster(currentPred, filename="Coihue/.asc", overwrite=T, bylayer=TRUE, suffix='names')
plot(currentPred)

currentPred_bin <- stack("/proj_Coihue/proj_Coihue_TSSbin.grd")
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
