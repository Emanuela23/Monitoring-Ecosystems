#Species distribution modelling

#no setwd: all the data are based directly on the library 

install.packages("sdm") #extensible framework for developing species distribution models using individual and community-based approaches, generate ensembles of models, evaluate the models, and predict species potential distributions in space and time

library(sdm)
library(raster) #ecological variables: predictors for species distribution
library(rgdal) #import vector layers like species data

file <- system.file("external/species.shp", package="sdm")
species <- shapefile(file)

plot(species)
species
species$Occurrence

plot(species[species$Occurrence == 1,],col='blue',pch=16) #the quadratic parenthesis puts a condition and in this case "equal" is == and comma is put when the condition is finished

#add to the existing plot the absences
points(species[species$Occurrence == 0,],col='red',pch=16) 

# predictors
path <- system.file("external", package="sdm")

lst <- list.files(path=path,pattern='asc$',full.names = T) # 
#asc is a type of file, an extension, called also ascii
lst

preds <- stack(lst)
plot(preds)

cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)

plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# model
d <- sdmData(train=species, predictors=preds)

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm")

p1 <- predict(m1, newdata=preds)

plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)

s1 <- stack(preds,p1)
plot(s1, col=cl)




