# R_code_snow.r

setwd("C:/lab/")
library(raster)
install.packages("ncdf4") # to import ncdf files (most of the Copernicus data are based on this raster data format)
library(ncdf4)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)

# Exercise: plot the snow cover with the color ramp palette
plot(snowmay, col=cl)

# create a folder called "snow" to import all the "snow" data contained in it
# Slow manner to import the set
setwd("C:/lab/snow/")

snow2000 <- raster("snow2000r.tif")
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))

plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

############################

# fast version of import and plot many data for lazy people
# make the list of the files we are going to import with the list.files function 
# pattern: the files have the same name until "snow20" -> pattern means the common letters
rlist <- list.files(pattern="snow")
rlist

#lapply function imports all of the raster
import <- lapply(rlist, raster)

# we take different layers and put it all together in an object called stack
snow.multitemp <- stack(import)
snow.multitemp
plot(snow.multitemp, col=cl)

# let's make a prediction 
# run the scripts with the function "source"
source("prediction.r")

plot(predicted.snow.2025.norm, col=cl)

#################################
# day 2nd

setwd("C:/lab/snow/")

# Exercise: import the snow cover images altogether
library(raster)

rlist <- list.files(pattern="snow")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl)

load("snow.RData")

prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# export the output
# you made the calculation and you want to send the output to a colleague 
# writeRaster function -> writes the entire Raster object to a file 
writeRaster(prediction, "final.tif")

# final stack (let's make a graph with the data altogether)
final.stack <- stack(snow.multitemp, prediction)
plot(final.stack, col=cl)

# export the R graph 
pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

png("my_final_exciting_graph.png")
plot(final.stack, col=cl)
dev.off()




