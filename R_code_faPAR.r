# R_code_faPAR.r
# How to look at chemical cycling from satellites

# t(copNDVI)
#
# setwd("C:/Downloads/")
#
#..?
#

setwd("C:/lab/")
library(raster)
library(rasterVis)
library(rasterdiv)

plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # removing the data from 253 to 255 (water)
levelplot(copNDVI)

faPAR10 <- raster("C:/lab/faPAR10.tif")
levelplot(faPAR10)

# save the images as PDF
pdf("copNDVI.pdf")
levelplot(copNDVI)
dev.off()

pdf("faPAR10.pdf")
levelplot(faPAR10)
dev.off()

################################

setwd("C:/lab/")

load("faPAR.RData")
ls()

library(raster)
library(rasterdiv)
library(rasterVis)
faPAR10

# the original faPAR from Copernicus is 2 GB
# let's see how much space is needed for an 8-bit set 

writeRaster(copNDVI, "copNDVI.tif") # write the file in our computer
# 5.28 MB

# faPAR levelplot this set # you need rasterVis
faPAR10 <- raster("C:/lab/faPAR10.tif")
levelplot(faPAR10) 
dev.off




