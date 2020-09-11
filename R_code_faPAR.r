# R_code_faPAR.r
# How to look at chemical cycling from satellites

# levelplot(copNDVI)
#
# setwd("C:/Downloads/")
#
# faPAR <- raster("~/Downloads/c_gls_FAPAR300-RT0_202004200000_GLOBE_PROBAV_V1.0.1.nc")
#
# plot(faPAR)
#
# faPAR10 <- aggregate(faPAR, fact=10)

# install.packages("raster")

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

########################################################################

# regression model between faPAR and NDVI
erosion <- c(12, 14, 16, 24, 26, 40, 55, 67) #amount in kg of erosion 
hm <- c(30, 100, 150, 200, 260, 340, 460, 600) #amount of heavy metals in ppm

plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals")

model1 <- lm(hm ~ erosion)
summary(model1)
abline(model1)

#######
setwd("C:/lab/")
library(raster)
library(rasterdiv)
faPAR10 <- raster("faPAR10.tif")

plot(faPAR10)
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

faPAR10 


library(sf) # to call st_* functions
random.points <- function(x,n)
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}

plot(faPAR10)
pts <- random.points(faPAR10,1000)
copNDVIp <- extract(copNDVI, pts)
faPAR10p <- extract(faPAR10, pts)

# photosynthesis vs. biomass
model2 <- lm(faPAR10p ~ copNDVIp)
plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")


