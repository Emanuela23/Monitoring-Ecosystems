# R_code_faPAR.r
# How to look at chemical cycling from satellites

# levelplot(copNDVI)
#
# setwd("C:/Downloads/")
#
# faPAR <- raster("C:/Downloads/c_gls_FAPAR300-RT0_202004200000_GLOBE_PROBAV_V1.0.1.nc")
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

################################ day 2

setwd("C:/lab/")

load("faPAR.RData")
ls() #When invoked with no argument at the top level prompt, ls shows what data sets and functions a user has defined. When invoked with no argument inside a function, ls returns the names of the functions local variables.

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

model1 <- lm(hm ~ erosion) #lm is used to fit linear models
summary(model1)
abline(model1) #adds one or more straight lines through the current plot

############ faPAR vs. NDVI model 
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
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') #as manages the relations that allow coercing an object to a given class 
# st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random') #spsample sample point locations within a square area, a grid, a polygon, or on a spatial line
}

plot(faPAR10)
pts <- random.points(faPAR10,1000)
copNDVIp <- extract(copNDVI, pts) #Extract values from a Raster* object at the locations of other spatial data.
faPAR10p <- extract(faPAR10, pts)

# photosynthesis vs. biomass
model2 <- lm(faPAR10p ~ copNDVIp)
summary(model2)
plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")


