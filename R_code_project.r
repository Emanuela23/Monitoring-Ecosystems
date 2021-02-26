# 16. R code project 

#Greenland monitoring

#Ice sheet melt extent
setwd("C:/lab/Greenland/")
library(raster)
library(rasterVis)

rlist <- list.files(pattern="melt")
rlist

import <- lapply(rlist,raster)
Gr <- stack(import)
levelplot(Gr)

#Ice sheet difference 2007 - 1979
rast <- raster("1979_melt.tif")
rast2 <- raster("2007_melt.tif")
dif <- rast2 - rast
clb <- colorRampPalette(c("blue","white","red"))(100)
plot(dif,col=clb,main="Melt area variation 1979 - 2007")
dev.off()

#Scatterplot cumulative annual melt area
setwd("C:/lab/Greenland/")
library(ggplot2)
library(readxl)
data <- read_excel("greenland_annual.xlsx")
head(data)
attach(data)
scatter.smooth(x=Year,y=Melt_area,cex=0.7,pch=18,col="red",xlab="Year",ylab="Melt area (km2)")
text(x=2012, y=45000000, "2012 Heat Wave",col="black", font=1.5, cex=0.7)
dev.off()

#Scatterplot annual maximum melt area
setwd("C:/lab/Greenland/")
library(ggplot2)
library(readxl)
data2 <- read_excel("greenland_annual_max.xlsx")
head(data2)
attach(data2)

scatter.smooth(x=year,y=Max_Value,cex=0.7,pch=17,col="purple",xlab="Year",ylab="Melt area max (km2)")
text(x=2012, y=1400000, "2012 Heat Wave",col="black", font=1.5, cex=0.7)
dev.off()

###Albedo

#let's convert HDF file to TIFF
setwd("C:/lab/Greenland/Albedo_Greenland/")
library(gdalUtils)
library(raster)
june2000 <- get_subdatasets("MCD43C3.A2000167.006.2016190203420.hdf")
june2000

#let's choose the layer, convert it and create raster
gdal_translate(june2000[1], dst_dataset = "15june2000_albedo.tif")

rast <- brick("15june2000_albedo.tif")
ext <- c(-0.022,-0.00,0.015,0.025)
zoom(rast,ext=ext)
june2000 <- crop(ext=ext)
writeRaster(june2000,"albedojune2000.tif")
dev.off()

#Albedo variation during summer (example 2000)

setwd("C:/lab/Greenland/Albedo_Greenland")
library(raster)
library(rasterVis)

rlist <- list.files(pattern="albedo")
rlist

import <- lapply(rlist,raster)
Al <- stack(import)

cl <- colorRampPalette(c("red","white","blue"))(100)

par(mfrow=c(1,3))
plot(Al$albedojune2005,col=cl,main="June 2000")
plot(Al$albedojuly2005,col=cl,main="July 2000")
plot(Al$albedoaugust2005,col=cl,main="August 2000")
dev.off()

#boxplot albedo variation 2000 - 2012
setwd("C:/lab/Greenland/Albedo_Greenland/comparison")
library(raster)
library(rasterVis)

rlist <- list.files(pattern="albedo")
rlist
import <- lapply(rlist,raster)
Gr2 <- stack(import)
boxplot(Gr2,main="Albedo comparison 15 august 2000 - 2012")

#plot albedo difference 2012 - 2000
dif <- Gr$albedoaugust2012 - Gr2$albedoaugust2000
cl <- colorRampPalette(c("purple","red","white","blue","black"))(100)
plot(dif,col=cl,main="Albedo variation 15 august 2012-2000")
dev.off()

###Land surface temperature

#Download MODIS data in R 
install.packages("MODIStsp")
library(MODIStsp)
MODIStsp()

#choose the dataset, temporal coverage, projection and data format 

#land surface temperature in summer 2000 - 2020
setwd("C:/lab/Greenland/Temperature/3 months/")
library(raster)
library(rasterVis)
rlist <- list.files(pattern="tif")
rlist
import <- lapply(rlist,raster)
TGr <- stack(import)

cl <- colorRampPalette(c("blue","light blue","pink","red"))(100)
levelplot(TGr,col.regions=cl,main="Summer land surface temperature",names.attr=c("June 2000","July 2000","August 2000","June 2005","July 2005","August 2005","June 2010","July 2010","August 2010","June 2015","July 2015","August 2015","June 2020","July 2020","August 2020"))

#boxplot land surface temperature in july 2000 - 2020
setwd("C:/lab/Greenland/Temperature/July/")
library(raster)
rlist <- list.files(pattern="tif")
rlist
import <- lapply(rlist,raster)
July <- stack(import)
boxplot(July,outline=F, horizontal=FALSE, axes=T, main="Surface temperature variation july 2000-2020",col="red")

#let's see how the temperature changes during the year from 2000 to 2020
#let's calculate the mean of every year

#Total mean temperature in a year (ex. 2000) 
setwd("C:/lab/Greenland/Temperature/Surf_Temp_2000")
library(raster)
rlist <- list.files(pattern="tif")
r_brick2000 <- brick(sapply(rlist,raster))
r_mean2000 <- calc(r_brick2000,mean)
writeRaster(r_mean2000,"r_mean2008.tif")

#Boxplot mean temperature variation
setwd("C:/lab/Greenland/Temperature/Surf_Temp_2000_2020/")
rlist <- list.files(pattern="tif")
rlist
import <- lapply(rlist,raster)

T <- stack(import)
boxplot(T,main="Mean surface temperature 2000-2020")

##all plots were exported as png
png("boxplot_mean_Temp_2000_2020.png")
boxplot(T,main="Mean surface temperature 2000-2020")
dev.off()
