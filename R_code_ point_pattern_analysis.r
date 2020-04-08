# Point pattern analysis: Density map

install.packages("spatstat")
library(spatstat)

attach(covid)
head(covid)

covids <- ppp(lon, lat, c(-180,180), c(-90,90))
#?ppp put ? to understand 

# duccio <- c(12,34,55,77,88,89) cluster all together

#without attaching the covid set
#covids <- ppp(covid$lon, covid$lat, c(-180,180), c(-90,90))

#build the density map
d <- density(covids)

plot(d)
points(covids)

#----

setwd("C:/lab/")
load("spatial.RData")
ls()

# covids: point pattern 
# d: density map
library(spatstat)

plot(d)
points(covids)

# let's see where the points are located - let's add the coastlines 

install.packages("rgdal")
library("rgdal") 

# let's input vector lines (x0y0, x1y1, x2,y2..)
coastlines <- readOGR("ne_10m_coastline.shp")

# let's add the lines to the previous image
plot(coastlines, add=T)

# change the color and make the map beautiful 

cl <- colorRampPalette(c("yellow", "orange", "red")) (100)
# cl <- colorRampPalette(c('yellow', 'orange', 'red'))(100) #
# 100: all colors from yellow to red (red for the hotspots)

plot(d, col=cl) # cl: is the colorRampPalette
points(covids)
plot(coastlines, add=T) 

# Exercise: new colorRampPalette

cl <- colorRampPalette(c("blue","green","yellow")) (100)
plot(d, col=cl, main="Densities of covid-19") 
points(covids)
plot(coastlines, add=T)

# example of prof: 
# clr <- colorRampPalette(c("light green", "yellow", "orange", "violet")) (100)

# export
pdf("covid_density.pdf")
cl <- colorRampPalette(c("blue","green","yellow")) (100)
plot(d, col=cl, main="Densities of covid-19") 
points(covids)
plot(coastlines, add=T)
dev.off()

# export png
png("covid_density.png")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()

