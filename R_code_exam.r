R_code_exam.r

# 1. R code first
# 2. R code Multipanel 
# 3. R code Spatial 
# 4. R code Point pattern analysis: Density map 
# 5. R code for multivariate analysis 
# 6. R code for Remote Sensing analysis 
# 7. R code PCA Remote Sensing
# 8. R code Ecosystem functions 
# 9. R code Radiance
# 10. R code faPAR 
# 11. R code Essential Biodiversity Variables
# 12. R code Climate change indicator monitoring - Snow cover from Copernicus 
# 13. R code Monitoring air pollution (NO2) changes in time
# 14. R code Interpolation 
# 15. R code Species distribution modelling 

# 1. R code first

install.packages("sp") #Classes and methods for spatial data; the classes document where the spatial location information resides, for 2D or 3D data. Utility functions are provided, e.g. for plotting data as maps, spatial selection, as well as methods for retrieving coordinates, for subsetting, print, summary, etc.

library(sp) #loads the package
data(meuse) #loads specified data sets, or list the available data sets
# there is a data set available named meuse

# Let's see how the meuse dataset is structured:
meuse

# let's look at the first row of the set 
head(meuse)

# let's plot two variables
# let's see if the zinc
concentration is related to that of copper
attach(meuse) #the database is attached to the R search path. This means that the database is searched by R when evaluating a variable, so objects in the database can be accessed by simply giving their names.
plot(zinc,copper)
plot(zinc,copper,col="green")
plot(zinc,copper,col="green",pch=19) #pch: different point symbols
plot(zinc,copper,col="green",pch=19,cex=2) #cex: dimension of points

##########################################################################
##########################################################################
##########################################################################

# 2. R code Multipanel 

### Seeing correlation among ecological variables

install.packages("sp")
install.packages("GGally") #The R package 'ggplot2' is a plotting system based on the grammar of graphics. 'GGally' extends 'ggplot2' by adding several functions to reduce the complexity of combining geometric objects with transformed data. Some of these functions include a pairwise plot matrix, a two group pairwise plot matrix, a parallel coordinates plot, a survival plot, and several functions to plot networks.
#Ggally is used for the function ggpairs()

library(sp) # require(sp) will also do the job
library(GGally)

data(meuse) 
attach(meuse)

# Exercise: see the name of the variables and plot cadmium versus zinc
# There are two ways to see the names of the variables:
names(meuse)
head(meuse) 

# exercise: make all the possible pairwise plots of the dataset
#plot(x,cadmium)
#plot(x,zinc)
#plot...
#plot is not a good idea!
pairs(meuse)

#in case you receive the error ''the size is too large" reshape with the mouse graph 

#group the variables to see them better: 
pairs(~cadmium+copper+lead+zinc,data=meuse)

#we can obtain the same result in this way: 
pairs(meuse[,3:6])

# Exercise: prettify this graph
pairs(meuse[,3:6],col="purple",cex=2,pch=13)

#GGally package will prettify the graph
ggpairs(meuse[,3:6])

###########################################################
###########################################################
###########################################################

# 3. R code Spatial

# R code for spatial view of points

library(sp)

data(meuse)

# see only some lines

head(meuse)

# coordinates - we are thinking spatially 
coordinates(meuse)=~x+y    #alt+126 

plot(meuse)

spplot(meuse,"zinc") #Lattice (trellis) plot methods for spatial data with attributes

# Exercise: plot the spatial amount of copper
spplot(meuse, "copper")

# change the title
spplot(meuse, "copper", main="Copper concentration") 

#plot with different size of points
bubble(meuse, "zinc")
bubble(meuse, "zinc", main="Zinc concentration")

# Exercise: bubble copper in red 
bubble(meuse, "copper",main="Copper concentration", col="red")

### Importing new data 

# download covid_agg-csv from our teaching site and builds a folder called lab into C:
# put the covid_agg.csv file into the folder lab

# setting the working directory: lab
# Windows 
setwd("C:/lab/")

# name of the data set and the first row of the table is not a datum, there is a header
# TRUE can be written T 
covid <- read.table("covid_agg.csv", head=TRUE)
head(covid)

# let's plot considering the number of cases per country (not spatial)
attach(covid) 
plot(country,cases)

# if you do not attach covid plot(covid$country,covid$cases)

#change orientation of the axis to see all the countries - vertical 
plot(country, cases, las=0)

#nothing changed: parallel labels 
plot(country, cases, las=1) # horizontal labels
plot(country, cases, las=2) # perpendicular labels 
plot(country, cases, las=3) # vertical labels

plot(country, cases, las=3, cex.axis=0.5)

# let's plot them spatially 

# ggplot2 package
install.packages("ggplot2")
library(ggplot2) # require(ggplot2)

# save the .RData under the menu File

# load the previously saved .RData

# setting the working directory: lab 
setwd("C:/lab/")
load("spatial.RData")

# let's see if the file was uploaded correctly, ls is a list of objects
ls()
# covid

library(ggplot2) #require(ggplot2)

data(mpg)
head(mpg)
#key components: data, aes, geometry 
ggplot(mpg, aes(x=displ, y=hwy)) + geom_point()

# change the geometry
ggplot(mpg, aes(x=displ, y=hwy)) + geom_line()
ggplot(mpg, aes(x=displ, y=hwy)) + geom_polygon()

head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point() #size to show the points with higher size with higher number of cases

###############################################################
###############################################################
###############################################################

# 4. R code Point pattern analysis: Density map 

install.packages("spatstat") #toolbox for analysing Spatial Point Patterns
library(spatstat)

attach(covid)
head(covid)

# Create a point pattern dataset in the two-dimensional plane
covids <- ppp(lon, lat, c(-180,180), c(-90,90)) 

# duccio <- c(12,34,55,77,88,89) cluster all together

#without attaching the covid set
#covids <- ppp(covid$lon, covid$lat, c(-180,180), c(-90,90))

#build the density map
d <- density(covids)

plot(d)

# add points to a plot
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

install.packages("rgdal") #Provides bindings to the 'Geospatial' Data Abstraction Library 
library("rgdal") 

# let's input vector lines (x0y0, x1y1, x2,y2..)
coastlines <- readOGR("ne_10m_coastline.shp") #reads an OGR data source and layer into a suitable Spatial vector object

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

# export PDF
pdf("covid_density.pdf")
cl <- colorRampPalette(c("blue","green","yellow")) (100)
plot(d, col=cl, main="Densities of covid-19") 
points(covids)
plot(coastlines, add=T)

# close the specified plot
dev.off()

# export png
png("covid_density.png")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()

####################################################
####################################################
####################################################

# 5. R code for multivariate analysis 

setwd("C:/lab/") #specify the path to the desired folder of your working directory

install.packages("vegan") #The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. Most of its multivariate tools can be used for other data types as well.

library(vegan)

biomes <- read.table("biomes.csv", head=T, sep=",") 
#head=T or F:a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: head is set to TRUE if and only if the first row contains one fewer field than the number of columns.
#sep=",": values on each line of the file are separated by the character ,

head(biomes) # View(biomes), biomes

# DEtrended CORrespondence ANAlysis (decorana)
multivar <- decorana(biomes)
plot(multivar)

multivar

plot(multivar)

# biomes types
biomes_types <- read.table("biomes_types.csv", head=T, sep=",")
head(biomes_types)

attach(biomes_types)

# ordiellipse: elipse connecting all the points
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3) #lwd: line width 
# col=c("green", "blue","red", "black")

#points connected in a spider shape
ordispider(multivar, type, col=1:4, label = T)

##############################################################
##############################################################
##############################################################

# 6. R code for Remote Sensing analysis 

setwd("C:/lab/")

# install.packages("raster") #Reading, writing, manipulating, analyzing and modeling of spatial data. The package implements basic and high-level functions for raster data and for vector data operations such as intersections.
install.packages("RStoolbox") #Toolbox for remote sensing image processing and analysis such as calculating spectral indices, principal component transformation, unsupervised and supervised classification or fractional cover analyses.

# let's install faster 
# install.packages(c("raster","RStoolbox"))

library(raster)

# import image and create a multi-layer raster object
p224r63_2011 <- brick("p224r63_2011_masked.grd")

plot(p224r63_2011)
cl <- colorRampPalette(c('black','grey','light grey'))(100) 

# Exercise: plot the image with the new color ramp palette
plot(p224r63_2011, col=cl)

#bands of Landsat
# B1: blue
# B2: green
# B3: red
# B4: NIR

#multiframe of different plots
par(mfrow=c(2,2))

# let's see all the colors of R
colors()

# B1: blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) 
plot(p224r63_2011$B1_sre, col=clb)

# B2 green
clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B2_sre, col=clg)

#B3 red
clr <- colorRampPalette(c('dark red','red','pink'))(100)
plot(p224r63_2011$B3_sre, col=clr)

#B4 NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

# let's have 4 rows and 1 column
# par(mfrow=c(4,1))

# let's close the window
dev.off()

# RGB 
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

# Exercise: mount NIR on top of the G component of RGB
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin")

# NIR infrared into the blue component
plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin")

#####################################################

setwd("C:/lab/")
load("rs.RData")
ls()

library(raster)

# 1988 image
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plot(p224r63_1988)

# Exercise: plot in visible RGB 321 both images (1988 and 2011)
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

# Exercise: plot in false colour RGB 432 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

#enhance the noise!
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist")

#plotRGB
#bands of Landsat
# B1: blue
# B2: green
# B3: red: B3_sre
# B4: NIR: B4_sre

dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
cl <- colorRampPalette(c("darkorchid3","light blue","lightpink4"))(100)
plot(dvi2011,col=cl)

#Exercise: dvi for 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
cl <- colorRampPalette(c("darkorchid3","light blue","lightpink4"))(100)
plot(dvi1988,col=cl)

# difference from 1988 to 2011
diff <- dvi2011 - dvi1988
plot(diff)

# changing the grain
# resampling
p224r63_2011res <- aggregate(p224r63_2011, fact=10)
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)

par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")

# let's see the resolution of p224r63_2011 -> 30 m
p224r63_2011

#and the resolution of p224r63_2011res100 -> 3 km
p224r63_2011res100

##############################################################
##############################################################
##############################################################

# 6. Ecosystem functions in R 

# R code to view biomass over the world and calculate changes in ecosystem functions
# energy
# chemical cycling
# proxies

install.packages("rasterdiv") #Providing functions to calculate indices of diversity on numerical matrices based on information theory
install.packages("rasterVis") #Methods for enhanced visualization and interaction with raster data. It implements visualization methods for quantitative data and categorical data, both for univariate and multivariate rasters. It also provides methods to display spatiotemporal rasters, and vector fields

library(rasterVis)
library(rasterdiv)

data(copNDVI)
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) #Reclassify values of a Raster* object. The function (re)classifies groups of values to other values. For example, all values between 1 and 10 become 1, and all values between 11 and 15 become 2
#cbind:Take a sequence of vector, matrix or data-frame arguments and combine by columns or rows, respectively.

levelplot(copNDVI)

copNDVI10 <- aggregate(copNDVI, fact=10) #Splits the data into subsets, computes summary statistics for each, and returns the result in a convenient form
levelplot(copNDVI10)

copNDVI100 <- aggregate(copNDVI, fact=100)
levelplot(copNDVI100)

##################################### gift
#library(ggplot2)

#myPalette <- colorRampPalette(c('white','green','dark green'))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

# ggR(copNDVI, geom_raster = TRUE) +
#scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
#labs(x="Longitude",y="Latitude", fill="")+
#   theme(legend.position = "bottom") +
#  NULL
# +
# ggtitle("NDVI")
######################################

setwd("C:/lab/")

library(raster)

defor1 <- brick("defor1_.jpg")
defor2 <- brick("defor2_.jpg")

# see info about the image
defor1

#band1: NIR, defor1_.1
#band2: red, defor1_.2
#band3: green

plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

dvi1 <- defor1$defor1_.1 - defor1$defor1_.2

# defor2
# band1: NIR, defor2_.1
# band2: red, defor2_.2

dvi2 <- defor2$defor2_.1 - defor2$defor2_.2

cl <- colorRampPalette(c("darkblue", "yellow", "red", "black"))(100) # specifying a color scheme 

par(mfrow=c(1,2))
plot(dvi1, col=cl)
plot(dvi2, col=cl)

difdvi <- dvi1 - dvi2
### Warning in dvi1 - dvi2: Raster objects have different extents. Result for their intersection is returned

dev.off() #This function closes the specified plot (by default the current device) 

cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col=cld)

hist(difdvi) #The generic function hist computes a histogram of the given data values

#################################################################
#################################################################
#################################################################

# 7. R code PCA Remote Sensing

# install.packages("RStoolbox")

setwd("C:/lab/")

library(raster)
library(RStoolbox)
library(ggplot2)

#now we use the function brick that is used to import the whole image of the satellite. 
p224r63_2011 <- brick("p224r63_2011_masked.grd")

# b1 blue
# b2 green 
# b3 red
# b4 NIR
# b5 SWIR 
# b6 thermal infrared
# b7 SWIR
# b8 panchromatic

# RGB
# now we plot the image in the RGB space
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin") #stretch: Provide the desired output range (minv and maxv) and the lower and upper bounds in the original data
ggRGB(p224r63_2011, 5, 4, 3) #ggRGB: Calculates RGB color composite raster for plotting with ggplot2

# do the same with the 1988 image 
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")

par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011)
# "B1_sre" "B2_sre" "B3_sre" "B4_sre" "B5_sre" "B6_bt"  "B7_sre"

plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)

# PCA
# decrease the resolution 
p224r63_2011_res <-aggregate(p224r63_2011, fact=10)

#library(RStoolbox)
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)

cl <- colorRampPalette(c("dark grey", "grey", "light grey")) (100)
plot(p224r63_2011_pca$map, col=cl)

summary(p224r63_2011_pca$model)
# PC1 99.83% of the whole variation 

pairs(p224r63_2011)
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin")

# 1988
p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res)
plot(p224r63_1988_pca$map, col=cl)

summary(p224r63_1988_pca$model)
# 99.56% by PC1

pairs(p224r63_1988)

# now we can make a difference between the 1988 and 2011 and then plot the difference. 
# we are making the difference of every pixel
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)

# since the PC1 contains most of the information, we can only plot this one, only 1 layer
cldif <- colorRampPalette(c('blue','black','yellow'))(100)
plot(difpca$PC1, col=cldif)
# we see the areas that have changed most
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map

plot(difpca)
cldif <- colorRamppalette(c("blue", "black", "yellow")) (100)

#####################################################################
#####################################################################
#####################################################################

# 8. R code Ecosystem functions 

setwd("C:/lab/")
library(raster)

snt <- brick("snt_r10.tif")
snt

plot(snt)

# B1 blue
# B2 green
# B3 red
# B4 NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin")

plotRGB(snt, 4, 3, 2, stretch="Lin")
pairs(snt)

# PCA analysis
library(RStoolbox)

sntpca <- rasterPCA(snt)
sntpca

summary(sntpca$model)
# 70% of information
plot(sntpca$map)

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")

# set the moving window - diversity measurement
window <- matrix(1, nrow = 5, ncol = 5)
window

# focal function - calculation in a certain extent 
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)
cl <- colorRampPalette(c('dark blue','green','orange','red'))(100)
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt, 4, 3, 2, stretch="Lin", main= "original image")
plot(sd_snt, col=cl, main="diversity")

################################# day 2 Cladonia example

setwd("C:/lab/")
library(raster)

clad <- brick("cladonia_stellaris_calaita.JPG")

plotRGB(clad, 1,2,3, stretch="lin")
 
window <- matrix(1, nrow = 3, ncol = 3)
window

### PCA analysis 
cladpca <- rasterPCA(clad)

cladpca

summary(cladpca$model)
# 98% Comp 1 (because it's a visible spectrum)

 plotRGB(cladpca$map, 1, 2, 3, stretch="lin")

# set the moving window 
window <- matrix(1, nrow = 5, ncol = 5)
window

sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd)

PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)
# graph left: original cladonia set - graph righ: aggregated set 

# plot the calculation 
par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)
# plot(sd_clad_agg, col=cl)

#############################################################
#############################################################
#############################################################

# 9. R code Radiance
 
library(raster)

toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)

# put the values into the raster toy 
# put each data into every pixel
values(toy) <- c(1.13,1.44,1.55,3.4)

plot(toy)
text(toy, digits=2)

# 4 potential values (2 bits) 
toy2bits <- stretch(toy,minv=0,maxv=3)
 
# avoid decimals, only integers to gain space
storage.mode(toy2bits[]) = "integer"
plot(toy2bits)
text(toy2bits, digits=2)

# range of 4 bits (16 values from 0 to 15)
toy4bits <- stretch(toy,minv=0,maxv=15)
storage.mode(toy4bits[]) = "integer" 
plot(toy4bits)
text(toy4bits, digits=2)

# 8 bits (256 potential values)
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer" 
plot(toy8bits)
text(toy8bits, digits=2)

par(mfrow=c(1,4))

plot(toy)
text(toy, digits=2)

plot(toy2bits)
text(toy2bits, digits=2)

plot(toy4bits)
text(toy4bits, digits=2)

plot(toy8bits)
text(toy8bits, digits=2)

dev.off
library(rasterdiv)
plot(copNDVI)
copNDVI # values : 0, 255  (min, max) # they are using 8 bits

###########################################################
###########################################################
###########################################################

# 10. R code faPAR 

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

######################################################################
######################################################################
######################################################################

# 11. R code Essential Biodiversity Variables

setwd("C:/lab/")
library(raster)

snt <- brick("snt_r10.tif")
snt

plot(snt)

# B1 blue
# B2 green
# B3 red
# B4 NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin")

plotRGB(snt, 4, 3, 2, stretch="Lin")
pairs(snt)

# PCA analysis
library(RStoolbox)

sntpca <- rasterPCA(snt)
sntpca

summary(sntpca$model)
# 70% of information
plot(sntpca$map)

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")

# set the moving window - diversity measurement
window <- matrix(1, nrow = 5, ncol = 5)
window

# focal function - calculation in a certain extent 
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)
cl <- colorRampPalette(c('dark blue','green','orange','red'))(100)
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt, 4, 3, 2, stretch="Lin", main= "original image")
plot(sd_snt, col=cl, main="diversity")

################################# day 2 Cladonia example

setwd("C:/lab/")
library(raster)

clad <- brick("cladonia_stellaris_calaita.JPG")

plotRGB(clad, 1,2,3, stretch="lin")
 
window <- matrix(1, nrow = 3, ncol = 3)
window

### PCA analysis 
cladpca <- rasterPCA(clad)

cladpca

summary(cladpca$model)
# 98% Comp 1 (because it's a visible spectrum)

 plotRGB(cladpca$map, 1, 2, 3, stretch="lin")

# set the moving window 
window <- matrix(1, nrow = 5, ncol = 5)
window

sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd)

PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)
# graph left: original cladonia set - graph righ: aggregated set 

# plot the calculation 
par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)
# plot(sd_clad_agg, col=cl)

########################################################################
########################################################################
########################################################################
 
# 12. R code Climate change indicator monitoring - Snow cover from Copernicus 

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

### [1] "snow2000r.tif" "snow2005r.tif" "snow2010r.tif" "snow2015r.tif" "snow2020r.tif"

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

## [1] "snow2000r.tif" "snow2005r.tif" "snow2010r.tif" "snow2015r.tif" "snow2020r.tif"

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

# export 2
png("my_final_exciting_graph.png")
plot(final.stack, col=cl)
dev.off()

###########################################################
###########################################################
###########################################################

# 13. R code Monitoring air pollution (NO2) changes in time

setwd("C:/lab/NO2/")
library(raster)

rlist <- list.files(pattern="EN")
rlist
import <- lapply(rlist, raster)

EN <- stack(import)
cl <- colorRampPalette(c('yellow','green','blue'))(100)
plot(EN, col=cl)

# january and march
par(mfrow=c(1,2))
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)
dev.off

# RGB space
# r(red)= EN_0001 - higher values in january; g= EN_0007 - higher values at mid time; b= EN_0013 - higher values in march)
plotRGB(EN, r=1, g=7, b=13, stretch="lin")

# difference map
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) # 
plot(dif, col=cld) # red= higher differences between the 2 images (first and final one)

# quantitative estimate
boxplot(EN)
# remove the ouliers
boxplot(EN,outline=F)

boxplot(EN,outline=F, horizontal=T, axes=T)
 
# plot values of the pixels
plot(EN$EN_0001, EN$EN_0013)
abline(0,1, col="red")    # y= a + bx -> a passes through zero, b is 1

###################### 1:1 line with snow data

setwd("C:/lab/snow/")

rlist <- list.files(pattern="snow20")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
 
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1, col=red)

#################################################################
#################################################################
#################################################################

# 14. R code Interpolation

## Interpolation: spatstat library
setwd("C:/lab/")

# library(dbmss)
library(spatstat)

# Beech forests Martina
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T)
attach(inp)
plot(X,Y)

inppp <- ppp(x=X,y=Y,c(716000,718000),c(4859000,4861000))
marks(inppp) <- Canopy.cov
canopy <- Smooth(inppp)

plot(canopy)

marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)

par(mfrow=c(1,2))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)

#########
# Dati psammofile Giacomo
inp.psam <- read.table("dati_psammofile.csv", sep=";", head=T)
attach(inp.psam)
summary(inp.psam)

plot(E,N)
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))

marks(inp.psam.ppp) <- C_org
C <- Smooth(inp.psam.ppp)

plot(C)
points(inp.psam.ppp)

#########################################################
#########################################################
#########################################################

# 15. R code Species distribution modelling 

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































