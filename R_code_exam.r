R_code_exam.r

# 1. R code first
# 2. Multipanel in R 
# 3. Spatial R 
# 4. R code for multivariate analysis 
# 5. R code for Remote sensing analysis 
# 6. Ecosystem functions in R 

# 1. R code first

install.packages("sp") #Classes and methods for spatial data; the classes document where the spatial location information resides, for 2D or 3D data. Utility functions are provided, e.g. for plotting data as maps, spatial selection, as well as methods for retrieving coordinates, for subsetting, print, summary, etc.

library(sp) #loads the package
data(meuse) #loads specified data sets, or list the available data sets.
# there is a data set available named meuse

# Let's see how the meuse dataset is structured:
meuse

# let's look at the first row of the set 
head(meuse)

# let's plot two variables
# let's see if the zin concentration is related to that of copper
attach(meuse) #the database is attached to the R search path. This means that the database is searched by R when evaluating a variable, so objects in the database can be accessed by simply giving their names.
plot(zinc,copper)
plot(zinc,copper,col="green")
plot(zinc,copper,col="green",pch=19) #pch: different point symbols
plot(zinc,copper,col="green",pch=19,cex=2) #cex: dimension of points

##########################################################################
##########################################################################
##########################################################################

### 2. Multi panel in R

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

# 3. Spatial R

# R code for spatial view of points

library(sp)

data(meuse)

# see only some lines

head(meuse)

# coordinates - we are thinking spatially 
coordinates(meuse)=~x+y    #alt+126 

plot(meuse)

spplot(meuse,"zinc")

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

# change the geom
ggplot(mpg, aes(x=displ, y=hwy)) + geom_line()
ggplot(mpg, aes(x=displ, y=hwy)) + geom_polygon()

head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point() #size to show the points with higher size with higher number of cases

###############################################################
###############################################################
###############################################################

# 4. R code for multivariate analysis 

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
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3)
# col=c("green", "blue","red", "black")

ordispider(multivar, type, col=1:4, label = T)

##############################################################
##############################################################
##############################################################

# 5. R code for Remote sensing analysis 

setwd("C:/lab/")

# install.packages("raster") #Reading, writing, manipulating, analyzing and modeling of spatial data. The package implements basic and high-level functions for raster data and for vector data operations such as intersections.
install.packages("RStoolbox") #Toolbox for remote sensing image processing and analysis such as calculating spectral indices, principal component transformation, unsupervised and supervised classification or fractional cover analyses.

# let's install faster 
# install.packages(c("raster","RStoolbox"))

library(raster)

# import image
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















