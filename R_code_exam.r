R_code_exam.r

# 1. R code first
# 2. Spatial R 

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

# 2. Spatial R

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




