# R_code_ecosystem_functions.r

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





