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
