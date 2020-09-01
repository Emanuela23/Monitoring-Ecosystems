### Multi panel in R: the second lecture of Monitoring ecosystem 

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


