# R code for multivariate analysis 

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

