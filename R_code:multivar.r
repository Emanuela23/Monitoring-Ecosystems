# R code for multivariate analysis 

setwd("C:/lab/")

install.packages("vegan")
library(vegan)

biomes <- read.table("biomes.csv", head=T, sep=",")
head(biomes) # View(biomes). biomes

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

