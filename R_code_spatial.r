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





