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

