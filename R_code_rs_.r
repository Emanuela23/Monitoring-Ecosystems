# R code for Remote sensing analysis 

setwd("C:/lab/")

# install.packages("raster")
install.packages("RStoolbox")

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
plot(dvi2011)

##################################23
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre

# DVI for the two years: compare with a difference in time
# NIR - RED
# NDVI = (NIR - RED) / (NIR - RED)

dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre

par(mfrow=c(2,1))
plot(dvi1988)
plot(dvi2011)


par(mfrow=c(2,1))
cldvi <- colorRampPalette(c('red','orange','yellow'))(100) # 
plot(dvi1988, col=cldvi)
plot(dvi2011, col=cldvi)


# difference in time
difdvi <- dvi2011 - dvi1988
cldif <- colorRampPalette(c('blue','white','red'))(100) #
plot(difdvi, col=cldif)



# install.packages("RStoolbox")
library(RStoolbox)

# PCA
p224r63_2011res <- aggregate(p224r63_2011, fact=10)

par(mfrow=c(2,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")

p224r63_2011_pca <- rasterPCA(p224r63_2011res)

summary(p224r63_2011_pca$model) 
plotRGB(p224r63_2011_pca$map, r=4, g=3, b=2, stretch="Lin")

plot(p224r63_2011_pca$map)

# land cover
p224r63_2011c <- unsuperClass(p224r63_2011, nClasses=5)
clclass <- colorRampPalette(c('red', 'green', 'yellow', 'blue', 'black'))(100) 
plot(p224r63_2011c$map, col=clclass)



################
setwd("C:/lab/")

load("rs.RData")
ls()

p224r63_1988_masked

p224r63_1988 <- brick("p224r63_1988_masked.grd") 

# Exercise: plot in visible RGB 321 both images

# B1: blue
# B2: green
# B3: red
# B4: NIR

par(nfor=c(2,1))

plotRGB(p224r63, r=3, g=2, b=1, stretch="lin")



