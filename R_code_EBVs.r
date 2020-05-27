# R_code_EBVs.r

setwd("C:/lab/")
library(raster)

snt <- brick("snt_r10.tif")
snt

plot(snt)

# B1 blue
# B2 green
# B3 red
# B4 NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin")

plotRGB(snt, 4, 3, 2, stretch="Lin")
pairs(snt)

# PCA analysis
library(RStoolbox)

sntpca <- rasterPCA(snt)
sntpca

summary(sntpca$model)
# 70% of information
plot(sntpca$map)

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")

# set the moving window - diversity measurement
window <- matrix(1, nrow = 5, ncol = 5)
window

# focal function - calculation in a certain extent 
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)
cl <- colorRampPalette(c('dark blue','green','orange','red'))(100)
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt, 4, 3, 2, stretch="Lin", main= "original image")
plot(sd_snt, col=cl, main="diversity")

################################# day 2 Cladonia example

setwd("C:/lab/")
library(raster)

clad <- brick("cladonia_stellaris_calaita.JPG")

plotRGB(clad, 1,2,3, stretch="lin")
 
window <- matrix(1, nrow = 3, ncol = 3)
window

### PCA analysis 
cladpca <- rasterPCA(clad)

cladpca

summary(cladpca$model)
# 98% Comp 1 (because it's a visible spectrum)

 plotRGB(cladpca$map, 1, 2, 3, stretch="lin")

# set the moving window 

sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd)

PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)
# graph left: original cladonia set - graph righ: aggregated set 

# plot the calculation 
par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)
# plot(sd_clad_agg, col=cl)


 

