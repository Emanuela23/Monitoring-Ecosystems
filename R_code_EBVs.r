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



