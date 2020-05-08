# R_code_faPAR.r
# How to look at chemical cycling from satellites

# t(copNDVI)
#
# setwd("C:/Downloads/")
#
#..?
#

setwd("C:/lab/")
library(raster)
library(rasterVis)
library(rasterdiv)

plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # removing the data from 253 to 255 (water)
levelplot(copNDVI)

faPAR10 <- raster("C:/lab/faPAR10.tif")
levelplot(faPAR10)

# save the images as PDF
pdf("copNDVI.pdf")
levelplot(copNDVI)
dev.off()

pdf("faPAR10.pdf")
levelplot(faPAR10)
dev.off()



