# R_code_no2.r

setwd("C:/lab/NO2/")
library(raster)

rlist <- list.files(pattern="EN")
import <- lapply(rlist, raster)

EN <- stack(import)
cl <- colorRampPalette(c('yellow','green','blue'))(100)
plot(EN, col=cl)

# january and march
par(mfrow=c(1,2))
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)
dev.off

# RGB space
# r(red)= EN_0001 - higher values in january; g= EN_0007 - higher values at mid time; b= EN_0013 - higher values in march)
plotRGB(EN, r=1, g=7, b=13, stretch="lin")

# difference map
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) # 
plot(dif, col=cld) # red= higher differences between the 2 images (first and final one)

# quantitative estimate
boxplot(EN)
# remove the ouliers
boxplot(EN,outline=F)

boxplot(EN,outline=F, horizontal=T, axes=T)
 
# plot values of the pixels
plot(EN$EN_0001, EN$EN_0013)
abline(0,1, col="red")    # y= a + bx -> a passes through zero, b is 1

###################### 1:1 line with snow data

setwd("C:/lab/snow/")

rlist <- list.files(pattern="snow20")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
 
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1, col=red)




