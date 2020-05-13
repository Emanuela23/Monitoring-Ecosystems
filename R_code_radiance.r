# R_code_radiance.r

library(raster)

toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)

# put the values into the raster toy 
# put each data into every pixel
values(toy) <- c(1.13,1.44,1.55,3.4)

plot(toy)
text(toy, digits=2)

# 4 potential values (2 bits) 
toy2bits <- stretch(toy,minv=0,maxv=3)
 
# avoid decimals, only integers to gain space
storage.mode(toy2bits[]) = "integer"
plot(toy2bits)
text(toy2bits, digits=2)

# range of 4 bits (16 values from 0 to 15)
toy4bits <- stretch(toy,minv=0,maxv=15)
storage.mode(toy4bits[]) = "integer" 
plot(toy4bits)
text(toy4bits, digits=2)

# 8 bits (256 potential values)
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer" 
plot(toy8bits)
text(toy8bits, digits=2)

par(mfrow=c(1,4))

plot(toy)
text(toy, digits=2)

plot(toy2bits)
text(toy2bits, digits=2)

plot(toy4bits)
text(toy4bits, digits=2)

plot(toy8bits)
text(toy8bits, digits=2)

dev.off
library(rasterdiv)
plot(copNDVI)
copNDVI # values : 0, 255  (min, max) # they are using 8 bits





