rm(list=ls())

library(maptools)
library(rgdal)
library(rgl)
library(sp)
library(foreign)
library(raster)
library(MBA)
library(RColorBrewer)
#library(fields)
## raster ##
slope <-    raster("/media/jason/Elements/Alaska_30m_data/DEM_Slope_degree_30m_Projected.img")
temp <-     raster("/media/jason/Elements/Alaska_30m_data/ak_tavg_14_30m.img")
heatload <- raster("/media/jason/Elements/Alaska_30m_data/dem_heatload_30m_Projected.img")
cti <-      raster("/media/jason/Elements/Alaska_30m_data/dem_cti_30m_projected.img")
texture <-  raster("/media/jason/Elements/Alaska_30m_data/AK_STexture_Jorgenson.img")

## Alaska shapefile
AK <- getData("GADM", country = "USA", level = 1)
AK <- AK[which(AK$NAME_1 == "Alaska"),]
## ## give projection information to Alaska data ##
AK <- spTransform(AK,CRS(proj4string(heatload)))
e <- extent(-1100000, 1541652, 399008.5, 2403609)
AK.crop <- crop(AK,e)

AK.pf <- read.csv("MegaAlaska.txt")
## reproject ##
coordinates(AK.pf) <- ~long_site+lat_site
proj4string(AK.pf) <- proj4string(cti)
coords.pf <- coordinates(AK.pf)
coords.pf <- project(coords.pf, proj4string(cti))
colnames(coords.pf) <- c("x","y")
AK.data <- slot(AK.pf,"data")
alaska.permafrost <- cbind(AK.data,coords.pf)
coordinates(alaska.permafrost) <- ~x+y
proj4string(alaska.permafrost) <- proj4string(cti)

## setEPS()
## postscript('akpf.eps', fonts="Courier")
png('png/AKpf.png', width=585, height=585, family="Courier")
postscript('png/AKpf.eps', width=8, height=8,horizontal = FALSE, onefile = FALSE, paper = "special",pointsize=18, family = "ComputerModern", encoding = "TeXtext.enc")
plot(AK.crop, xaxt="n", yaxt="n")
title(xlab="Easting (km)", ylab="Northing (km)")#, family="Courier")
points(alaska.permafrost, col=ifelse(alaska.permafrost$pfSCORE>0,"blue","red"))
legend(830000,2300000, pch=c(1,1),
       col=c("blue","red"),box.lwd=0,box.col="white",
       legend=c("Presence","Absence"))
box(lty="solid")
xticks=c(-1150000, -650000, -150000, 350000, 850000, 1350000)
xlab=c(0, 500, 1000, 1500, 2000, 2500)
yticks=c(180000, 680000, 1180000, 1680000, 2180000,2635000)
ylab=c(0, 500, 1000, 1500, 2000,2500)
axis(1, at=xticks, labels=xlab)
axis(2, at=yticks, labels=ylab)
dev.off()
## cti.1 <- cti[AK.crop,]
## plot(cti.1)

## dev.off()
png("png/slope.png", width=585, height=585, family="Courier")
postscript('png/slope.eps', width=8, height=8,horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=18, family = "ComputerModern", encoding = "TeXtext.enc")
plot(slope, xaxt="n", yaxt="n", col=tim.colors(100))
title(xlab="Easting (km)", ylab="Northing (km)")
xticks=c(-1000000, -500000, 0, 500000, 1000000, 1500000)
xlab=c(0, 500, 1000, 1500, 2000, 2500)
yticks=c(150000,650000, 1150000, 1650000, 2150000, 2650000)
ylab=c(0, 500, 1000, 1500, 2000, 2500)
axis(1, at=xticks, labels=xlab)
axis(2, at=yticks, labels=ylab)
dev.off()


png("png/temp.png", width=585, height=585, family="Courier")
postscript('png/temp.eps', width=8, height=8,horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=18, family = "ComputerModern", encoding = "TeXtext.enc")
plot(temp, xaxt="n", yaxt="n", col=tim.colors(100))
title(xlab="Easting (km)", ylab="Northing (km)")
xticks=c(-1000000, -500000, 0, 500000, 1000000, 1500000)
xlab=c(0, 500, 1000, 1500, 2000, 2500)
yticks=c(150000,650000, 1150000, 1650000, 2150000, 2650000)
ylab=c(0, 500, 1000, 1500, 2000, 2500)
axis(1, at=xticks, labels=xlab)
axis(2, at=yticks, labels=ylab)
dev.off()

rnb <- reclassify(heatload, cbind(255,NA))

rna <- reclassify(cti, cbind(94, NA))
rc <- reclassify(rna, c(150,255,NA))

## png("png/heatload.png", width=585, height=585, family="Courier")
postscript('png/heatload.eps', width=8, height=8, ,horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=18, family = "ComputerModern", encoding = "TeXtext.enc")
plot(rnb, xaxt="n", yaxt="n")
title(xlab="Easting (km)", ylab="Northing (km)")
xticks=c(-1000000, -500000, 0, 500000, 1000000, 1500000)
xlab=c(0, 500, 1000, 1500, 2000, 2500)
yticks=c(150000,650000, 1150000, 1650000, 2150000, 2650000)
ylab=c(0, 500, 1000, 1500, 2000, 2500)
axis(1, at=xticks, labels=xlab)
axis(2, at=yticks, labels=ylab)
dev.off()

## png("png/cti.png", width=585, height=585, family="Courier")
postscript('png/cti.eps', width=8, height=8,horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=18, family = "ComputerModern", encoding = "TeXtext.enc")
plot(rc, xaxt="n", yaxt="n")
title(xlab="Easting (km)", ylab="Northing (km)")
xticks=c(-1000000, -500000, 0, 500000, 1000000, 1500000)
xlab=c(0, 500, 1000, 1500, 2000, 2500)
yticks=c(150000,650000, 1150000, 1650000, 2150000, 2650000)
ylab=c(0, 500, 1000, 1500, 2000, 2500)
axis(1, at=xticks, labels=xlab)
axis(2, at=yticks, labels=ylab)
dev.off()


  
## extract values ##
## ##
## slope <- as.data.frame(extract(slope, coords.pf))
## temp <- as.data.frame(extract(temp, coords.pf))
## heatload <- as.data.frame(extract(heatload, coords.pf))
## cti <- as.data.frame(extract(cti, coords.pf))
## texture <- as.data.frame(extract(texture, coords.pf))
## ##
## colnames(slope) <- c("slope")
## colnames(temp) <- c("temp")
## colnames(heatload) <- c("heatload")
## colnames(cti) <- c("cti")
## colnames(texture) <- c("texture")
## permafrost <- as.data.frame(AK.pf$pfSCORE)
## colnames(permafrost) <- c("permafrost")
##
## texture <- as.numeric(as.matrix(texture))
## texture <- as.data.frame(factor(texture, labels=c("Sandy","Silty","Rocky","Water")))
## colnames(texture) <- c("texture")
## write csv
## write.csv(slope, "covariates/slope.csv", row.names=FALSE)
## write.csv(temp, "covariates/temp.csv", row.names=FALSE)
## write.csv(heatload, "covariates/heatload.csv", row.names=FALSE)
## write.csv(cti, "covariates/cti.csv", row.names=FALSE)
## write.csv(texture, "covariates/texture.csv", row.names=FALSE)
## write.csv(permafrost, "covariates/permafrost.csv", row.names=FALSE)
## write.csv(coords.pf, "covariates/coords.csv", row.names=FALSE)

## ################
## ## start here ##
## ################
## ## Not run: 
## rmvn <- function(n, mu=0, V = matrix(1)){
##   p <- length(mu)
##   if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
##   D <- chol(V)
##   t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
## }

## load data ##
slope <- read.csv("covariates/slope.csv")
temp <- read.csv("covariates/temp.csv")
heatload <- read.csv("covariates/heatload.csv")
cti <- read.csv("covariates/cti.csv")
permafrost <- read.csv("covariates/permafrost.csv")
coords <- read.csv("covariates/coords.csv")

## get complete cases ##
dat <- cbind(permafrost,slope,temp,heatload,cti,coords)
dat <- dat[complete.cases(dat),]

## get y ##
y <- as.data.frame(dat[,"permafrost"])
names(y) <- c("permafrost")

## get and scale x ##
x <- dat[,2:5]
#x <- scale(x)
x <- cbind(1,x)
colnames(x) <- c("Intercept","slope","temp","heatload","cti")
coords <- dat[,6:7]

## remove duplicate coordinates
data.frame <- as.data.frame(cbind(x,coords))
coordinates(data.frame) <- data.frame[,5:6]
zd <- zerodist(data.frame)
x <- x[-zd[,2],]
coords <- coords[-zd[,2],]
y <- y[-zd[,2],]

## get n nd holdout rows ##
n <- nrow(x)
ho <- sample(1:n, floor(0.7849*n))
x.dat <- as.data.frame(x)

## slope.x <- cbind(coords,x.dat$slope)
## temp.x <- cbind(coords,x.dat$temp)
## heatload.x <- cbind(coords,x.dat$heatload)
## cti.x <- cbind(coords,x.dat$cti)

## colors <- brewer.pal(5, "Spectral")
## pal <- colorRampPalette(colors)


## mba.slope <- mba.surf(slope.x, 300,300, extend=TRUE)$xyz.est
## image(mba.slope, xaxs="r", yaxs="r", col=pal(50))
## plot(AK.crop, axes=TRUE, add=TRUE)

## mba.temp <- mba.surf(temp.x, 300,300, extend=TRUE)$xyz.est
## image(mba.temp, xaxs="r", yaxs="r", col=pal(50))
## plot(AK.crop, axes=TRUE, add=TRUE)


## mba.heat <- mba.surf(heatload.x, 30,30, extend=TRUE)$xyz.est
## image(mba.heat, xaxs="r", yaxs="r", col=pal(50))
## plot(AK.crop, axes=TRUE, add=TRUE)

## mba.cti <- mba.surf(cti.x, 30,30, extend=TRUE)$xyz.est
## image(mba.cti, xaxs="r", yaxs="r", col=pal(50))
## lot(AK.crop, axes=TRUE, add=TRUE)

## write data
write.table(y[ho], "./raw/y.ho", row.names=F, col.names=F, sep="\t")
write.table(x[ho,], "./raw/x.ho", row.names=F, col.names=F, sep="\t")
write.table(coords[ho,], "./raw/coords.ho", row.names=F, col.names=F, sep="\t")

write.table(y[-ho], "./raw/y.mod", row.names=F, col.names=F, sep="\t")
write.table(x[-ho,], "./raw/x.mod", row.names=F, col.names=F, sep="\t")
write.table(coords[-ho,], "./raw/coords.mod", row.names=F, col.names=F, sep="\t")




####################
## Plot some BINS ##
####################

covs <- cbind(y,x)
slopes <- covs$slope
temps <- covs$temp
heat <- covs$heatload
ctis <- covs$cti
pf <- covs$y

#binned.slopes=cut(slopes, breaks=c(0.00, 0.07, 0.14, 0.21, 0.28, 0.35, 0.42, 0.49, 0.56, 0.63, 0.7,
#                                   1.02, 1.34, 1.66, 1.98, 2.30, 2.62, 2.94, 3.26, 3.58, 3.90,
#                                   6.075,  8.250, 10.425, 12.600, 14.775, 16.950, 19.125, 21.300, 23.475, 25.650, 27.825, 30, Inf), labels=c(as.character(0:31),'30+'))

binned.slopes=cut(slopes, breaks=c(0:35), labels=c(as.character(1:35)))
binned.temps=cut(temps, breaks=c(-15:15), labels=c(as.character(-15:14)))
binned.heatload=cut(heat, breaks=c(seq(0,200,10)), labels=c(as.character(seq(10,200,10))))
binned.cti=cut(ctis, breaks=c(seq(0,105,5)), labels=c(as.character(0:20)))
covs=cbind(covs,binned.slopes, binned.temps, binned.heatload, binned.cti)

heights <- tapply(covs$y,binned.temps,mean)

# Change Heatload to HLI
png("png/tempbins.png", width=585, height=585, family="Courier")
#postscript('png/slopebins.eps', width=8, height=8,horizontal = FALSE, onefile = FALSE, paper = "special",pointsize=18, family = "ComputerModern", encoding = "TeXtext.enc")
barplot(heights, ylim=c(0,1),
        xlim=c(-25,25),
        ylab="Probability of permafrost",
        xlab="Mean annual temperature",     
        col="lightgrey")

                                        #, xaxt="n")
 #axis(1, lwd.ticks=0,at=0:20, labels=as.character(seq(0,100,5)))
 #text(5, 0.8, "Flat")
 #text(5.1, 0.75, "0 to 0.7%")
 #text(18, 0.8, "Low Slope")
 #text(18.1, 0.75, "0.7 to 3.9%")
 #text(32.5, 0.8, "High Slope")
 #text(32.5, 0.75, "3.9% and above")
box(lty="solid")
rug(temps, side=3)
#lines(x=mat[,2],y=mat[,1],col='red')
dev.off()
