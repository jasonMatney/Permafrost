####
rm(list=ls())
options(width = 200)
library(sp)
library(rgdal)
library(maptools)
library(rgl)
library(foreign)
library(raster)
library(gstat)
library(rgeos)
library(spBayes)

source("/home/jason/Documents/pf/functions/functions.R")

landcover <- raster("/home/jason/Documents/pf/Alaska_data/AK_Landcover_nlcd.img")

## Alaska shapefile
AK <- getData("GADM", country = "USA", level = 1)
AK <- AK[which(AK$NAME_1 == "Alaska"),]
## ## give projection information to Alaska data ##
AK <- spTransform(AK,CRS(proj4string(landcover)))

## Get Alaska data ##
AK.pf <- read.csv("/home/jason/Documents/pf/pf-data.txt")
## reproject ##
coordinates(AK.pf) <- ~long_site+lat_site
proj4string(AK.pf) <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"
coords.pf <- coordinates(AK.pf)
coords.pf <- project(coords.pf, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0")
colnames(coords.pf) <- c('x','y')
AK.data <- slot(AK.pf,"data")
AK.pf <- cbind(AK.data,coords.pf)
coordinates(AK.pf) <- ~x+y
proj4string(AK.pf) <- proj4string(AK)
plot(AK)
points(AK.pf, col="red")
plot(landcover, add=TRUE)

################################
## read in data ##
landcover.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_landcover.csv"))
texture.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_texture.csv"))
maxndvi.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_maxndvi.csv"))
maxndvi.vi.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_maxndvi_vi.csv"))
flowDirect.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_flowDirect.csv"))
flowAcc.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_flowAcc.csv"))
curvature.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_curvature.csv"))
cti.dat <- as.matrix(read.csv("Alaska_data/covariates/AK_cti.csv"))

###################################
AK.dat <- slot(AK.pf,"data")
AK.dat <- as.data.frame(cbind(AK.dat,seq(1,nrow(AK.dat),1)))
colnames(AK.dat) <- c("contributor",
                      "profilename",
                      "pf",
                      "pfSCORE",
                      "ald",
                      "smartslope",
                      "slpperc1000",
                      "ak_tavg_14",
                      "heatload_1000",
                      "aws0100wta",
                      "X1kmPlot",
                      "X1kmPlot_5obs",
                      "nlcd_30m",
                      "nlcd_1km",
                      "aspect_1000",
                      "olthick",
                      "index")

pfscore <- AK.dat[,"pfSCORE"]

## set seed
set.seed(1)

## add index to y  ##
pfscore <- as.data.frame(cbind(pfscore,seq(1,length(pfscore),1)))
colnames(pfscore) <- c("pfscore","index")

## create holdout data for y ##
y.index <- as.data.frame(sample(as.matrix(pfscore[,2]), round(0.2 * nrow(pfscore))))
y.index <- as.data.frame(sort(y.index[,1]))
colnames(y.index) <- "index"
y.hold <- pfscore[pfscore[,2] %in% y.index[["index"]],1:2]

## add indicies to covariates ##
landcover.dat <- as.data.frame(cbind(landcover.dat,seq(1,length(landcover.dat),1)))
colnames(landcover.dat) <- c("landcover","index")
texture.dat <- as.data.frame(cbind(texture.dat,seq(1,length(texture.dat),1)))
colnames(texture.dat) <- c("texture","index")
maxndvi.dat <- as.data.frame(cbind(maxndvi.dat,seq(1,length(maxndvi.dat),1)))
colnames(maxndvi.dat) <- c("maxndvi","index")
maxndvi.vi.dat <- as.data.frame(cbind(maxndvi.vi.dat,seq(1,length(maxndvi.vi.dat),1)))
colnames(maxndvi.vi.dat) <- c("maxndvi.vi","index")
flowDirect.dat <- as.data.frame(cbind(flowDirect.dat,seq(1,length(flowDirect.dat),1)))
colnames(flowDirect.dat) <- c("flowDirect","index")
flowAcc.dat <- as.data.frame(cbind(flowAcc.dat,seq(1,length(flowAcc.dat),1)))
colnames(flowAcc.dat) <- c("flowAcc","index")
curvature.dat <- as.data.frame(cbind(curvature.dat,seq(1,length(curvature.dat),1)))
colnames(curvature.dat) <- c("curvature","index")
cti.dat <- as.data.frame(cbind(cti.dat,seq(1,length(cti.dat),1)))
colnames(cti.dat) <- c("cti","index")

#########################
## separate covariates ##
#########################
## ald: The Active Layer Depth, the depth beneath
## the soil surface that the permafrost was observed. 
ald <- as.data.frame(cbind(AK.dat["ald"],seq(1,nrow(AK.dat["ald"]),1)))
colnames(ald) <- c("ald","index")
## too much missing ald data for inclusion!
## smartslope: the slope class according to response in the permafrost probability.
## 1 - flat,0<x<0.7%, 2 - low_slope,0.7<=x<3.89%, 3 - high_slope,>=3.89%
smartslope <- as.data.frame(cbind(AK.dat["smartslope"], seq(1,nrow(AK.dat["smartslope"]),1)))
newcolumn <- matrix(NA,nrow(smartslope))
fac1 <-smartslope[1,1]
fac2 <-smartslope[4,1]
fac3 <-smartslope[5,1]
for(i in 1:nrow(smartslope)) {
    if (as.character(smartslope[i,1]) == as.character(fac1)) {
        newcolumn[i] = 1
    } else if(as.character(smartslope[i,1]) == as.character(fac2)) {
        newcolumn[i] = 2
    } else if(as.character(smartslope[i,1]) == as.character(fac3)) {
        newcolumn[i] = 3
    } else {
        newcolumn[i] = NA }
}
smartslope <- as.data.frame(cbind(AK.dat["smartslope"], newcolumn, seq(1,nrow(AK.dat["smartslope"]),1)))
colnames(smartslope) <- c("smartslope","smartslopeclass","index")
## slpperc1000: the slope percent at 1km, resampled from 45m DEM
slpperc1000 <- as.data.frame(cbind(AK.dat["slpperc1000"],seq(1,nrow(AK.dat["slpperc1000"]),1)))
colnames(slpperc1000) <- c("slpperc1000","index")
## ak_tavg_14_1km: mean annual temperature from PRISM, degrees C, snapped to topographic layers and resampled to 1km from original 2km grid
ak_tavg_14 <- as.data.frame(cbind(AK.dat["ak_tavg_14"],seq(1,nrow(AK.dat["ak_tavg_14"]),1)))
colnames(ak_tavg_14) <- c("ak_tavg_14","index")
## heatload_1000: the aspect adjusted so that north slopes are highest, south slopes lowest, and west and east slope are equal. Exp(Cosine(:aspect_1000 * (Pi() / 180)))
heatload_1000 <- as.data.frame(cbind(AK.dat["heatload_1000"],seq(1,nrow(AK.dat["heatload_1000"]),1)))
colnames(heatload_1000) <- c("heatload_1000","index")
## aws0100wta_1kmsnap: available water storage at 1km, a measure of the coarseness of the soils 
## in the parent material and how well permafrost persists.  Lower values should reflect low potential 
## for permafrost. NRCS polygons were converted to raster cells (cubic convolution) at 1km cell size 
## and snapped to topographic layers.
aws0100wta <-as.data.frame(cbind(AK.dat["aws0100wta"],seq(1,nrow(AK.dat["aws0100wta"]),1)))
colnames(aws0100wta) <- c("aws0100wta","index")
## ??
X1kmPlot <- as.data.frame(cbind(AK.dat["X1kmPlot"],seq(1,nrow(AK.dat["X1kmPlot"]),1)))
colnames(X1kmPlot) <- c("X1kmPlot","index")
##  nlcd_30m
nlcd_30m <-as.data.frame(cbind(AK.dat["nlcd_30m"],seq(1,nrow(AK.dat["nlcd_30m"]),1)))
colnames(nlcd_30m) <- c("nlcd_30m","index")
## nlcd_1km_snap: nlcd_30m resampled to 1km scale and snapped to 1km topographic variables
nlcd_1km <-as.data.frame(cbind(AK.dat["nlcd_1km"],seq(1,nrow(AK.dat["nlcd_1km"]),1)))
colnames(nlcd_1km) <- c("nlcd","index")
## aspect_1000: Slope aspect derived from 45m DEM, used to calculate 'heatload_1000'.
aspect_1000 <-as.data.frame(cbind(AK.dat["aspect_1000"],seq(1,nrow(AK.dat["aspect_1000"]),1)))
colnames(aspect_1000) <- c("aspect_1000","index")
## organic layer thickness
olthick <- as.data.frame(cbind(AK.dat["olthick"],seq(1,nrow(AK.dat["olthick"]),1)))
colnames(olthick) <- c("olthick","index")

########################################
## create holdout data for covariates ##
########################################
x.hold <- cbind(landcover.dat[  landcover.dat[,2]   %in% y.index[["index"]],1],
                texture.dat[    texture.dat[,2]     %in% y.index[["index"]],1],
                maxndvi.dat[    maxndvi.dat[,2]     %in% y.index[["index"]],1],
                maxndvi.vi.dat[ maxndvi.vi.dat[,2]  %in% y.index[["index"]],1],
                flowDirect.dat[ flowDirect.dat[,2]  %in% y.index[["index"]],1],
                flowAcc.dat[    flowAcc.dat[,2]     %in% y.index[["index"]],1],
                curvature.dat[  curvature.dat[,2]   %in% y.index[["index"]],1],
                cti.dat[        cti.dat[,2]         %in% y.index[["index"]],1],
                smartslope[     smartslope[,3]      %in% y.index[["index"]],2],
                slpperc1000[    slpperc1000[,2]     %in% y.index[["index"]],1],
                ak_tavg_14[     ak_tavg_14[,2]      %in% y.index[["index"]],1],
                heatload_1000[  heatload_1000[,2]   %in% y.index[["index"]],1],
                aws0100wta[     aws0100wta[,2]      %in% y.index[["index"]],1],
                X1kmPlot[       X1kmPlot[,2]        %in% y.index[["index"]],1],
                nlcd_30m[       nlcd_30m[,2]        %in% y.index[["index"]],1],
                nlcd_1km[       nlcd_1km[,2]        %in% y.index[["index"]],1],
                olthick[        olthick[,2]         %in% y.index[["index"]],1:2])
                              
colnames(x.hold) <- c("landcover","texture","maxndvi","maxndvi.vi","flowDirect","flowAcc","curvature","cti","smartslope","slpperc1000","ak_tavg_14","heatload_1000","aws0100wta","X1kmPlot","nlcd_30m","nlcd_1km","olthick","index")

## create holdout data for coordinates ##
coords <- coordinates(AK.pf)
coords <- as.data.frame(cbind(coords, seq(1,nrow(coords),1)))
colnames(coords) <- c("x","y","index")
coords.hold <- coords[coords[,3] %in% y.index[["index"]],] 

## create model data ##
## y.mod ##
y.mod <- as.data.frame(pfscore[!pfscore[,2] %in% y.index[["index"]],1:2])
colnames(y.mod) <- c("pfscore","index")

## x.mod ##
x.mod <- cbind(landcover.dat[ !landcover.dat     [["index"]] %in% y.index[["index"]],1],
               texture.dat[   !texture.dat       [["index"]] %in% y.index[["index"]],1],
               maxndvi.dat[   !maxndvi.dat       [["index"]] %in% y.index[["index"]],1],
               maxndvi.vi.dat[!maxndvi.vi.dat    [["index"]] %in% y.index[["index"]],1],
               flowDirect.dat[!flowDirect.dat    [["index"]] %in% y.index[["index"]],1],
               flowAcc.dat[   !flowAcc.dat       [["index"]] %in% y.index[["index"]],1],
               curvature.dat[ !curvature.dat     [["index"]] %in% y.index[["index"]],1],
               cti.dat[       !cti.dat           [["index"]] %in% y.index[["index"]],1],
               smartslope[    !smartslope        [["index"]] %in% y.index[["index"]],2],
               slpperc1000[   !slpperc1000       [["index"]] %in% y.index[["index"]],1],
               ak_tavg_14[    !ak_tavg_14        [["index"]] %in% y.index[["index"]],1],
               heatload_1000[ !heatload_1000     [["index"]] %in% y.index[["index"]],1],
               aws0100wta[    !aws0100wta        [["index"]] %in% y.index[["index"]],1],
               X1kmPlot[      !X1kmPlot          [["index"]] %in% y.index[["index"]],1],
               nlcd_30m[      !nlcd_30m          [["index"]] %in% y.index[["index"]],1],
               nlcd_1km[      !nlcd_1km          [["index"]] %in% y.index[["index"]],1],
               olthick[       !olthick           [["index"]] %in% y.index[["index"]],1:2])

colnames(x.mod) <- c("landcover","texture","maxndvi","maxndvi.vi","flowDirect","flowAcc","curvature","cti","smartslope","slpperc1000","ak_tavg_14","heatload_1000","aws0100wta","X1kmPlot","nlcd_30m","nlcd_1km","olthick","index") 

## create  model coordinates ##
coords.mod <- coords[!coords[["index"]] %in% y.index[["index"]],]

x.mod <- x.mod[-2454,]
y.mod <- y.mod[-2454,]
coords.mod <- coords.mod[-2454,]
##
write.csv(y.hold, "/home/jason/Documents/pf/model_data/y.hold.csv", row.names=FALSE)
write.csv(x.hold, "/home/jason/Documents/pf/model_data/x.hold.csv", row.names=FALSE)
write.csv(coords.hold, "/home/jason/Documents/pf/model_data/coords.hold.csv", row.names=FALSE)
write.csv(y.mod, "/home/jason/Documents/pf/model_data/y.mod.csv", row.names=FALSE)
write.csv(x.mod, "/home/jason/Documents/pf/model_data/x.mod.csv", row.names=FALSE)
write.csv(coords.mod, "/home/jason/Documents/pf/model_data/coords.mod.csv", row.names=FALSE)

## the rasters themselves ##
## sTexture <- raster("/media/Elements/Alaska_30m_data/AK_STexture_Jorgenson.img")
## maxndvi <- raster("/media/Elements/Alaska_30m_data/weld_maxndvi_annual_2008_to_2010_project.img")
## maxndvi.vi <- raster("/media/Elements/Alaska_30m_data/weld_maxndvi_annual_2008_to_2010_vi_project.img")
## flowDirect <- raster("/media/Elements/Alaska_30m_data/DEM_FlowDirection_30m_Projected.img")
## flowAcc <- raster("/media/Elements/Alaska_30m_data/DEM_FlowAcc_30m_Projected.img")
## heatload <- raster("/media/Elements/Alaska_30m_data/dem_heatload_30m_Projected.img")
## cti <- raster("/media/Elements/Alaska_30m_data/dem_cti_30m_projected.img")
## curvature <- raster("/media/Elements/Alaska_30m_data/DEM_Curvature_30m_Projected.img")
## slope <- raster("/media/Elements/Alaska_30m_data/DEM_Slope_degree_30m_Projected.img")

## needs reprojection ##
## yrb.prob <- raster("/media/Elements/Alaska_30m_data/YRB_PProbability.img")
########################

## AK.landcover <- extract(landcover,AK.pf)
## AK.texture <- extract(sTexture,AK.pf)
## AK.maxndvi <- extract(maxndvi,AK.pf)
## AK.maxndvi.vi <- extract(maxndvi.vi,AK.pf)
## AK.flow.dir <- extract(flowDirect,AK.pf)
## AK.flow.acc <- extract(flowAcc,AK.pf)
## AK.heatload <- extract(heatload,AK.pf)
## AK.cti <- extract(cti,AK.pf)
## AK.curvature <- extract(curvature,AK.pf)
## AK.slope <- extract(slope,AK.pf)
## ##
## write.csv(AK.maxndvi.vi,"AK_maxndvi_vi.csv",row.names=FALSE)
## write.csv(AK.flow.dir,"AK_flowDirect.csv",row.names=FALSE)
## write.csv(AK.flow.acc,"AK_flowAcc.csv",row.names=FALSE)
## write.csv(AK.heatload,"AK_heatload.csv",row.names=FALSE)
## write.csv(AK.cti,"AK_cti.csv",row.names=FALSE)
## write.csv(AK.curvature,"AK_curvature.csv",row.names=FALSE)
## write.csv(AK.slope,"AK_slope.csv",row.names=FALSE)


## if you want to stack ## 
## dem.rasters <- stack(dem.slope, dem.curve, dem.flow.acc, dem.flow.dir, dem.heatload, dem.cti)
## r <- landcover
## newproj <- proj4string(AK)
## x <- projectRaster(r,crs=newproj)

## permafrost data from eons ago ##
## coords <- read.table("/media/Elements/pf-data/new_data/coords")
## coords.ho <- read.table("/media/Elements/pf-data/new_data/coords.ho")
## coords.pred <- read.table("/media/Elements/pf-data/new_data/coords.pred")
## x.mod <- read.table("/media/Elements/pf-data/new_data/x")
## x.hold <- read.table("/media/Elements/pf-data/new_data/x.ho")
## x.pred <- read.table("/media/Elements/pf-data/new_data/x.pred")
## y.mod <- read.table("/media/Elements/pf-data/new_data/y")
## y.hold <- read.table("/media/Elements/pf-data/new_data/y.ho")

## dbf files :: metadata ##
## land.dbf <- read.dbf("AK_Landcover_nlcd.img.vat.dbf")
## txtr.dbf <- read.dbf("AK_STexture_Jorgenson.img.vat.dbf")
## flow.dbf <- read.dbf("DEM_FlowDirection_30m_Projected.img.vat.dbf")
## dem.dbf <- read.dbf("dem_heatload_30m_Projected.img.vat.dbf")
## yrb.dbf <- read.dbf("YRB_PProbability.img.vat.dbf")


############
## pf-data ##
#############

permafrost <- read.csv("ForAndyFeb2013.txt")

pf.coords <- read.table("pf-data/new_data/coords")
colnames(pf.coords) <- c("x","y")
pf.coords.ho <- read.table("pf-data/new_data/coords.ho")
colnames(pf.coords.ho) <- c("x","y")
pf.coords.pred <- read.table("pf-data/new_data/coords.pred")
pf.x <- read.table("pf-data/new_data/x")
colnames(pf.x) <- c("slpperc1000","ak_tavg_14","heatload_1000","aws0100wta")
pf.x.ho <- read.table("pf-data/new_data/x.ho")
colnames(pf.x.ho) <- c("slpperc1000","ak_tavg_14","heatload_1000","aws0100wta")
pf.x.pred <- read.table("pf-data/new_data/x.pred")
pf.y <- read.table("pf-data/new_data/y")
colnames(pf.y) <- c("pfscore")
pf.y.ho <- read.table("pf-data/new_data/y.ho")
colnames(pf.y.ho) <- c("pfscore")
    
dsn <- "model_data"
ogrListLayers(dsn)
ogrInfo(dsn,"akecoregions_project")
ecoregions <- readOGR(dsn=dsn, layer="akecoregions_project")

ak.tavg <- readGDAL("model_data/ak_tavg_14_1km.tif")
aspect <- readGDAL("model_data/aspect_1000.tif")
aws <- readGDAL("model_data/aws0100_1kmsnap.tif")
heatload <- readGDAL("model_data/heatload_1000.tif")
ncld <- readGDAL("model_data/ncld_1km_snap.tif")
slpperc <- readGDAL("model_data/slpperc1000.tif")
smartslope <- readGDAL("model_data/smartslope_proj1.tif")

image(ak.tavg, col = terrain.colors(20))
image(aspect, col=terrain.colors(20))
image(aws, col=terrain.colors(20))
image(heatload, col=terrain.colors(20))
image(ncld, col=terrain.colors(20))
image(slpperc, col=terrain.colors(20))
image(smartslope)
