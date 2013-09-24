rm(list=ls())
library(sp)
library(rgdal)
library(maptools)
library(stats)
library(raster)
library(spsurvey)
library(spBayes)
library(RColorBrewer)
library(rgeos)
library ("arm")
setwd('/home/jason/Documents/NDVI/scripts')
fpar <- read.csv("fpar.csv", header=TRUE)
lai  <- read.csv("lai.csv", header=TRUE)
alaska <- readShapeSpatial("tl_2008_02_county")
plot(alaska, xlim=c(-180,-130), ylim = c(50,75), col = "red")
points(lai$x, lai$y, col = "darkblue")
plot(1:length(5:ncol(fpar)),fpar[1,5:ncol(fpar)] )

############################
### populate fpar matrix ###
############################

fpar.mat <- matrix(NA, 30, 7)
colnames(fpar.mat) <- c('april','may','june','july','august','september','october')
years  <- c(1982:2011)
months <- c('apr','may','jun','jul','aug','sep','oct')

for(i in 1:length(years)) {
  for(j in 1:length(months)) {
    time <- paste(years[i],months[j], sep = "")
    columns <- grep(time, names(fpar))
    suma = sum(fpar[,columns[1]], na.rm=TRUE)
    sumb = sum(fpar[,columns[2]], na.rm=TRUE)
    fpar.mat[i,j] = (suma + sumb)
  }
}

d = data.frame('year'=1982:2011)
fpar.mat <- cbind(d, fpar.mat)
colnames(fpar.mat) <- c('year','april','may','june','july','august','september','october')

###########################
### populate lai matrix ###
###########################

lai.mat <- matrix(NA, 30, 7)
colnames(lai.mat) <- c('april','may','june','july','august','september','october')
years <- c(1982:2011)
months <- c('apr','may','jun','jul','aug','sep','oct')

for(i in 1:length(years)) {
  for(j in 1:length(months)) {
    time <- paste(years[i], months[j], sep = "")
    columns <- grep(time, names(lai))
    suma = sum(lai[,columns[1]], na.rm=TRUE)
    sumb = sum(lai[,columns[2]], na.rm=TRUE)
    lai.mat[i,j] = (suma + sumb)
  }
}

d = data.frame('year'=1982:2011)
lai.mat <- cbind(d, lai.mat)
colnames(lai.mat) <- c('year','april','may','june','july','august','september','october')

####################################
### principle component analysis ###
####################################

arc.pca1 <- princomp(fpar.mat[,2:8], scores = TRUE, cor = TRUE)
summary(arc.pca1)

plot(arc.pca1)
biplot(arc.pca1)

arc.pca1$loadings
arc.pca1$scores

lai.lm <- lm(log(october) ~ year, data = lai.mat)
summary(lai.lm)

#######################
### Read test files ###
#######################

setwd('/home/jason/Documents/NDVI/test')
filenames <- list.files(".", pattern="*.asc", full.names=TRUE)
fpar.list <- lapply(filenames, readGDAL)
names(fpar.list) <- substr(filenames, 3, 30)

fpar.over <- over(AK, fpar.list$fpar_2011sepb.asc)
colnames(fpar.over) <- c("fpar_2011sepb")
## isolate permafrost presence / absence ##
permafrost <- AK$PF_yes_no

character.to.binary <- function(x) {
  switch(x,
         "Yes <1" = 0,
         "No"  = 1,
         NA)
}

permafrost <- as.data.frame(sapply(permafrost, character.to.binary))
coordinates(permafrost) <- cbind(AK$long_site, AK$lat_site)
names(permafrost@data)  <- c('permafrost')






#############
### MODEL ###
#############

y <- as.matrix(permafrost@data$permafrost)
x <- fpar.over$fpar_2011sepb

dat <- cbind(y,x)
dat <- dat[!apply(dat, 1, function(x) any(is.na(x))),]
y <- dat[,1]
x <- dat[,2]

fit <- glm(y ~ x, family = binomial)
summary(fit)

y.hat <- fitted(fit) > 0.5

out <- cbind(y.hat, y)

100*sum(apply(out, 1, function(x) x[1]==x[2]))/length(y.hat)


#################
### sim model ###
#################
attach (lai.mat)

fall.clean <- as.data.frame(cbind(year, october, april))
n <- nrow(fall.clean)
fall.jitter.add <- runif (n, -.2, .2)

## Model fit
lm.fall <- lm (october ~ april)
display (lm.fall)
sim.fall <- sim (lm.fall)
beta.hat <- coef (lm.fall)
plot (april + fall.jitter.add, october, xlab="spring", ylab="fall", pch=20,col="gray10", main="Fitted linear model")

for (i in 1:50) {
  curve (sim.fall@coef[i,1] + sim.fall@coef[i,2]*x, lwd=.5, col="gray", add=TRUE)
}
curve (beta.hat[1] + beta.hat[2]*x, add=TRUE, col="red")

## Regression line as a function of one input variable
fit.2 <- glm (october ~ april)
plot (april, october, xlab="april", ylab="october")
curve (coef(fit.2)[1] + coef(fit.2)[2]*x, add=TRUE)

 # alternately
curve (cbind(1,x) %*% coef(fit.2), add=TRUE)

### Two fitted regression lines

 ## model with no interaction
 fit.3 <- lm (october ~ april + year)
#colors <- ifelse (mom_hs==1, "black", "grey")
plot (april, october, xlab="april", ylab="october",
  pch=20)
curve (cbind (1, 1, x) %*% coef (fit.3), add=TRUE, col="black")
curve (cbind (1, 0, x) %*% coef (fit.3), add=TRUE, col="grey")

 ## model with interaction
fit.4 <- lm (october ~ april + year + april:year)
#colors <- ifelse (mom_hs==1, "black", "gray")
plot (april + year, october, xlab="april", ylab="october",
  pch=20)
curve (cbind (1, 1, x, 1*x) %*% coef(fit.4), add=TRUE, col="black")
curve (cbind (1, 0, x, 0*x) %*% coef(fit.4), add=TRUE, col="gray")

### Displaying uncertainty in the fitted regression (Figure 3.10)
fit.2 <- lm (october ~ april)
display(fit.2)

fit.2.sim <- sim (fit.2)
plot (april, october, xlab="april", ylab="october", pch=20)
for (i in 1:10){
  curve (fit.2.sim@coef[i,1] + fit.2.sim@coef[i,2]*x, add=TRUE,col="gray")
}
curve (coef(fit.2)[1] + coef(fit.2)[2]*x, add=TRUE, col="red")

 # alternatively
plot (april, october, xlab="april", ylab="october", pch=20)
Oneline <- function (beta) {curve (beta[1]+beta[2]*x, add=TRUE, col="gray")}
apply (fit.2.sim$coef, 1, Oneline)
curve (coef(fit.2)[1] + coef(fit.2)[2]*x, add=TRUE, col="red")


### Displaying using one plot for each input variable (Figure 3.11)
fit.3 <- lm (october ~ april + year)
beta.hat <- coef (fit.3)
beta.sim <- sim (fit.3)@coef

october.jitter <- jitter(october)

par (mfrow=c(1,1))
plot (april, october, xlab="april", ylab="october", 
  pch=20)
#axis (1, c(80,100,120,140))
#axis (2, c(20,60,100,140))
for (i in 1:10){
  curve (cbind (1, mean(year), x) %*% beta.sim[i,], lwd=.5, col="gray",
    add=TRUE)
}
curve (cbind (1, mean(year), x) %*% beta.hat, col="black", add=TRUE)
        
plot (jitter.mom_hs, kidscore.jitter, xlab="Mother completed high school",
  ylab="Child test score", pch=20, xaxt="n", yaxt="n")
axis (1, seq(0,1))
axis (2, c(0,50,100,150))
for (i in 1:10){
  curve (cbind (1, x, mean(mom_iq)) %*% beta.sim[i,], lwd=.5, col="gray",
    add=TRUE)
}
curve (cbind (1, x, mean(mom_iq)) %*% beta.hat, col="black", add=TRUE)


