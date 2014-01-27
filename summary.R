rm(list=ls())
## library(maptools)
library(geoR)
library(spBayes)
library(Matrix)
library(MBA)
library(fields)
#library(data.table)
library(classInt)
library(RColorBrewer)

y.mod <- as.matrix(read.table("synthetic/y.mod"))[,1]
x.mod <- as.matrix(read.table("synthetic/x.mod"))
w.mod <- as.matrix(read.table("synthetic/w.mod"))
coords <- as.matrix(read.table("synthetic/coords.mod"))

y.ho <- as.matrix(read.table("synthetic/y.ho"))[,1]
x.ho <- as.matrix(read.table("synthetic/x.ho"))
coords.ho <- as.matrix(read.table("synthetic/coords.ho"))

samps <- t(matrix(scan("samples-chain1"), nrow=9, byrow=TRUE))
nms <- c("beta.1","beta.2","beta.3",
         "phi.1","phi.2","phi.3",
         "sigma.1","sigma.2","sigma.3")

colnames(samps) <- nms

plot(mcmc(samps), density=FALSE)

##after thinning
y <- matrix(scan("samples-chain1-ySamples"), nrow=200, byrow=TRUE)
w <- matrix(scan("samples-chain1-wSamples"), nrow=600, byrow=TRUE)


n.samples <- nrow(samps)
burn.in <- floor(0.75*n.samples)
w.hat <- apply(w[,burn.in:n.samples], 1, mean)
w.1.hat <- w.hat[seq(1,length(w.hat),3)]
w.2.hat <- w.hat[seq(2,length(w.hat),3)]
w.3.hat <- w.hat[seq(3,length(w.hat),3)]

par(mfrow=c(1,3))
plot(w.1.hat, w.mod[,1])
plot(w.2.hat, w.mod[,2])
plot(w.3.hat, w.mod[,3])

y.hat <- apply(y[,burn.in:n.samples], 1, mean)

surf <- mba.surf(cbind(coords, y.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf)
points(coords)

surf <- mba.surf(cbind(coords, y.mod), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf)
points(coords)

1
## quant <- function(x){
##   quantile(x, prob=c(0.5,0.05,0.975))
## }

## rmvn <- function(n, mu=0, V = matrix(1)){
##   p <- length(mu)
##   if(any(is.na(match(dim(V),p))))
##     stop("Dimension problem!")
##   D <- chol(V)
##   t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
## }

## pred <- function(beta, tau.sq, phi, sigma.sq, w, x.ho, coords, coords.ho){

##   p <- ncol(x.ho)
  
##   n.samples <- nrow(phi)
##   n.pred <- nrow(coords.ho)
##   n <- nrow(coords)
  
##   D <- iDist(coords)
##   dD <- iDist(coords.ho, coords)
##   dd <- iDist(coords.ho)

##   w.pred <- matrix(0, n.pred*p, n.samples)
##   y.pred <- matrix(0, n.pred, n.samples)

 
##   for(s in 1:n.samples){
    
##     for(i in 1:p){

##       C <- sigma.sq[s,i]*exp(-phi[s,i]*D)
##       cC <- sigma.sq[s,i]*exp(-phi[s,i]*dD)
##       cc <- sigma.sq[s,i]*exp(-phi[s,i]*dd)
      
##       C.inv <- chol2inv(chol(C))

##       a <- cC%*%C.inv
##       mu <- a%*%w[s, seq(i, n*p, p)]
##       S <- cc-a%*%t(cC)
##       w.pred[seq(i, n.pred*p, p), s] <- rmvn(1, mu, S)
      
##     }

##     mu.y <- apply(t(matrix(rep(beta[s,], n.pred) + w.pred[,s], nrow=p))*x.ho, 1, sum)
    
##     y.pred[,s] <- rnorm(n.pred, mu.y, sqrt(tau.sq[s]))
##     print(s)
##   }
  
##   list("y.pred"=y.pred, "w.pred"=w.pred)
  
## }

## sub <- sample(floor(0.75*nrow(samps)):nrow(samps), 100)

## beta <- samps[sub,grep("beta", nms)]
## phi <- samps[sub,grep("phi", nms)]
## sigma.sq <- samps[sub,grep("sigma", nms)]
## tau.sq <- samps[sub,grep("tau", nms)]
## w.sub <- t(w[,sub])

## out <- pred(beta, tau.sq, phi, sigma.sq, w.sub, x.ho, coords, coords.ho)

## y.mu <- apply(out$y.pred, 1, quant)

## plot(y.ho, y.mu[1,])

## ## arrows(y.ho, y.mu[1,], y.ho, y.mu[3,], length=0.02, angle=90)
## ## arrows(y.ho, y.mu[1,], y.ho, y.mu[2,], length=0.02, angle=90)

## ## b <- summary(mcmc(beta))$quantiles[,3]

## ## plot(y.mod, fitted(lm(y.mod ~ x.mod-1)))

## mean((y.mu[1,]-y.ho)^2)

## ## y.hat <- apply(y[,sub], 1, mean)
## ## plot(y.mod, y.hat)
