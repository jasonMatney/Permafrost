####
rm(list=ls())
options(width = 200)
library(sp)
##library(rgdal)
library(spBayes)
library(MBA)
library(geoR)
## source("/home/jason/Documents/pf/functions/functions.R")
rmvn <- function(n, mu=0, V = matrix(1)){
p <- length(mu)
if(any(is.na(match(dim(V),p))))
stop("Dimension problem!")
D <- chol(V)
t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
}
## ## Holdout Data ##
## y.hold <- as.vector(read.csv("covariates/y_hold.csv", header=TRUE))
## x.hold <- as.matrix(read.csv("covariates/x_hold.csv", header=TRUE))
## coords.hold <- as.matrix(read.csv("covariates/coords_hold.csv", header=TRUE))
## Model Data ##
y.mod <- as.vector(read.csv("covariates/y_mod.csv", header=TRUE))
x.mod <- as.data.frame(read.csv("covariates/x_mod.csv", header=TRUE))
coords.mod <- as.matrix(read.csv("covariates/coords_mod.csv", header=TRUE))
## CTI stands for composite topography index ##

#####################
## binomial  model ##
#####################
X <-  as.matrix(cbind(
                      x.mod["slope"],
                      x.mod["temp"],
                      x.mod["heatload"],
                      x.mod["cti"]
                      ))

coords <- coords.mod[,1:2]
Y <- as.vector(y.mod[,"permafrost"])

## remove duplicate coordinates
dat <- as.data.frame(cbind(X,coords))
coordinates(dat) <- dat[,5:6]
zd <- zerodist(dat)
X <- X[-zd[,2],]
coords <- coords[-zd[,2],]/1000
Y <- as.matrix(Y)
Y <- Y[-zd[,2],]

##############################
## Generalized Linear Model ##
##############################

mod <- glm(Y ~ X, family = binomial("logit"))
summary(mod)
sigma.sq <-1
phi <- 3/(0.5*max(iDist(coords)))

## collect samples
beta.starting <- coefficients(mod)
beta.tuning <-t(chol(vcov(mod)))

n.batch <- 1000
batch.length <- 50
n.samples <- n.batch*batch.length


starting <- list("beta"=beta.starting,
                 "phi"=phi,"sigma.sq"=sigma.sq, "w"=0)
tuning <- list("beta"=beta.tuning, "phi"=0.001, "sigma.sq"=0.001, "w"=0.001)
priors <- list("beta.Flat",
               "phi.Unif"=c(3/(0.75*max(iDist(coords))),
                 3/(0.01*max(iDist(coords)))),
               "sigma.sq.IG"=c(2, sigma.sq))

cov.model <- "exponential"
verbose <- TRUE
n.report <- 10
n.samples <- 50000
km <- kmeans(coords, 300)
knots <- km$centers

######################
## starting values  ##
######################
set.seed(1)
## permafrost.starting <- spGLM(Y~X, family="binomial", coords=coords,
##                             starting=starting, tuning=tuning,
##                             priors=priors, cov.model=cov.model,
##                             knots=knots,
##                              amcmc=list("n.batch"=n.batch,
##                                 "batch.length"=batch.length,
##                                 "accept.rate"=0.43),
##                             n.samples=n.samples,
##                             verbose=verbose,
##                             n.report=n.report)

## save(permafrost.starting, file="pf_start.RData")

load("pf_start.RData")

burn.in <- as.integer(0.75 * n.samples)
w.starting <- rowMeans(permafrost.starting$p.w.samples[ ,burn.in:n.samples])

starting <- list("beta"=beta.starting,
                 "phi"=phi,"sigma.sq"=sigma.sq, "w"=w.starting)
tuning <- list("beta"=c(0.001, 0.0001, 0.0001, 0.0001, 0.0001),
               "phi"=0.001,
               "sigma.sq"=0.001,
               "w"=0.001)

m.1 <- spGLM(Y~X, family="binomial", coords=coords,
              starting=starting, tuning=tuning,
              ##  amcmc=list("n.batch"=n.batch,
              ##  "batch.length"=batch.length,
              ##  "accept.rate"=0.43),   
              priors=priors, cov.model=cov.model,
              n.samples=n.samples, verbose=verbose,
              n.report=n.report)

load("pf_m1.RData")
summary(mcmc(m.1$p.beta.theta.samples))$quantiles
pdf("trace_collections/heatload_slope_cti_temp.pdf")
plot(m.1$p.beta.theta.samples)
dev.off()

burn.in <- floor(0.75*n.samples)
sub.samps <- burn.in:n.samples
print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))
beta.hat <- m.1$p.beta.theta.samples[sub.samps,
                                     c("(Intercept)",
                                       "Xslope",
                                       "Xtemp",
                                       "Xheatload",
                                       "Xcti")]
w.hat <- m.1$p.w.samples[,sub.samps]
p.hat <- 1/(1+exp(-(m.1$X%*%t(beta.hat)+w.hat)))

n <- nrow(coords)
R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- rmvn(1, rep(0,n), R)
x <- as.matrix(rep(1,n))
beta <- 0.1
p <- 1/(1+exp(-(x%*%beta+w)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n, size=1,prob=p)})
y.hat.mu <- apply(y.hat, 1, mean)
y.hat.w <- apply(w.hat, 1, mean)

##Take a look
par(mfrow=c(1,3))
surf <- mba.surf(cbind(coords,Y),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf, main="Interpolated mean of observed\n(observed rate)")
surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf, main="Interpolated mean of posterior rate\n(observed rate)")
#contour(surf, add=TRUE)
#text(coords, label=paste("(",Y,")",sep=""))
surf <- mba.surf(cbind(coords,y.hat.w),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image.plot(surf, main="Interpolated random effects of posterior rate\n(observed #
of trials)")
#contour(surf, add=TRUE)
#text(coords, label=paste("(",weights,")",sep=""))


