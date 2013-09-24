####
rm(list=ls())
options(width = 200)
library(sp)
library(rgdal)
library(spBayes)
library(MBA)
library(geoR)
source("/home/jason/Documents/pf/functions/functions.R")

## Holdout Data ##
y.hold <- as.vector(read.csv("model_data/y.hold.csv", header=TRUE))
x.hold <- as.matrix(read.csv("model_data/x.hold.csv", header=TRUE))
coords.hold <- as.matrix(read.csv("model_data/coords.hold.csv", header=TRUE))
## Model Data ##
y.mod <- as.vector(read.csv("model_data/y.mod.csv", header=TRUE))
x.mod <- as.data.frame(read.csv("model_data/x.mod.csv", header=TRUE))
coords.mod <- as.matrix(read.csv("model_data/coords.mod.csv", header=TRUE))

#####################
## binomial  model ##
#####################
X <-  as.matrix(cbind(x.mod["landcover"],
                      x.mod["texture"],
                      x.mod["maxndvi"],
                      x.mod["maxndvi.vi"],
                      x.mod["flowDirect"],
                      x.mod["flowAcc"],
                      x.mod["curvature"],
                      x.mod["cti"],
                     ## x.mod["smartslope"],
                      x.mod["slpperc1000"],
                      x.mod["ak_tavg_14"],
                      x.mod["heatload_1000"],
                      x.mod["aws0100wta"]
                      ))

coords <- coords.mod[,1:2]
Y <- as.vector(y.mod[,"pfscore"])

mod <- glm(Y ~ X, family = binomial("logit"))
summary(mod)
vario <- variog(coords = coords,
                data = resid(mod),
                max.dist=max(iDist(coords))/2)
variomod <- variofit(vario, cov.model="exponential")
plot(vario);lines(variomod, col = "red")

tau.sq <- variomod$nugget
sigma.sq <- variomod$cov.pars[1]
phi <- 3/variomod$cov.pars[2]

#############################
## spatial starting values ##
#############################
n.samples <- 100000
knots <- c(10,10)
beta.starting <- coefficients(mod)
beta.tuning <- t(chol(vcov(mod))) 
starting <- list("beta"=beta.starting,
                 "phi"=phi,"sigma.sq"=sigma.sq, "tau.sq"=tau.sq, "w"=0)
tuning <- list("beta"=beta.tuning, "phi"=0, "sigma.sq"=0, "tau.sq"=0, "w"=0.001)
priors <- list("beta.Flat", "phi.Unif"=c(3/max(iDist(coords)), 3/1),
               "sigma.sq.IG"=c(2,sigma.sq),
               "tau.sq.IG"=c(2,tau.sq))

cov.model <- "exponential"
verbose <- TRUE
n.report <- 10

######################
## starting values  ##
######################
set.seed(1)
permafrost.starting <- spGLM(Y~X, family="binomial", coords=coords,
                             knots=knots, starting=starting, tuning=tuning,
                             priors=priors, cov.model=cov.model,
                             n.samples=n.samples, verbose=verbose,
                             n.report=n.report)

burn.in <- as.integer(0.75 * n.samples)
w.starting <- rowMeans(permafrost.starting$p.w.knots.samples[, burn.in:n.samples])

starting <- list("beta"=beta.starting,
                 "phi"=phi,"sigma.sq"=sigma.sq, "tau.sq"=tau.sq, "w"=w.starting)
tuning <- list("beta"=beta.tuning, "phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01, "w"=0.01)

m.1 <- spGLM(Y~X, family="binomial", coords=coords,
             knots=knots, starting=starting, tuning=tuning,
             priors=priors, cov.model=cov.model,
             n.samples=n.samples, verbose=verbose,
             n.report=n.report)

summary(mcmc(m.1$p.beta.theta.samples))$quantiles

pdf("trace_collections/multitude_of_covariates.pdf")
plot(m.1$p.beta.theta.samples, density=FALSE)
dev.off()

## pdf("trace_collections/landcover.pdf")
## dev.off()
##
colMeans(m.1$p.beta.theta.samples)
##
burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples
print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

