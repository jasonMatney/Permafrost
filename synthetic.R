rm(list=ls())
library(MASS)
library(fields)
library(MBA)
library(spBayes)
set.seed(1)

## Not run: 
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

##Make some data with spatially-varying regression coefficients. We'll
##assume an intercept and two slope parameters vary.
n <- 400

coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))
D <- as.matrix(dist(coords))

p <- 3

theta <- c(3/0.5, 3/0.5, 3/0.5)

sigma.sq <- c(0.1,0.1,0.1)

C.1 <- sigma.sq[1]*exp(-theta[1]*D)
C.2 <- sigma.sq[2]*exp(-theta[2]*D)
C.3 <- sigma.sq[3]*exp(-theta[3]*D)

w.1 <- rmvn(1, rep(0, n), C.1)
w.2 <- rmvn(1, rep(0, n), C.2)
w.3 <- rmvn(1, rep(0, n), C.3)

B <- c(0.001, 0.001, -0.001)

x <- cbind(1, rnorm(n), rnorm(n))

mu <- x[,1]*(B[1]+w.1) + x[,2]*(B[2]+w.2) + x[,3]*(B[3]+w.3)

p <- 1/(1+exp(-mu))

y <- rbinom(n, size=1, prob=p)

ho <- sample(1:n, floor(0.5*n))

write.table(y[ho], "y.ho", row.names=F, col.names=F, sep="\t")
write.table(x[ho,], "x.ho", row.names=F, col.names=F, sep="\t")
write.table(cbind(w.1,w.2,w.3)[ho,], "w.ho", row.names=F, col.names=F, sep="\t")
write.table(coords[ho,], "coords.ho", row.names=F, col.names=F, sep="\t")

write.table(y[-ho], "y.mod", row.names=F, col.names=F, sep="\t")
write.table(x[-ho,], "x.mod", row.names=F, col.names=F, sep="\t")
write.table(cbind(w.1,w.2,w.3)[-ho,], "w.mod", row.names=F, col.names=F, sep="\t")
write.table(coords[-ho,], "coords.mod", row.names=F, col.names=F, sep="\t")
write.table(rep(1,nrow(coords[-ho,])), "weights.mod", row.names=F, col.names=F, sep="\t")


##Make some data with spatially-varying regression coefficients. Now we'll
##assume only the intercept varies.
n <- 400

coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))
D <- as.matrix(dist(coords))

p <- 3

theta <- c(3/0.5, 3/0.5, 3/0.5)

sigma.sq <- c(0.1,0.1,0.1)

C.1 <- sigma.sq[1]*exp(-theta[1]*D)

w.1 <- rmvn(1, rep(0, n), C.1)

B <- c(0.001, 0.001, -0.001)

x <- cbind(1, rnorm(n), rnorm(n))

mu <- x[,1]*(B[1]+w.1) + x[,2]*(B[2]) + x[,3]*(B[3])

p <- 1/(1+exp(-mu))

y <- rbinom(n, size=1, prob=p)

ho <- sample(1:n, floor(0.5*n))

write.table(y[ho], "y.ho", row.names=F, col.names=F, sep="\t")
write.table(x[ho,], "x.ho", row.names=F, col.names=F, sep="\t")
write.table(cbind(w.1)[ho,], "w.ho", row.names=F, col.names=F, sep="\t")
write.table(coords[ho,], "coords.ho", row.names=F, col.names=F, sep="\t")

write.table(y[-ho], "y.mod", row.names=F, col.names=F, sep="\t")
write.table(x[-ho,], "x.mod", row.names=F, col.names=F, sep="\t")
write.table(cbind(w.1)[-ho,], "w.mod", row.names=F, col.names=F, sep="\t")
write.table(coords[-ho,], "coords.mod", row.names=F, col.names=F, sep="\t")
write.table(rep(1,nrow(coords[-ho,])), "weights.mod", row.names=F, col.names=F, sep="\t")
