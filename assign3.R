## Last update Sat Nov 13 00:20:32 2010 GONG-YI LIAO

require(MCMCpack)
require(arm)
require(R2jags)


## Q1, q2 on page 95, BDA textbook

ξ <- c(1,1,1)
y1 <- c(294, 307, 38)
y2 <- c(288, 332, 19)
π.1 <- rdirichlet(10000, y1+ξ)
π.2 <- rdirichlet(10000, y2+ξ)
α.1 <- π.1[,1]/rowSums(π.1[,1:2])
α.2 <- π.2[,1]/rowSums(π.2[,1:2])
pdf("1-1.pdf")
hist(α.2 - α.1, prob=TRUE, main=expression("Histogram of " * alpha[2] - alpha[1]), xlab=expression(alpha[2] - alpha[1]))
lines(density(α.2 - α.1), type="l", col="red")
dev.off()

rm(ξ, α.1, α.2, π.1, π.2, y1, y2)

## Q2, q8 on page 98, BDA textbook

## original data set
## bicycles <- c(16, 9, 10, 13, 19, 20, 18, 17, 35, 55,
##               12, 1, 2, 4, 9, 7, 9, 8, NA, NA,
##               8, 35, 31, 19, 38, 47, 44, 44, 29, 18,
##               10, 43, 5, 14, 58, 15, 0, 47, 51, 32,
##               60, 51, 58, 59, 53, 68, 68, 60, 71, 63,
##               8, 9, 6, 9, 19, 61, 31, 75, 14, 25)
##
## vehicles <- c(58, 90, 48, 57, 103, 57, 86, 112, 273, 64,
##               113, 18, 14, 44, 208, 67, 29, 154, NA, NA,
##               29, 415, 425, 42, 180, 675, 620, 437, 47, 462,
##               557, 1258, 499, 601, 1163, 700, 90, 1093, 1459, 1086,
##               1545, 1499, 1598, 503, 407, 1494, 1558, 1706, 476, 752,
##               1248, 1246, 1596, 1765, 1290, 2498, 2346, 3101, 1918, 2318)
##
## traffic.counting <- data.frame(type=c(rep("Residential", 20), rep("Fairly busy", 20), rep("Busy", 20)),
##                                bike.route=rep(c(rep(TRUE, 10), rep(FALSE, 10)), 3),
##                                bicycles=bicycles, vehicles=vehicles)
##
## write.csv(traffic.counting, 'bike.csv', row.names=FALSE)
## rm(bicycles, vehicles, traffic.counting)

bike.data <- read.csv('bike.csv', header=TRUE)


y1 <- bike.data[1:10,3]/rowSums(bike.data[1:10, 3:4])
z1 <- bike.data[11:18,3]/rowSums(bike.data[11:18, 3:4])


data.3.1.y <- list(y=y1, N=10, p1=.5, r1=10, p2=.5, r2=1000)

data.3.1.z <- list(y=z1, N=8, p1= .5, r1=10, p2=.5, r2=1000)

jags.3.1.y <- jags(data=data.3.1.y, inits=list("alpha"=rnegbin(1, 30, theta=1), "beta"=rnegbin(1, 700, theta=1)),
                   parameters.to.save=c("alpha", "beta"), model.file="2-1-1.bug",
                   n.iter=1000, n.chains=2)

jags.3.1.z <- jags(data=data.3.1.z, inits=list("alpha"=rnegbin(1, 30, theta=1), "beta"=rnegbin(1, 700, theta=1)),
                   parameters.to.save=c("alpha", "beta"), model.file="2-1-1.bug",
                   n.iter=1000, n.chains=2)


gelman.diag(as.mcmc(jags.3.1.y))
geweke.diag(as.mcmc(jags.3.1.y))
plot(as.mcmc(jags.3.1.y))


gelman.diag(as.mcmc(jags.3.1.z))
geweke.diag(as.mcmc(jags.3.1.z))
plot(as.mcmc(jags.3.1.z))


## second approach of modeling

y1.logit <- as.numeric(logit(y1))
z1.logit <- logit(z1)
  
data.3.2.y <- list(y=y1.logit, mu0=0, alpha0=3, beta0=1, tau2=2, N=10)
data.3.2.z <- list(y=z1.logit, mu0=0, alpha0=3, beta0=1, tau2=2, N=8)

inits.1 <- function() {
  list("mu"=rnorm(1), "invsgm2"=rinvgamma(1,3,1))
}

jags.3.2.y <- jags(data=data.3.2.y, inits=inits.1, 
                   parameters.to.save=c("mu", "sgm2"),
                   model.file="3-2.bug", n.iter=1000, n.chains=2)

gelman.diag(as.mcmc(jags.3.2.y))
geweke.diag(as.mcmc(jags.3.2.y))
xyplot(as.mcmc(jags.3.2.y))

jags.3.2.z <- jags(data=data.3.2.y, inits=inits.1, 
                   parameters.to.save=c("mu", "sgm2"),
                   model.file="3-2.bug", n.iter=1000, n.chains=2)

gelman.diag(as.mcmc(jags.3.2.z))
geweke.diag(as.mcmc(jags.3.2.z))
xyplot(as.mcmc(jags.3.2.z))




## Q3, q12 on

dose.fatal <- data.frame(dose=c(-0.86, -0.30, -0.05, 0.73),
                         anomials=rep(5,4),
                         death=c(0, 1, 3, 5))


## Q4, designed by Prof. Chen

## General functions
post.inte.two.sample <- function(theta, obs) {
  def1 <- 1
  ## the reason that we need to use the inefficient loop
  ## is that the function integration needs the exact form
  ## of the function, but, if we use "sum" or "prod" then the
  ## exact form of function will be invisible to
  ## the function "integration" and error will occur at the second
  ## stage.
  for (i in 1:length(obs))
    def1 <- def1/(1+(theta - obs[i])^2)
  def1
}

post.two.sample.cauchy <- function(theta, obs) {
  norm.post.0 <- integrate(post.inte.two.sample, lower=-Inf, upper=Inf, obs=obs)$value
  post.inte.two.sample(theta, obs)/norm.post.0
}

logfun1 <- function(theta, obs) {
  ## usage:
  log(post.two.sample.cauchy(theta, obs))
}

post.mean.fun.1 <- function(theta, obs) {
  theta*post.two.sample.cauchy(theta, obs)
}


## plot for (4.b.i)
pdf("4-b-i.pdf")

plot(seq(-5, 9, length=10000),
     post.two.sample.cauchy(seq(-5, 9, length=10000), obs=c(1.5, 2.5)),
     type="l", col="red", main=expression("posterior density of " * theta * "conditional on y = (1.5, 2.5)" ),
     xlab=expression(theta * "'s support"), ylab="posterior density")

dev.off()

## Metropolis-Hasting sampled θ for y = (1.5, 2.5)
post.sample.b.ii <- MCMCmetrop1R(logfun1, theta.init=rnorm(1), obs=c(1.5, 2.5),
                                 thin=1, mcmc=40000, burnin=500, logfun=TRUE) 

mean(post.sample.b.ii)
var(post.sample.b.ii)
mean.nc <- integrate(post.mean.fun.1, lower=-Inf, upper=Inf, obs=c(1.5, 2.5))$value

## plot for (4.b.ii)
pdf("4-b-ii.pdf")
plot(seq(-5, 9, length=10000),
     post.two.sample.cauchy(seq(-5, 9, length=10000), obs=c(1.5, 2.5)),
     type="l", col="red", main=expression("posterior density of " * theta * "conditional on y = (1.5, 2.5)" ),
     xlab=expression(theta * "'s support"), ylab="posterior density", lty=5)

lines(density(post.sample.b.ii), col="blue", lty=2)
dev.off()


## plot (4.c.i)
pdf("4-c-i.pdf")

plot(seq(-9, 13, length=20000),
     post.two.sample.cauchy(seq(-9, 13, length=20000), obs=c(-3, -2, 1.5, 2.5)),
     type="l", col="red", main=expression("posterior density of " * theta * "conditional on y = (-3, -2, 1.5, 2.5)"),
     xlab=expression(theta * "'s support"), ylab="posterior density")

dev.off()

## Metropolis-Hasting sampled θ for y = (-3, -2, 1.5, 2.5)
post.sample.c.ii <- MCMCmetrop1R(logfun1, theta.init=rnorm(1), obs=c(-3, -2, 1.5, 2.5),
                                 thin=1, mcmc=80000, burnin=500, logfun=TRUE) 

mean(post.sample.c.ii)
var(post.sample.c.ii)
mean.nc <- integrate(post.mean.fun.1, lower=-Inf, upper=Inf, obs=c(-3, -2, 1.5, 2.5))$value

## plot for (4.b.ii)
pdf("4-c-ii.pdf")

plot(seq(-9, 13, length=20000),
     post.two.sample.cauchy(seq(-9, 13, length=20000), obs=c(-3, -2, 1.5, 2.5)),
     type="l", col="red", main=expression("posterior density of " * theta * "conditional on y = (-3, -2, 1.5, 2.5)" ),
     xlab=expression(theta * "'s support"), ylab="posterior density", lty=5)

lines(density(post.sample.c.ii), col="blue", lty=2)
dev.off()

## Q5, q11 on page 155, BDA textbook.

y.Js <- bike.data[1:10,3]
n.Js <- rowSums(bike.data[1:10, 3:4])

αβ.lkhd <- function(α, β, ys=y.Js, ns=n.Js) {
  ## usage:
  ref.svl <- 0
  for (j in 1:length(ys))
    ref.svl <-ref.svl + ns[j]*log(2) + dbinom(ys[j], size=ns[j], prob=.5, log=TRUE)
                      + lgamma(ys[j]+α) + lgamma(ns[j]-ys[j]+β) - lgamma(ns[j]+α+β)
  ref.svl <- ref.svl + (lgamma(α+β)-lgamma(α)-lgamma(β))*length(ys)
  ref.svl <- ref.svl - 2.5*log(α+β)
  exp(ref.svl)
}




α.1 <- seq(.1, 10, length=1000)
β.1 <- seq(.1, 30, length=1000)
αβ.lkhd.val <- outer(α.1, β.1, "αβ.lkhd")

pdf("5-1.pdf") 
image(α.1, β.1, αβ.lkhd.val, col=heat.colors(10))
contour(α.1, β.1, αβ.lkhd.val, col="brown", levels=seq(8*1E-20, 2*1E-19, by=1E-20), add=TRUE)
dev.off()


α.lkhd.marg <- function(α, β.seq, ys=y.Js, ns=n.Js) {
  slice <- (max(β.seq) - min(β.seq))/length(β.seq)
  for (i in 1:length(β.seq))
    chop <- αβ.lkhd(α,β.seq[i], ys=ys, ns=ns)
  chop <- slice*chop
  chop
}

β.lkhd.marg <- function(α.seq, β, ys=y.Js, ns=n.Js) {
  slice <- (max(α.seq) - min(α.seq))/length(α.seq)
  for (i in 1:length(α.seq))
    chop <- αβ.lkhd(α.seq[i], β, ys=ys, ns=ns)
  chop <- slice*chop
  chop
}


mar.α.numerical <- sapply(seq(.1, 12, length=200), "α.lkhd.marg" , β.seq=seq(.1, 24, length=1000))
plot(seq(.1, 12, length=200), mar.α.numerical, type="l", col="blue", xlab=expression(alpha), ylab="marginal density")

mar.β.numerical <- sapply(seq(.1, 24, length=800), "β.lkhd.marg" , α.seq=seq(.1, 6, length=200))
plot(seq(.1, 30, length=400), mar.β.numerical, type="l", col="red", xlab=expression(beta), ylab="marginal density")










