## Last update Tue Nov  9 22:50:03 2010 GONG-YI LIAO

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

## Q4

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










