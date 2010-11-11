## Last update Wed Nov 10 20:05:05 2010 GONG-YI LIAO

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









