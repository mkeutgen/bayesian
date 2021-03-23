### CHAPTER 2: FROM PRIOR TO POSTERIOR DISTRIBUTION ###

# Disclaimer : not my own code
#############################################################
## BINOMIAL CASE
#############################################################

## data specification
s<-9
n<-12

## prior specification
alpha<-2
beta<-2


# METHOD 1: using bayes theorem, numerically calculating denominator (AUC)

library(DescTools)

theta<-seq(0,1,0.01)

likelihood<-dbinom(s,n,theta)
prior<-dbeta(theta,alpha,beta)
#prior<-dbeta(theta,alpha,beta)*(theta>0.5)*2
posterior<-likelihood*prior/AUC(theta,likelihood*prior,method="trapezoid")

plot(theta,likelihood,type="l",ylim=c(0,5),col="blue",ylab="")
lines(theta,likelihood/AUC(theta,likelihood,method="trapezoid"),col="blue")
lines(theta,prior,col="red")
lines(theta,posterior)


# METHOD 2: using analytical derivations for binomial setting

theta<-seq(0,1,0.01)

likelihood<-dbeta(theta,s+1,n-s+1)
prior<-dbeta(theta,alpha,beta)
posterior<-dbeta(theta,alpha+s,beta+n-s)

plot(theta,likelihood,type="l",ylim=c(0,5),col="blue",ylab="")
lines(theta,prior,col="red")
lines(theta,posterior)



#############################################################
## GAUSSIAN CASE
#############################################################

## data specification
n<-50       # sample size
ybar<-318     # sample mean
sd<-119.5      # sample sd

## prior specification
mu0<-328     # prior mean
sd0<-5       # prior sd


# METHOD 1: using bayes theorem, numerically calculating denominator (AUC)

library(DescTools)

mu<-seq(250,400,0.01)

likelihood<-dnorm(mu,mean=ybar,sd=sd/sqrt(n))
prior<-dnorm(mu,mean=mu0,sd=sd0)
posterior<-likelihood*prior/AUC(mu,likelihood*prior,method="trapezoid")

plot(mu,likelihood,ylim=c(0,0.1),type="l",col="blue",ylab="")
lines(mu,prior,col="red")
lines(mu,posterior)


# METHOD 2: using analytical derivations for gaussian setting

mu<-seq(250,400,0.01)

likelihood<-dnorm(mu,mean=ybar,sd=sd/sqrt(n))
prior<-dnorm(mu,mean=mu0,sd=sd0)
prec<-1/sd0^2+n/sd^2
varbar<-1/prec
mubar<-(mu0/sd0^2 + ybar*n/sd^2)*varbar
posterior<-dnorm(mu,mean=mubar,sd=sqrt(varbar))

plot(mu,likelihood,type="l",ylim=c(0,0.10),col="blue",ylab="")
lines(mu,prior,col="red")
lines(mu,posterior)




#############################################################
## POISSON CASE
#############################################################

## data specification
n<-4351       # sample size
sum.y<-9758     # sum of the counts

## prior specification
alpha0<-3     # prior shape parameter
beta0<-1       # prior rate parameter


# METHOD 1: using bayes theorem, numerically calculating denominator (AUC)

library(DescTools)

mu<-seq(0,10,0.001)

likelihood<-dgamma(mu,shape=sum.y+1,rate=n)
prior<-dgamma(mu,shape=alpha0,rate=beta0)
posterior<-likelihood*prior/AUC(mu,likelihood*prior,method="trapezoid")

plot(mu,likelihood,type="l",col="blue",ylab="",ylim=c(0,20))
lines(mu,prior,col="red")
lines(mu,posterior)

# METHOD 2: using analytical derivations for poisson setting

mu<-seq(0,10,0.001)

likelihood<-dgamma(mu,shape=sum.y+1,rate=n)
prior<-dgamma(mu,shape=alpha0,rate=beta0)
posterior<-dgamma(mu,shape=alpha0+sum.y,rate=beta0+n)

plot(mu,likelihood,type="l",col="blue",ylab="",ylim=c(0,20))
lines(mu,prior,col="red")
lines(mu,posterior)

