#first Rstan program
#point at your working directory
setwd("C:/Users/marti/OneDrive/Documents/job/bayesian inference/ONLINE")
rm(list=ls())
library(rstan)

x<-c(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0)
y<-c(1.1,1.84,2.6,4.2,5.0,6.4,8.1,9.1) 
N<-length(y)

data_reg<-list(x=x,y=y,N=N)
fit <- stan(file = 'reg.stan', data = data_reg,chains = 3,iter = 3000,warmup = 1000)

fit

class(fit)
names(fit)

print(fit)

library(mcmcplots) #automatically loads coda package
as_list<-As.mcmc.list(fit)
mcmcplot(as_list)

geweke.diag(as_list)





summary(lm(y~x))
plot(y~x)
plot(fit)
