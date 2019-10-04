#----------------------------------------------------------------------------#
# Date: 3 October 2019
# Description code: Stan program for Bayesian inference on visit length at the 
#  Feeder as the next visit is less or equal than 60 seconds, the dataset
#  belong to the trial 1, the model include fixed effects: location and weight 
#---------------------------------------------------------------------------#
setwd("C:/Users/marti/OneDrive/Documents/job/bayesian inference/Minicourse/")
rm(list = ls())

library(tidyverse)
library(coda)
library(mcmcplots)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load data and filtered by trial 1
#------

load("trial1_wtcvoar.Rdata")
# Filter next visit less or equal than 60 sc, delete the 7 first days
tn<-60
trial1.event<-trial1.event%>%filter(Consumed > 0 & to_next<=tn & trial.day>=7 )
dim(trial1.event)
#-------------------------------------------------------------------------------------
# 1. Model with Eartag as random effect
#-------------------------------------------------------------------------------------
# Desing matrix with fixed effects
Xmatrix = model.matrix(~trial1.event$Location + trial1.event$wt_md -1)
# Eartag vector
Eartag<-as.numeric(as.factor(trial1.event$Ear_Tag))
# Stan data file
stanDat1 <- list(Nobs = nrow(trial1.event), 
                 Npreds=ncol(Xmatrix),
                 Neart=length(unique(Eartag)),
                 ET=Eartag,
                 x = Xmatrix,
                 y = trial1.event$visit.length)
#--------------------------------------------*
# 1.1. Sample from posterior distribution

system.time(
  {m1t1 <- stan(file =  "M1_trial1.stan", data = stanDat1,
                pars = c("beta", "var_eartag", "var_error", "prp_var_eartag","prp_var_error"),
                save_warmup=FALSE,iter = 6000, chains = 1, warmup = 1000)})

save(m1t1, file="M1_trial1_Stan.Rdata")


#----------------------------------------------------------------------------------
# 2. Model with Eartag plus Follower random effects
#----------------------------------------------------------------------------------

#  Followers vector
follower<-as.numeric(as.factor(trial1.event$follower))

# Stan data file
stanDat2 <- list(Nobs = nrow(trial1.event), 
                 Npreds=ncol(Xmatrix),
                 Neart=length(unique(Eartag)),
                 Nfoll=length(unique(follower)),
                 ET=Eartag,
                 Foll=follower,
                 x = Xmatrix,
                 y = trial1.event$visit.length)

#--------------------------------------------*
# 2.1. Sample from posterior distribution
#-------- 
system.time({
  m2t1<- stan(file = "M2_trial1.stan", data = stanDat2,
              pars = c("beta", "var_eartag", "var_follower","var_error", 
                       "prp_var_eartag","prp_var_follower","prp_var_error"),
              save_warmup=FALSE,iter = 12000, chains = 1, warmup = 2000)})


save(m2t1, file = "M2_trial1_Stan.Rdata")


#----------------------------------------------------------------------------------
# 3. Model with covariance between eartag and follower random effects
#----------------------------------------------------------------------------------
# Stan data file
stanDat3 <- list(Nobs = nrow(trial1.event), 
                 Npreds=ncol(Xmatrix),
                 Neart=length(unique(Eartag)),
                 ET=Eartag,
                 Foll=follower,
                 x = Xmatrix,
                 y = trial1.event$visit.length)

#--------------------------------------------*
# 3. 1. Sampling posterior distribution

system.time({
  m3t1 <- stan(file =  "M3_trial1.stan", data = stanDat3,
              pars = c("beta","rho", "var_eartag", "var_follower", "var_error",
                       "prp_var_eartag","prp_var_follower","prp_var_error"),
              save_warmup=FALSE,iter = 12000, chains = 1, warmup = 2000)})

save(m3t1, file = "M3_trial1_Stan.Rdata")

sessionInfo()


