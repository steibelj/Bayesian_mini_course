//
// This Stan program defines a fixed effect model, with a
// vector of values 'y' modeled as normally distributed
// with  and standard deviation 'sigma_e'(prior are default).
//

// The input data is a vector 'y' of length 'Nobs'.
data {
  int<lower=0> Nobs;                     // observation number
  int<lower=0> Npreds;                  // locations and wt covariate predictors 
  matrix[Nobs,Npreds] x;               // desing matrix fixed effects
  vector[Nobs] y;                     // visit lenght vector data
}

// The parameters accepted by the model. Our model
// accepts parameters  beta: locations(2) and wt, 'sigma_e'.
parameters {
  vector[Npreds] beta;                  // fixed effects vector 
  real<lower=0> sigma_e;               // standar deviation error
  }

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'locat+ wt cova'
// and standard deviation 'sigma_e', prior to sigma_e~ uniform(0,inf).

model {   
  // priors for beta(locations,weight) ~ uniform(-inf, inf) 
  // likelihood
for  (i in 1: Nobs)               
  y[i]~ normal(x[i]*beta , sigma_e);  
}

generated quantities{
 real var_error;               // error variance
 
 var_error = pow(sigma_e,2);       // Compute error variance 
 
  }