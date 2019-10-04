//
// This Stan program defines a mixed model effects, with a
// vector of values 'y' modeled as normally distributed
// with  and standard deviation 'sigma_e'(prior are default).
//

// The input data is a vector 'y' of length 'Nobs'.
data {
  int<lower=0> Nobs;                     // observation number
  int<lower=0> Npreds;                  // locations and wt covariate predictors 
  int<lower=0> Neart;                   // eartags random effects number
  int<lower=0> Nfoll;                   // follower random effects
  int<lower=1, upper=Neart> ET[Nobs];  // eartag identifier records
  int<lower=1,upper=Nfoll> Foll[Nobs]; // follower identifier records
  matrix[Nobs,Npreds] x;               // desing matrix fixed effects
  vector[Nobs] y;                     // visit lenght vector data
}

transformed data{
  matrix[Nobs,Npreds] Q_ast;
  matrix[Npreds,Npreds] R_ast;
  matrix[Npreds,Npreds] R_ast_inv;
  // thin and scale the QR decomposition
  
  Q_ast = qr_Q(x)[, 1:Npreds]*sqrt(Nobs - 1);
  R_ast = qr_R(x)[1:Npreds, ]/sqrt(Nobs - 1);
  R_ast_inv = inverse(R_ast);
  
}

// The parameters  of the model= beta: locations(2) and wt,eartag (24),"sigma_eartag",'sigma_e'.
parameters {
  vector[Npreds] theta;                           // fixed effects vector 
  vector[Neart] Eartag;                         // vector random effects eartag
  vector[Nfoll] Follower;                       // vector random effects follower
  real<lower=0, upper=10> sigma_eartag;         // sd eartag effet
  real<lower=0, upper=10> sigma_follower;       // sd follower effect
  real<lower=0> sigma_e;                       // standar deviation error
  }
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'locat+ wt cova'
// and standard deviation 'sigma_e', prior to sigma_e~ uniform(0,inf).
model {   
  // priors for beta(locations,weight) ~ uniform(-inf, inf), 
  // Eartag~ Normal(0,sigma_eartag), Follower ~ Normal(0,sigma_follower),
  // sigma_eartag and sigma_follower ~ unif(0,10), defaults stan
Eartag ~ normal(0,sigma_eartag); 
Follower ~ normal(0, sigma_follower);
  // likelihood
for  (i in 1: Nobs)                       
  y[i]~ normal(Q_ast[i]*theta + Eartag[ET[i]] + Follower[Foll[i]], sigma_e);  

}

generated quantities{
  vector[Npreds] beta;
 real var_eartag;           //  eartag variance
 real var_follower;          // follower variance
 real var_error;               // error variance
 real prp_var_eartag;        // proportion eartag variance
 real prp_var_follower;     // proportion follower variance
 real prp_var_error;       // proportion error variance

 beta = R_ast_inv *theta;
 var_eartag = pow(sigma_eartag,2);       // compute the eartag variance
 var_follower = pow(sigma_follower,2);   // compute the follower variance
 var_error = pow(sigma_e,2);            // compute the error variance
 
 prp_var_eartag = var_eartag/(var_eartag + var_follower + var_error);  // compute the proportion of variance eartag
 prp_var_follower = var_follower/(var_eartag + var_follower + var_error);  // compute the proportion of variance eartag
 prp_var_error = var_error/(var_eartag + var_follower + var_error);   // compute the proportion of variance error
 
}

