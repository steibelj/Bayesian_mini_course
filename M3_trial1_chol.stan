//
// This Stan program mixed model effects to estimating variance components on 
// visit length time at the feeder. 'y' is modeled as normally distributed, with 
// fixed effects location, and ,median weight, and random effects od Eartag and follower,
// and the correaltion parameter between Eartag and follower effects
// 

// 1. The input data is a vector 'y' of length 'Nobs', locations and median weight,
// Number of individuasl and.
data {
  int<lower=0> Nobs;                     // observation number
  int<lower=0> Npreds;                  // locations and wt covariate predictors 
  int<lower=0> Neart;                   // eartags random effects number
  int<lower=1, upper=Neart> ET[Nobs];  // eartag identifier records
  int<lower=1,upper=Neart> Foll[Nobs]; // follower identifier records
  matrix[Nobs,Npreds] x;               // desing matrix fixed effects
  vector[Nobs] y;                     // visit lenght vector data
  
}

transformed data{

  matrix[Nobs,Npreds] Q_ast;
  matrix[Npreds,Npreds] R_ast;
  matrix[Npreds,Npreds] R_ast_inv;
  vector[2] mu;
  
// thin and scale the QR decomposition
  
  Q_ast = qr_Q(x)[, 1:Npreds]*sqrt(Nobs - 1);
  R_ast = qr_R(x)[1:Npreds, ]/sqrt(Nobs - 1);
  R_ast_inv = inverse(R_ast);

 
  for(i in 1:2) mu[i]=0;
}
// 2. Model paramenters: Fixed: location, median weight,
// Random: sd error, sd Eartag, sd Follower, correlation EF
parameters{
  vector[Npreds] theta;              // fixed effects vector
  matrix[Neart,2] Animal;           // matrix
  vector<lower=0.1, upper=10>[2] sigma_uf;      // sd eartag effect, sd follower vector
  real<lower=0> sigma_e;            // standar deviation error
  real<lower=-1,upper=1> rho;       // correlation(Eartag,Follower) 
  }


// 3. Model specification: Priors to parameters and likelihood function
model{
  
  matrix[2,2] Sigma;
  
  Sigma[1,1]= pow(sigma_uf[1],2);
  Sigma[2,2]= pow(sigma_uf[2],2);
  Sigma[1,2]= (sigma_uf[1]*sigma_uf[2]*rho);
  Sigma[2,1]= Sigma[1,2];
  
  // Priors
  for (i in 1:Neart)
  Animal[i]~multi_normal(mu,Sigma);

  for  (i in 1: Nobs)                       
  y[i]~ normal(Q_ast[i]*theta + Animal[ET[i],1] + Animal[Foll[i],2], sigma_e);  
  
}
generated quantities{
  // generate estimate variance components of the model
 vector[Npreds] beta;
 real var_eartag;           //  eartag variance
 real var_follower;          // follower variance
 real var_error;               // error variance
 real prp_var_eartag;        // proportion eartag variance
 real prp_var_follower;     // proportion follower variance
 real prp_var_error;       // proportion error variance
 
 beta = R_ast_inv *theta;
 var_eartag = pow(sigma_uf[1],2);  ///
 var_follower = pow(sigma_uf[2],2);
 var_error = pow(sigma_e,2);

 prp_var_eartag = var_eartag/(var_eartag + var_follower + var_error);  // compute the proportion of variance eartag
 prp_var_follower = var_follower/(var_eartag + var_follower + var_error);  // compute the proportion of variance eartag
 prp_var_error = var_error/(var_eartag + var_follower + var_error);   // compute the proportion of variance error
 
}


