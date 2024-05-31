//Standard Hierarchical Linear Model
//Notes:
data {
  int<lower=0> N;   // number of observations items through time, across regions
  int<lower=0> NP;   // number of time periods
  int<lower=0> P[N];   // pointer vector for number of periods
  int<lower=0> K;   // number of predictors (only 1 for this analysis)
  matrix[N, K] x;   // predictor matrix - i.e., SST, etc.
  real y[N];      // outcome vector - first pass this is upwelling
}
parameters {
  real alphaP[NP];  // intercepts 
  real betaP[NP];  // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
    // priors
  for(p in 1:NP) {
  alphaP[p] ~ normal(0, 10);
  betaP[p] ~ normal(0, 10);
  }
  sigma ~ normal(0, 10);
  
  //LIKELIHOOD
  for(n in 1:N){
     for(k in 1:K){
    y[n] ~ normal(x[n, k] * betaP[P[n]]+ alphaP[P[n]], sigma);  // likelihood
     }
  }
  

}

generated quantities {
  real y_new[N];  
for(n in 1:N) {
      y_new[n] = normal_rng(x[n, K] * betaP[P[n]] + alphaP[P[n]], sigma);  // likelihood

}
}





