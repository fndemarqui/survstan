
functions{
#include /inst/stan/chunks/baselines.stan
#include /inst/stan/chunks/loglikelihoods.stan
}

data{
  int n;
  int p;
  vector[n] time;
  vector[n] event;
  matrix[p == 0 ? 0 : n, p] X;
  real tau;
  int baseline;
  int survreg;    // 1 - AFT; 2 - PH; 3 - PO; 4 - AH
}

transformed data{
  int is_alpha = 0;
  int is_gamma = 0;
  int is_lambda = 0;
  int is_mu = 0;
  int is_sigma = 0;

  vector[p == 0 ? n : 0] zeros;

  if(p==0){
    for(i in 1:n){
      zeros[i] = 0;
    }
  }

  if(baseline == 1){ //exponential
    is_lambda = 1;
  }else if(baseline == 2){ // weibull
    is_alpha = 1;
    is_gamma = 1;
  }else if(baseline == 3){ // lognormal
    is_mu = 1;
    is_sigma = 1;
  }else if(baseline == 4){ // loglogistic
    is_alpha = 1;
    is_gamma = 1;
  }


}

parameters{
  vector[p == 0 ? 0 : p] beta;
  real<lower=0> alpha[is_alpha == 0 ? 0 : 1];
  real<lower=0> gamma[is_gamma == 0 ? 0 : 1];
  real<lower=0> lambda[is_lambda == 0 ? 0 : 1];
  real mu[is_mu == 0 ? 0 : 1];
  real<lower=0> sigma[is_sigma ==  0 ? 0 : 1];
}


model{

  vector[n] lp;
  vector[n] y;
  vector[n] loglik;
  vector[n] lpdf;
  vector[n] lsurv;

  if(p>0){
    lp = X*beta;
  }else{
    lp = zeros;
  }

  if(survreg == 1){
    y = time ./ exp(lp);
  }else if(survreg == 4){
    y = time .* exp(lp);
  }else{
    y = time;
  }



  if(baseline == 1){ // exponential
    for(i in 1:n){
      lpdf[i] = exponential_lpdf(y[i]|lambda);
      lsurv[i] = exponential_lccdf(y[i]|lambda);
    }
  }else if(baseline == 2){ // Weibull
    for(i in 1:n){
          lpdf[i] = weibull_lpdf(y[i]|alpha, gamma);
          lsurv[i] = weibull_lccdf(y[i]|alpha, gamma);
    }
  }else if(baseline == 3){ // lognormal
    for(i in 1:n){
      lpdf[i] = lognormal_lpdf(y[i]|mu, sigma);
      lsurv[i] = lognormal_lccdf(y[i]|mu, sigma);
    }
  }else{
    for(i in 1:n){
      lpdf[i] = loglogistic_lpdf(y[i]|alpha[1], gamma[1]);
      lsurv[i] = loglogistic_lccdf(y[i]|alpha[1], gamma[1]);
    }
  }

  if(survreg ==  1){ //AFT model
    loglik = loglik_aft(lpdf, lsurv, event, lp, tau);
  }else if(survreg == 2){ //PH model
    loglik = loglik_ph(lpdf, lsurv, event, lp, tau);
  }else if(survreg == 3){ //PO model
    loglik = loglik_po(lpdf, lsurv, event, lp, tau);
  }else{ //AH model
    loglik = loglik_ah(lpdf, lsurv, event, lp, tau);
  }

  target += sum(loglik);

}



