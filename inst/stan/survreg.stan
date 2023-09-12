
functions{
#include /chunks/baselines.stan
#include /chunks/loglikelihoods.stan
}

data{
  int n;
  int p;
  vector[n] time;
  vector[n] event;
  matrix[p == 0 ? 0 : n, p] X;
  real tau;
  int baseline;
  int survreg;    // 1 - AFT; 2 - PH; 3 - PO; 4 - AH; 5 - YP
}

transformed data{
  int is_alpha = 0;
  int is_gamma = 0;
  int is_lambda = 0;
  int is_mu = 0;
  int is_sigma = 0;
  int is_phi = 0;

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
  }else if(baseline == 5){ // Birnbaum–Saunders (fatigue)
    is_alpha = 1;
    is_gamma = 1;
  }

  if(survreg == 5){
    is_phi = 1;
  }


}

parameters{
  vector[p == 0 ? 0 : p] beta;
  vector[is_phi == 0 ? 0 : p] phi;
  array[is_alpha == 0 ? 0 : 1] real<lower=0> alpha;
  array[is_gamma == 0 ? 0 : 1] real<lower=0> gamma;
  array[is_lambda == 0 ? 0 : 1] real<lower=0> lambda;
  array[is_mu == 0 ? 0 : 1] real mu;
  array[is_sigma ==  0 ? 0 : 1] real<lower=0> sigma;
}


model{

  vector[n] lp;
  vector[n] y;
  vector[n] loglik;
  vector[n] lpdf;
  vector[n] lsurv;
  vector[survreg == 5 ? n : 0] lp_long;
  vector[survreg == 5 ? n : 0] ratio;

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
  }else if(baseline == 4){
    for(i in 1:n){
      lpdf[i] = loglogistic2_lpdf(y[i]|alpha[1], gamma[1]);
      lsurv[i] = loglogistic2_lccdf(y[i]|alpha[1], gamma[1]);
    }
  }else if(baseline == 5){
    for(i in 1:n){
      lpdf[i] = fatigue_lpdf(y[i]|alpha[1], gamma[1]);
      lsurv[i] = fatigue_lccdf(y[i]|alpha[1], gamma[1]);
    }
  }

  if(survreg ==  1){ //AFT model
    loglik = loglik_aft(lpdf, lsurv, event, lp, tau);
  }else if(survreg == 2){ //PH model
    loglik = loglik_ph(lpdf, lsurv, event, lp, tau);
  }else if(survreg == 3){ //PO model
    loglik = loglik_po(lpdf, lsurv, event, lp, tau);
  }else if(survreg == 4){ //AH model
    loglik = loglik_ah(lpdf, lsurv, event, lp, tau);
  }else{
      if(p>0){
        lp_long = X*phi;
        ratio =  exp(X*(beta-phi));
      }else{
        lp_long = zeros;
        ratio = exp(zeros);
      }
    loglik = loglik_yp(event, lpdf, lsurv, lp, lp_long, ratio, tau);
  }

  target += sum(loglik);

}



