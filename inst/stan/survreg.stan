
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
  vector[n] offset;
  real tau;
  int baseline;
  int survreg;    // 1 - AFT; 2 - PH; 3 - PO; 4 - AH; 5 - YP; 6 - EH
}

transformed data{
  int is_alpha = 0;
  int is_gamma = 0;
  int is_kappa = 0;
  int is_lambda = 0;
  int is_mu = 0;
  int is_sigma = 0;
  int is_varphi = 0;
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
  }else if(baseline == 6){ // gamma
    is_alpha = 1;
    is_lambda = 1;
  }else if(baseline == 7){ // Rayleigh
    is_sigma = 1;
  }else if(baseline == 8){ // Gompertz
    is_alpha = 1;
    is_gamma = 1;
  }else if(baseline == 9){ // Generalized Gamma (Stacy)
    is_alpha = 1;
    is_gamma = 1;
    is_kappa = 1;
  }else if(baseline == 10){ // Generalized Gamma (Prentice)
    is_mu = 1;
    is_sigma = 1;
    is_varphi = 1;
  }


  //------------------------------------------------------------
  if(survreg > 4 ){
    is_phi = 1;
  }


}

parameters{
  vector[p == 0 ? 0 : p] beta;
  vector[is_phi == 0 ? 0 : p] phi;
  array[is_alpha == 0 ? 0 : 1] real<lower=0> alpha;
  array[is_gamma == 0 ? 0 : 1] real<lower=0> gamma;
  array[is_kappa ==  0 ? 0 : 1] real<lower=0> kappa;
  array[is_lambda == 0 ? 0 : 1] real<lower=0> lambda;
  array[is_mu == 0 ? 0 : 1] real mu;
  array[is_sigma ==  0 ? 0 : 1] real<lower=0> sigma;
  array[is_varphi ==  0 ? 0 : 1] real varphi;
}


model{

  vector[n] lp;
  vector[n] y;
  vector[n] loglik;
  vector[n] lpdf;
  vector[n] lsurv;
  vector[survreg > 4 ? n : 0] lp2;
  vector[survreg > 4 ? n : 0] K;

  if(p>0){
    lp = X*beta + offset;
    if(survreg > 4){
      if(survreg == 5){
        lp2 = X*phi + offset;
        K =  exp(X*(beta-phi));
      }else{
        lp2 = X*phi + offset;
        K =  exp(lp+lp2);
      }
    }
  }else{
    lp = zeros + offset;
    if(survreg > 4){
      lp2 = zeros + offset;
      K = exp(zeros);
    }
  }

  if(survreg == 2 || survreg == 3 || survreg == 5){
    y = time;
  }else{
    y = time ./ exp(lp);
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
  }else if(baseline == 6){ // Gamma
    for(i in 1:n){
          lpdf[i] = gamma_lpdf(y[i]|alpha, lambda);
          lsurv[i] = gamma_lccdf(y[i]|alpha, lambda);
    }
  }else if(baseline == 7){ // Rayleigh
    for(i in 1:n){
          lpdf[i] = rayleigh_lpdf(y[i]|sigma);
          lsurv[i] = rayleigh_lccdf(y[i]|sigma);
    }
  }else if(baseline == 8){ // Gompertz
    for(i in 1:n){
          lpdf[i] = gompertz_lpdf(y[i]|alpha[1], gamma[1]);
          lsurv[i] = gompertz_lccdf(y[i]|alpha[1], gamma[1]);
    }
  }else if(baseline == 9){ // Generalized Gamma (Stacy)
    for(i in 1:n){
          lpdf[i] = ggstacy_lpdf(y[i]|alpha[1], gamma[1], kappa[1]);
          lsurv[i] = ggstacy_lccdf(y[i]|alpha[1], gamma[1], kappa[1]);
    }
  }else if(baseline == 10){ // Generalized Gamma (Prentice)
    for(i in 1:n){
          lpdf[i] = ggprentice_lpdf(y[i]|mu[1], sigma[1], varphi[1]);
          lsurv[i] = ggprentice_lccdf(y[i]|mu[1], sigma[1], varphi[1]);
    }
  }


  if(p == 0){
    loglik = event .* lpdf + (1-event) .* lsurv;
  }else{
    if(survreg ==  1){ //AFT model
      loglik = loglik_aft(lpdf, lsurv, event, lp, tau);
    }else if(survreg == 2){ //PH model
      loglik = loglik_ph(lpdf, lsurv, event, lp, tau);
    }else if(survreg == 3){ //PO model
      loglik = loglik_po(lpdf, lsurv, event, lp, tau);
    }else if(survreg == 4){ //AH model
      loglik = loglik_ah(lpdf, lsurv, event, lp, tau);
    }else if(survreg == 5){ // YP model
      loglik = loglik_yp(event, lpdf, lsurv, lp, lp2, K, tau);
    }else{ // EH model
      loglik = loglik_eh(event, lpdf, lsurv, lp2, K, tau);
    }
  }

  target += sum(loglik);

}



