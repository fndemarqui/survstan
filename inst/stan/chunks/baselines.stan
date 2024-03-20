

// *****************************************************************
// Loglogistic distribution:

real loglogistic2_lpdf(real x, real alpha, real gamma){
  real aux = log1p(pow(x/gamma, alpha));
  real lpdf = log(alpha) - log(gamma) + lmultiply(alpha-1, x) -
  lmultiply(alpha-1, gamma) - 2*aux;
  return lpdf;
}

real loglogistic2_lccdf(real x, real alpha, real gamma){
  real lsurv = - log1p(pow(x/gamma, alpha));
  return lsurv;
}

// *****************************************************************
// Birnbaum–Saunders (fatigue) distribution

real fatigue_lpdf(real x, real alpha, real gamma){
  real xg = sqrt(x/gamma);
  real gx = sqrt(gamma/x);
  // real xg = exp(lmultiply(0.5, x) - lmultiply(0.5, gamma));
  // real gx = exp(lmultiply(0.5, gamma) - lmultiply(0.5, x));
  real z = (xg - gx)/alpha;
  real lpdf = log(xg + gx) - log2() - log(alpha) - log(x) + normal_lpdf(z|0, 1);
  return lpdf;
}

real fatigue_lccdf(real x, real alpha, real gamma){
  real xg = sqrt(x/gamma);
  real gx = sqrt(gamma/x);
  // real xg = exp(lmultiply(0.5, x) - lmultiply(0.5, gamma));
  // real gx = exp(lmultiply(0.5, gamma) - lmultiply(0.5, x));
  real z = ( xg - gx )/alpha;
  return normal_lcdf(-z|0, 1);
}

// *****************************************************************
// Gompertz distribution:

real gompertz_lpdf(real x, real alpha, real gamma){
  real lpdf = log(alpha) + log(gamma) + gamma*x - alpha*expm1(gamma*x);
  return lpdf;
}

real gompertz_lccdf(real x, real alpha, real gamma){
  real lsurv = - alpha*expm1(gamma*x);
  return lsurv;
}

// *****************************************************************
// Generalized Gamma distribution (Stacy's original parametrization):

real ggstacy_lpdf(real time, real alpha, real gamma, real kappa){
  real lpdf;
  lpdf = log(kappa) - lmultiply(alpha, gamma) - lgamma(alpha/kappa) +
         lmultiply(alpha-1, time) - pow(time/gamma, kappa);
  return(lpdf);
}

real ggstacy_lccdf(real time, real alpha, real gamma, real kappa){
  real lccdf;
  real x = pow(time/gamma, kappa);
  real shape = alpha/kappa;
  lccdf = gamma_lccdf(x|shape, 1);
  return(lccdf);
}


// *****************************************************************
// Generalized Gamma distribution (Prentice's parametrization):

real ggprentice_lpdf(real x, real mu, real sigma, real varphi){
  real lpdf = 0;
  real w;
  real abs_q;
  real q2;
  real qw;
  w = (log(x) - mu) / sigma;
  abs_q = abs(varphi);
  q2 = pow(varphi, -2);
  qw = varphi * w;
  if(varphi != 0){
    lpdf += - log(sigma) - log(x) + lmultiply(1 - 2 * q2, abs_q) + q2 * (qw - exp(qw)) - lgamma(q2) ;
  }else{
    lpdf += lognormal_lpdf(x|mu, sigma);
  }

  return(lpdf);
}

real ggprentice_lccdf(real x, real mu, real sigma, real varphi){
  real lccdf;
  real w;
  real q2;
  real aux;
  if(varphi == 0){
    lccdf = lognormal_lccdf(x|mu, sigma);
  }else{
    w = (log(x) - mu)/sigma;
    q2 = pow(varphi, -2);
    aux = q2*exp(varphi*w);
    if(varphi > 0){
      lccdf = gamma_lccdf(aux|q2, 1);
    }else if(varphi < 0){
      lccdf = gamma_lcdf(aux|q2, 1);
    }
  }
  return(lccdf);
}
