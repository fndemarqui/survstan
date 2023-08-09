
// *****************************************************************
// Loglogistic distribution:


real loglogistic_lpdf(real x, real alpha, real gamma){
  real aux = log1p(pow(x/gamma, alpha));
  real lpdf = log(alpha) - log(gamma) + lmultiply(alpha-1, x) -
  lmultiply(alpha-1, gamma) - 2*aux;
  return lpdf;
}

real loglogistic_lccdf(real x, real alpha, real gamma){
  real lsurv = - log1p(pow(x/gamma, alpha));
  return lsurv;
}


// real loglogistic_lpdf(real y, real mu, real sigma){
//     return logistic_lpdf(log(y)|mu, sigma) - log(y);
// }
//
// real loglogistic_lcdf(real y, real mu, real sigma){
//     return logistic_lcdf(log(y)|mu, sigma);
// }
//
// real loglogistic_lccdf(real y, real mu, real sigma){
//     return logistic_lccdf(log(y)|mu, sigma);
// }
//
// real loglogistic_rng(real mu, real sigma){
//     return exp(logistic_rng(mu, sigma));
// }
//


// *****************************************************************
// Birnbaum–Saunders (fatigue) distribution

  real fatigue_lpdf(real x, real alpha, real gamma){
    real aux1 = sqrt(x/gamma) + sqrt(gamma/x);
    real aux2 = ( sqrt(x/gamma) - sqrt(gamma/x) )/alpha;
    real lpdf = log(aux1) - log(2*alpha*x) + normal_lpdf(aux2|0, 1);
    return lpdf;
  }

  real fatigue_lccdf(real x, real alpha, real gamma){
    real aux = ( sqrt(x/gamma) - sqrt(gamma/x) )/alpha;
    return normal_lcdf(-aux|0, 1);
  }


