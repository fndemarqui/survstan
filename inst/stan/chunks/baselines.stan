

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

