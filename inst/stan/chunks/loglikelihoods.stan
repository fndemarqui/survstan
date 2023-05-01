

vector loglik_aft(vector lpdf, vector lsurv, vector event, vector lp, real tau){
  int n = num_elements(lpdf);
  vector[n] loglik;
  loglik = event  .* (lpdf - lp - log(tau)) + (1-event) .* lsurv;
  return loglik;
}


vector loglik_ah(vector lpdf, vector lsurv, vector event, vector lp, real tau){
  int n = num_elements(lpdf);
  vector[n] loglik;
  loglik = event  .* (lpdf - lsurv - log(tau)) +  exp(-lp) .* lsurv ;
  return loglik;
}


vector loglik_ph(vector lpdf, vector lsurv, vector event, vector lp, real tau){
  int n = num_elements(lpdf);
  vector[n] loglik;
  vector[n] lht = lpdf - lsurv;
  loglik = event  .* (lht + lp - log(tau)) +  exp(lp) .* lsurv ;
  return loglik;
}


vector loglik_po(vector lpdf, vector lsurv, vector event, vector lp, real tau){
  int n = num_elements(lpdf);
  vector[n] lht = lpdf - lsurv;
  vector[n] Ht = - lsurv;
  vector[n] aux = exp(lp) .* expm1(Ht);
  vector[n] loglik = event .* (lht + lp + Ht - log(tau)) - (1+event) .* log1p(aux);
  return loglik;
}
