

vector loglik_aft(vector lpdf, vector lsurv, vector event, vector lp, real tau){
  int n = num_elements(lpdf);
  vector[n] loglik;
  loglik = event  .* (lpdf - lp - log(tau)) + (1-event) .* lsurv;
  return loglik;
}


vector loglik_ah(vector lpdf, vector lsurv, vector event, vector lp, real tau){
  int n = num_elements(lpdf);
  vector[n] loglik;
  // loglik = event  .* (lpdf - lsurv - log(tau)) +  exp(-lp) .* lsurv ;
  loglik = event  .* (lpdf - lsurv - log(tau)) +  exp(lp) .* lsurv ;
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


vector loglik_yp(vector status, vector lpdf, vector lsurv, vector lp_short, vector lp_long, vector K, real tau){

  int n = num_elements(lpdf);
  vector[n] Rt0;
  vector[n] log_ht;
  vector[n] log_St;
  vector[n] loglik;
  vector[n] theta;
  vector[n] aux;
  vector[n] lht0 = lpdf - lsurv - log(tau);
  vector[n] Ht0 = -lsurv;

  Rt0 = expm1(Ht0);
  theta = exp(lp_long);

  aux = K .* Rt0;
  log_ht = lp_short - log1p(aux) + lht0 + Ht0;
  log_St = -theta .* log1p(aux);
  loglik = status .* log_ht + log_St;

  return loglik;
}


vector loglik_eh(vector status, vector lpdf, vector lsurv, vector lp, vector K, real tau){

  int n = num_elements(lpdf);
  vector[n] loglik;
  loglik = status .* (lp + lpdf - lsurv - log(tau)) + lsurv .* K;

  return loglik;
}
