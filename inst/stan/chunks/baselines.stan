

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
  real lpdf = log(alpha) + gamma*x - alpha*expm1(gamma*x)/gamma;
  return lpdf;
}

real gompertz_lccdf(real x, real alpha, real gamma){
  real lsurv = - alpha*expm1(gamma*x)/gamma;
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
  if(varphi != 0){
    w = (log(x) - mu) / sigma;
    abs_q = abs(varphi);
    q2 = pow(varphi, -2);
    qw = varphi * w;
    lpdf += - log(sigma*x) + lmultiply(1 - 2 * q2, abs_q) + q2 * (qw - exp(qw)) - lgamma(q2) ;
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
  if(varphi != 0){
    w = (log(x) - mu)/sigma;
    q2 = pow(varphi, -2);
    aux = q2*exp(varphi*w);
    if(varphi > 0){
    lccdf = gamma_lccdf(aux|q2, 1);
  }else if(varphi < 0){
      lccdf = gamma_lcdf(aux|q2, 1);
    }
  }else{
    lccdf = lognormal_lccdf(x|mu, sigma);
  }
  return(lccdf);
}



// *****************************************************************
// Bernstein Polynomials:

real beta_pdf(real x, real a, real b){
  real pdf = 0;
  if(x == 0 && a == 1 || x == 1 && b == 1){
    pdf = exp(-lbeta(a, b));
  }else{
    pdf = exp(beta_lpdf(x| a, b));
  }
  return pdf;
}

// For compatibility with R dbeta function:
// real beta_pdf(real x, real a, real b){
//   real pdf = 0;
//   if(x == 0 && a == 1 || x == 1 && b == 1){
//     pdf = exp(-lbeta(a, b));
//   }else{
//     pdf = exp(beta_lpdf(x| a, b));
//   }
//   return pdf;
// }

// real bernstein_lpdf(real y, vector xi){
//   int m = num_elements(xi);
//   real lpdf;
//   real ht = 0;
//   real Ht = 0;
//   for(j in 1:m){
//     ht += beta_pdf(y, j, m - j + 1)*xi[j];
//     Ht += beta_cdf(y| j, m - j + 1)*xi[j];
//   }
//   lpdf = log(ht) - Ht;
//   return(lpdf);
// }
//
//
// real bernstein_lccdf(real y, vector xi){
//   int m = num_elements(xi);
//   real Ht = 0;
//   for(j in 1:m){
//     Ht += beta_cdf(y| j, m - j + 1)*xi[j];
//   }
//   return(-Ht);
// }

real bernstein_lpdf(real y, vector xi){
  int m = num_elements(xi);
  real lpdf;
  row_vector[m] g;
  row_vector[m] G;
  for(j in 1:m){
    g[j] = beta_pdf(y, j, m - j + 1);
    G[j] = exp( beta_lcdf(y| j, (m - j + 1)) );
  }
  lpdf = log(g*xi) - G*xi;
  return(lpdf);
}

real bernstein_lccdf(real y, vector xi){
  int m = num_elements(xi);
  row_vector[m] G;
  for(j in 1:m){
    G[j] = exp( beta_lcdf(y| j, (m - j + 1)) );
  }
  return(-G*xi);
}

vector bernstein_vlpdf(matrix g, matrix G, vector xi){
  return log(g*xi) - G*xi;
}

vector bernstein_vlccdf(matrix G, vector xi){
  return -G*xi;
}


// *****************************************************************
// Piecewise exponential distribution:

 matrix TTT(vector time, vector rho){
  int n = num_elements(time);
  int m = num_elements(rho) - 1;
  vector[2] aux;
  matrix[n, m] ttt;
  for(i in 1:n){
      for(j in 1:m){
        aux[1] = time[i];
        aux[2] = rho[j+1];
        ttt[i, j] = (min(aux) - rho[j])*(time[i] - rho[j] > 0);
      }
  }
  return(ttt);
}


int[] IDT(vector time, vector rho){
  int n = num_elements(time);
  int j;
  int idt[n];
    for(i in 1:n){
      j = 0;
      while(time[i] > rho[j+1]){
        j = j + 1;
      }
      idt[i] = j;
  }
  return(idt);
}

// version 1:
vector piecewise_vlpdf(matrix ttt, int[] idt, vector xi){
  return log(xi[idt]) - ttt*xi;
}

vector piecewise_vlccdf(matrix ttt, vector xi){
  return -ttt*xi;
}

