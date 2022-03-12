/*
 * This code computes the asymptotic variances of these VUS estimates.
 * We use the RcppArmadillo package to compute the transpose and inverse of 
 * matrix or vector. 
 * 
 * Note that we use logistic regressions for the 
 * verification and disease model, respectively.
 * The verification model is:
 * 					logit(lambda D + beta_0 + beta_1 T + beta_2 A)
 * The disease model is:
 * 					logit(alpha_0 + alpha_1 T  + alpha_2 A)
 * The number of parameters is 11 (default)
 */

//[[Rcpp::depends(RcppArmadillo)]]

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------//
static inline double indauc(double a, double b) {
  if(a > b) return 0;
  else{
    double d = 0;
    if(a == b) d += 1;
    return (8 - 3*d)/(8 + 2*d);
  }
}

//------------------------------------------------------------------------------//
/*
 * estimating function of AUC
 * Note:
 * dd_i, dd_j, dd_3k are estimated verison of D_1i, D_2j and D_3k, respectively.
 * E.g., FI: dd_1i, dd_2j and dd_3k correspond to rho_1i, rho_2j and rho_3k.
 * bias_ijk = I_ijk - mu_hat
 */

static inline double aucEstFunc(double dd1_i, double dd2_j, double bias_ij){
  return dd1_i*dd2_j*bias_ij;
}
//==============================================================================//
/*
 * derivative of estimating function of AUC
 * Note:
 * deri_dd_i and deri_dd_j are the derivatives of estimated version
 * of D_i, D_j.
 */

static inline vec aucDerEstFunc(double dd1_i, double dd2_j, vec deri_dd1_i, vec deri_dd2_j, double bias_ij){
  vec tem1 = deri_dd1_i*dd2_j;
  tem1 += dd1_i*deri_dd2_j;
  return tem1*bias_ij;
}

//==============================================================================//
// mu_hat == AUC estimated
//[[Rcpp::export]]
NumericVector asyVarAUC_C(NumericVector tt, NumericMatrix D_hat, double mu_hat, NumericMatrix Score, 
                          NumericMatrix Hess_inv, NumericMatrix Der_D1_hat, NumericMatrix Der_D2_hat){
  vec T = as<vec>(tt);
  mat Dhat = as<mat>(D_hat);
  mat SS = as<mat>(Hess_inv);
  mat U2 = as<mat>(Score);
  mat der_Dhat1 = as<mat>(Der_D1_hat);
  mat der_Dhat2 = as<mat>(Der_D2_hat);
  uword nn = T.n_elem;
  vec S1(nn, fill::zeros), U1(SS.n_rows, fill::zeros);
  for(uword i = 0; i < nn; ++i){
    double temp2 = 0.0;
    double t_i = T[i];
    rowvec Dhat_i = Dhat.row(i);
    for(uword j = 0; j < nn; ++j){
      if(j != i){
        double t_j = T[j];
        rowvec Dhat_j = Dhat.row(j);
        double i_ij = indauc(t_i, t_j) - mu_hat;
        double i_ji = indauc(t_j, t_i) - mu_hat;
        double temp1 = aucEstFunc(Dhat_i[0], Dhat_j[1], i_ij);
        temp1 += aucEstFunc(Dhat_j[0], Dhat_i[1], i_ji);
        temp2 += temp1;
        U1 += aucDerEstFunc(Dhat_i[0], Dhat_j[1], der_Dhat1.col(i), der_Dhat2.col(j), i_ij);
      }
    }
    S1[i] = temp2;
  }
  vec Q_auc = (S1 + U2*SS*U1)/nn;
  return wrap(Q_auc);
}

