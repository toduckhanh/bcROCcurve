/*
 Estimation Area Under ROC Curve in Presence of Nonignorable missing data.
 Author: Duc-Khanh To
 */

#include <Rcpp.h>

using namespace Rcpp;

//-------------------------------------------------------------------------//
static inline double indauc(double a, double b) {
  if(a > b) return 0;
  else{
    double d = 0;
    if(a == b) d += 1;
    return (8 - 3*d)/(8 + 2*d);
  }
}

//-------------------------------------------------------------------------//
// Computing FULL estimators.
// [[Rcpp::export]]
double aucFull_C(NumericVector tt, NumericVector dd) {
  int n = tt.size();
  double den = 0.0, num = 0.0, temp;
  //double ans = 0.0;
  int i, j;
  for(i = 0; i < n; ++i){
    for(j = 0; j < n; ++j){
      if(j != i){
        temp = (1 - dd[i])*dd[j];
        num += indauc(tt[i],tt[j])*temp;
        den += temp;
      }
    }
  }
  return num/den;
}

//-------------------------------------------------------------------------//
// Computing parallel FI, MSI, IPW and DR estimators.
// [[Rcpp::export]]
NumericVector aucMAR_C(NumericVector tt, NumericVector dd, NumericVector vv, NumericVector rho, NumericVector pi){
  int nn = tt.size();
  double den_fi = 0.0, num_fi = 0.0, temp_fi = 0.0;
  double den_msi = 0.0, num_msi = 0.0, temp_msi = 0.0;
  double den_ipw = 0.0, num_ipw = 0.0, temp_ipw = 0.0;
  double den_dr = 0.0, num_dr = 0.0, temp_dr = 0.0;
  NumericVector Vpi = vv/pi;
  // Define D_MSI
  NumericVector ddMSI(nn);
  // Define D_DR
  NumericVector ddDR(nn);
  // Computing D_MSI and D_DR
  ddMSI = dd*vv + (1.0 - vv)*rho;
  ddDR = dd*Vpi - (Vpi - 1.0)*rho;
  double I_ij = 0.0;
  for(int i = 0; i < nn; i++){
    for(int j = 0; j < nn; j++){
      if(j != i){
        I_ij = indauc(tt[i], tt[j]);
        temp_fi = rho[j]*(1 - rho[i]);
        num_fi += I_ij*temp_fi;
        den_fi += temp_fi;
        temp_msi = ddMSI[j]*(1 - ddMSI[i]);
        num_msi += I_ij*temp_msi;
        den_msi += temp_msi;
        temp_ipw = dd[j]*Vpi[j]*(1 - dd[i])*Vpi[i];
        num_ipw += I_ij*temp_ipw;
        den_ipw += temp_ipw;
        temp_dr = ddDR[j]*(1 - ddDR[i]);
        num_dr += I_ij*temp_dr;
        den_dr += temp_dr;
      }
    }
  }
  NumericVector out = NumericVector::create(
    Named("fi") = 0.0, Named("msi") = 0.0, Named("ipw") = 0.0, Named("dr") = 0.0);
  out("fi") = num_fi/den_fi;
  out("msi") = num_msi/den_msi;
  out("ipw") = num_ipw/den_ipw;
  out("dr") = num_dr/den_dr;
  return out;
}

//-------------------------------------------------------------------------//
// Computing FI estimators.
// [[Rcpp::export]]
double aucFI_C(NumericVector tt, NumericVector rho){
  int n = tt.size();
  double den = 0.0, num = 0.0, temp;
  int i, j;
  for(i = 0; i < n; ++i){
    for(j = 0; j < n; ++j){
      if(j != i){
        temp = rho[j]*(1 - rho[i]);
        num += indauc(tt[i],tt[j])*temp;
        den += temp;
      }
    }
  }
  return num/den;
}

//-------------------------------------------------------------------------//
// Computing MSI estimators.
// [[Rcpp::export]]
double aucMSI_C(NumericVector tt, NumericVector dd, NumericVector vv, NumericVector rho){
  int n = tt.size();
  double den = 0.0, num = 0.0, temp;
  NumericVector ddMSI(n);
  int i, j;
  // Creat D_MSI.
  ddMSI = dd*vv + (1.0 - vv)*rho;
  for(i = 0; i < n; ++i){
    for(j = 0; j < n; ++j){
      if(j != i){
        temp = ddMSI[j]*(1 - ddMSI[i]);
        num += indauc(tt[i],tt[j])*temp;
        den += temp;
      }
    }
  }
  return num/den;
}

//-------------------------------------------------------------------------//
// Computing IPW estimators.
// [[Rcpp::export]]
double aucIPW_C(NumericVector tt, NumericVector dd, NumericVector vv, NumericVector pi) {
  int n = tt.size();
  double den = 0.0, num = 0.0, temp;
  int i, j;
  NumericVector Vpi = vv/pi;
  // NumericVector ddIPW(n);
  // Create D_IPW.
  // ddIPW = dd*Vpi;
  for(i = 0; i < n; ++i){
    for(j = 0; j < n; ++j){
      if(j != i){
        temp = dd[j]*Vpi[j]*(1 - dd[i])*Vpi[i];
        num += indauc(tt[i],tt[j])*temp;
        den += temp;
      }
    }
  }
  return num/den;
}

//-------------------------------------------------------------------------//
// Computing DR estimators.
// [[Rcpp::export]]
double aucDR_C(NumericVector tt, NumericVector dd, NumericVector vv, NumericVector rho,
                  NumericVector pi) {
  int n = tt.size();
  double den = 0.0, num = 0.0, temp;
  NumericVector ddDR(n);
  NumericVector Vpi = vv/pi;
  int i, j;
  // Creat D_DR.
  ddDR = dd*Vpi - (Vpi - 1.0)*rho;
  for(i = 0; i < n; ++i){
    for(j = 0; j < n; ++j){
      if(j != i){
        temp = ddDR[j]*(1 - ddDR[i]);
        num += indauc(tt[i],tt[j])*temp;
        den += temp;
      }
    }
  }
  return num/den;
}
