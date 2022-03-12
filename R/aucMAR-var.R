## Estimate for asymptotic variances for AUC
#' @useDynLib bcROCcurve
#' @importFrom Rcpp evalCpp

#' @export
asyVarAUC_FI <- function(T, Matrix.dise, rho.hat, score, hess, npar.pi, auc.est){
  n <- length(T)
  der_rho <- t(derv_rho(Matrix.dise, rho.hat, npar.pi))
  der_1_rho <- (-1)*der_rho
  rho.hat_mat <- cbind(1 - rho.hat, rho.hat)
  Q_Fi <- asyVarAUC_C(T, rho.hat_mat, auc.est, score, solve(hess), der_1_rho, der_rho)
  return(n^2*sum(Q_Fi^2)/prod(colSums(rho.hat_mat)^2))
}

#' @export
asyVarAUC_MSI <- function(T, D, V, Matrix.dise, rho.hat, score, hess, npar.pi, auc.est){
  n <- length(T)
  D_MSI <- V*D + (1 - V)*rho.hat
  D_MSI_mat <- cbind(1 - D_MSI, D_MSI)
  der_D_MSI <- t((1 - V)*derv_rho(Matrix.dise, rho.hat, npar.pi))
  der_1_D_MSI <- (-1)*der_D_MSI
  Q_Msi <- asyVarAUC_C(T, D_MSI_mat, auc.est, score, solve(hess), der_1_D_MSI, der_D_MSI)
  return(n^2*sum(Q_Msi^2)/prod(colSums(D_MSI_mat)^2))
}

#' @export
asyVarAUC_IPW <- function(T, D, V, Matrix.veri, pi.hat, score, hess, npar.pi, npar.rho, auc.est){
  n <- length(T)
  pi_inv.der <- derv_pi_inv(Matrix.veri, pi.hat, npar.pi, npar.rho)
  D_IPW <- cbind(1 - D, D)*V/pi.hat
  der_D_IPW1 <- t(V*(1 - D)*pi_inv.der)
  der_D_IPW2 <- t(V*D*pi_inv.der)
  Q_Ipw <- asyVarAUC_C(T, D_IPW, auc.est, score, solve(hess), der_D_IPW1, der_D_IPW2)
  return(sum(Q_Ipw^2)/prod((colSums(D_IPW)/sum(V/pi.hat))^2)/n^2)
}

#' @export
asyVarAUC_DR <- function(T, D, V, Matrix.veri, Matrix.dise, rho.hat, pi.hat, score, hess, npar.pi, npar.rho, auc.est){
  n <- length(T)
  D_DR <- V*D/pi.hat - (V/pi.hat - 1)*rho.hat
  D_DR_mat <- cbind(1 - D_DR, D_DR)
  der_rho <- derv_rho(Matrix.dise, rho.hat, npar.pi)
  pi_inv.der <- derv_pi_inv(Matrix.veri, pi.hat, npar.pi, npar.rho)
  der_D_DR <- t(V*(D - rho.hat)*pi_inv.der - (V/pi.hat - 1)*der_rho)
  der_1_D_DR <- (-1)*der_D_DR
  Q_Dr <- asyVarAUC_C(T, D_DR_mat, auc.est, score, solve(hess), der_1_D_DR, der_D_DR)
  return(n^2*sum(Q_Dr^2)/prod(colSums(D_DR_mat)^2))
}

