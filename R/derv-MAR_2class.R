derv_pi_inv <- function(X, pi_hat, npar.pi, npar.rho){
  res1 <- -X*(1 - pi_hat)/pi_hat
  res2 <- matrix(0, nrow = nrow(X), ncol = npar.rho)
  res <- cbind(res1, res2)
  return(res)
}

derv_rho <- function(X, rho_hat, npar.pi){
  res_null <- matrix(0, nrow = nrow(X), ncol = npar.pi)
  res <- X*rho_hat*(1 - rho_hat)
  return(cbind(res_null, res))
}

## Estimating function for pi = Pr(V = 1| T, A).
piEstFunc_MAR <- function(vv, u_matrix, pi_hat, pi_cov,
                          model = c("logit", "probit", "threshold")){
  model <- match.arg(model)
  res <- switch(EXPR = model,
                logit = (vv - pi_hat)*u_matrix,
                probit = u_matrix*as.numeric(dnorm(u_matrix %*% pi_cov))*
                  (vv - pi_hat)/(pi_hat*(1 - pi_hat)),
                threshold = (vv/pi_hat - 1)*u_matrix/(1 - pi_hat))
  return(res)
}

## the first derivatives of 1/pi.
pi_inv_deriv_MAR <- function(u_matrix, pi_hat, pi_cov,
                             model = c("logit", "probit", "threshold")){
  # u_matrix: the design matrix is used to fit the verification model.
  model <- match.arg(model)
  res <- switch(EXPR = model,
                logit = -u_matrix*(1 - pi_hat)/pi_hat,
                probit = -u_matrix*as.numeric(dnorm(u_matrix %*% pi_cov))/pi_hat,
                threshold = -u_matrix/pi_hat^2)
  return(res)
}

## Estimating function for rho = Pr(D = 1 | T, A).
rhoEstFunc_MAR <- function(dd, vv, u_matrix, rho_hat){
  return(u_matrix*vv*(dd - rho_hat))
}

