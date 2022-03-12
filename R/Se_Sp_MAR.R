### Var for (Sp(c), Se(c)) at fixed cut-off c
h_deriv <- function(theta, n.nuis){
  res <- matrix(0, nrow = (3 + n.nuis), ncol = 2)
  res[,1] <- c(theta[2]/theta[1]^2, 1/theta[1], 0, rep(0,n.nuis))
  res[,2] <- c(theta[3]/(1 - theta[1])^2, 0, 1/(1 - theta[1]), rep(0,n.nuis))
  return(t(res))
}

## Estimating function of Imputation Estimator (FI and MSI).
EstFuncIE <- function(dd, vv, tt, cpt, rho_score, rho.hat, theta, m){
  n <- length(vv)
  n.par <- 3 + ncol(rho_score)
  res <- matrix(0, nrow = n, ncol = n.par)
  res[,1] <- vv*(m*dd - theta[1] + (1 - m)*rho.hat) + (1 - vv)*(rho.hat - theta[1])
  res[,2] <- vv*((tt > cpt)*(m*dd + (1 - m)*rho.hat) - theta[2]) + (1 - vv)*((tt > cpt)*rho.hat - theta[2])
  res[,3] <- vv*((tt <= cpt)*(m*(1 - dd) + (1 - m)*(1 - rho.hat)) - theta[3]) +
    (1 - vv)*((tt <= cpt)*(1 - rho.hat) - theta[3])
  res[,4:n.par] <- rho_score
  return(res)
}

## Derivative of Estimating function of Imputation Estimator (FI and MSI).
EstFuncIE_deriv <- function(vv, tt, X, cpt, rho.hat, rho.hes, theta, m){
  n.par <- 3 + ncol(rho.hes)
  der_rho <- X*rho.hat*(1 - rho.hat)
  temp <- (1 - m*vv)*der_rho
  A <- colSums(temp)
  B1 <- colSums((tt > cpt)*temp)
  B2 <- colSums((tt <= cpt)*(-temp))
  B <- rbind(B1, B2)
  res.left <- rbind(diag(-length(tt), 3, 3), matrix(0, nrow = ncol(rho.hes), ncol = 3))
  res.right <- rbind(A, B, -rho.hes)
  res <- cbind(res.left, res.right)
  return(res)
}

## Estimating function of IPW Estimator
EstFuncIPW <- function(dd, vv, tt, cpt, pi_score, pi.hat, theta){
  n.par <- 3 + ncol(pi_score)
  res <- matrix(0, nrow = length(tt), ncol = n.par)
  res[,1] <- vv*(dd - theta[1])/pi.hat
  res[,2] <- vv*((tt > cpt)*dd - theta[2])/pi.hat
  res[,3] <- vv*((tt <= cpt)*(1 - dd) - theta[3])/pi.hat
  res[,4:n.par] <- pi_score
  return(res)
}

## Derivative of Estimating function of IPW Estimator
EstFuncIPW_deriv <- function(dd, vv, tt, X, cpt, pi.hat, pi.hes, theta){
  n.par <- 3 + ncol(pi.hes)
  der_pi_inv <- -X*(1 - pi.hat)/pi.hat
  term <- vv*der_pi_inv
  A <- colSums(term*(dd - theta[1]))
  B1 <- colSums(term*((tt > cpt)*dd - theta[2]))
  B2 <- colSums(term*((tt <= cpt)*(1 - dd) - theta[3]))
  res.left <- rbind(diag(-sum(vv/pi.hat), 3, 3), matrix(0, nrow = ncol(pi.hes), ncol = 3))
  res.right <- rbind(A, B1, B2, -pi.hes)
  res <- cbind(res.left, res.right)
  return(res)
}

## Estimating function of SPE Estimator
EstFuncDR <- function(dd, vv, tt, cpt, rho.hat, rho_score, pi.hat, pi_score, theta){
  n.par <- 3 + ncol(pi_score) + ncol(rho_score)
  res <- matrix(0, nrow = length(tt), ncol = n.par)
  res[,1] <- vv*(dd - theta[1])/pi.hat - (rho.hat - theta[1])*(vv - pi.hat)/pi.hat
  res[,2] <- vv*((tt > cpt)*dd - theta[2])/pi.hat - (vv - pi.hat)*((tt > cpt)*rho.hat - theta[2])/pi.hat
  res[,3] <- vv*((tt <= cpt)*(1 - dd) - theta[3])/pi.hat - (vv - pi.hat)*((tt <= cpt)*(1 - rho.hat) - theta[3])/pi.hat
  res[,4:(3 + ncol(rho_score))] <- rho_score 
  res[,(4 + ncol(rho_score)):n.par] <- pi_score
  return(res)
}

## Derivative of Estimating function of SPE Estimator
EstFuncDR_deriv <- function(dd, vv, tt, X_dise, X_veri, cpt, rho.hat, rho.hes, pi.hat, pi.hes, theta){
  n.par <- 3 + ncol(rho.hes) + ncol(pi.hes)
  der_rho <- X_dise*rho.hat*(1 - rho.hat)
  temp.v <- (1 - vv/pi.hat)
  temp <- temp.v*der_rho
  der_pi_inv <- -X_veri*(1 - pi.hat)/pi.hat
  temp2 <- vv*der_pi_inv
  H <- colSums(temp)
  G1 <- colSums((tt > cpt)*temp)
  G2 <- colSums((tt <= cpt)*(-temp))
  G <- rbind(G1, G2)
  D <- colSums(temp2*(dd - rho.hat))
  E1 <- colSums(temp2*((tt > cpt)*(dd - rho.hat)))
  E2 <- colSums(temp2*((tt <= cpt)*(rho.hat - dd)))
  res.left <- rbind(diag(-length(tt), 3, 3), matrix(0, nrow = ncol(rho.hes) + ncol(pi.hes), ncol = 3))
  res.right.rho <- rbind(H, G, -rho.hes, diag(0, nrow = ncol(pi.hes), ncol = ncol(rho.hes)))
  res.right.pi <- rbind(D, E1, E2, diag(0, nrow = ncol(rho.hes), ncol = ncol(pi.hes)), -pi.hes)
  res <- cbind(res.left, res.right.rho, res.right.pi)
  return(res)
}

#' @export
Sen_Spe_FI <- function(out_dise_model, tt, dd, vv, cpt){
  rho_est <- out_dise_model$rho_est
  theta <- c(mean(rho_est), mean((tt > cpt)*rho_est), mean((tt <= cpt)*(1 - rho_est)))
  rho_hess <- out_dise_model$rho_hess
  rho_score <- out_dise_model$rho_score
  mat.dise <- out_dise_model$mat.dise
  npar.dise <- out_dise_model$npar.dise
  term1 <- EstFuncIE_deriv(vv, tt, mat.dise, cpt, rho_est, rho_hess, theta, m = 0)
  term2 <- EstFuncIE(dd, vv, tt, cpt = cpt, rho_score = rho_score, rho.hat = rho_est, theta, m = 0)
  term1_inv <- solve(term1)
  Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
  hh <- h_deriv(theta, npar.dise)
  var_se_sp <- hh %*% Sig %*% t(hh)
  std_se_sp <- sqrt(diag(var_se_sp))
  se_sp <- c(theta[2]/theta[1], theta[3]/(1 - theta[1]))
  return(list(se_sp = se_sp, std_se_sp = std_se_sp, var_se_sp = var_se_sp))
}

#' @export
Sen_Spe_MSI <- function(out_dise_model, tt, dd, vv, cpt){
  rho_est <- out_dise_model$rho_est
  D_MSI <- data_analy$V*data_analy$D + (1 - data_analy$V)*rho_est
  theta <- c(mean(D_MSI), mean((tt > cpt)*D_MSI), mean((tt <= cpt)*(1 - D_MSI)))
  rho_hess <- out_dise_model$rho_hess
  rho_score <- out_dise_model$rho_score
  mat.dise <- out_dise_model$mat.dise
  npar.dise <- out_dise_model$npar.dise
  term1 <- EstFuncIE_deriv(vv, tt, mat.dise, cpt, rho_est, rho_hess, theta, m = 1)
  term2 <- EstFuncIE(dd, vv, tt, cpt = cpt, rho_score = rho_score, rho.hat = rho_est, theta, m = 1)
  term1_inv <- solve(term1)
  Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
  hh <- h_deriv(theta, npar.dise)
  var_se_sp <- hh %*% Sig %*% t(hh)
  std_se_sp <- sqrt(diag(var_se_sp))
  se_sp <- c(theta[2]/theta[1], theta[3]/(1 - theta[1]))
  return(list(se_sp = se_sp, std_se_sp = std_se_sp, var_se_sp = var_se_sp))
}

#' @export
Sen_Spe_IPW <- function(out_veri_model, tt, dd, vv, cpt){
  pi_est <- out_veri_model$pi_est
  pi_score <- out_veri_model$pi_score
  pi_hess <- out_veri_model$pi_hess
  mat.veri <- out_veri_model$mat.veri
  npar.veri <- out_veri_model$npar.veri
  D_IPW <- vv*dd/pi_est
  D0_IPW <- vv*(1 - dd)/pi_est
  theta <- c(mean(D_IPW), mean((tt > cpt)*D_IPW), mean((tt <= cpt)*D0_IPW))
  term1 <- EstFuncIPW_deriv(dd, vv, tt, mat.veri, cpt, pi_est, pi_hess, theta)
  term2 <- EstFuncIPW(dd, vv, tt, cpt = cpt, pi_score = pi_score, pi.hat = pi_est, theta)
  term1_inv <- solve(term1)
  Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
  hh <- h_deriv(theta, npar.veri)
  var_se_sp <- hh %*% Sig %*% t(hh)
  std_se_sp <- sqrt(diag(var_se_sp))
  se_sp <- c(theta[2]/theta[1], theta[3]/(1 - theta[1]))
  return(list(se_sp = se_sp, std_se_sp = std_se_sp, var_se_sp = var_se_sp))
}

#' @export
Sen_Spe_DR <- function(out_dise_model, out_veri_model, tt, dd, vv, cpt){
  rho_est <- out_dise_model$rho_est
  rho_hess <- out_dise_model$rho_hess
  rho_score <- out_dise_model$rho_score
  mat.dise <- out_dise_model$mat.dise
  npar.dise <- out_dise_model$npar.dise
  pi_est <- out_veri_model$pi_est
  pi_score <- out_veri_model$pi_score
  pi_hess <- out_veri_model$pi_hess
  mat.veri <- out_veri_model$mat.veri
  npar.veri <- out_veri_model$npar.veri
  D_DR <- vv*dd/pi_est - (vv - pi_est)*rho_est/pi_est
  theta <- c(mean(D_DR), mean((tt > cpt)*D_DR), mean((tt <= cpt)*(1 - D_DR)))
  term1 <- EstFuncDR_deriv(dd, vv, tt, mat.dise, mat.veri, cpt, rho_est, rho_hess, pi_est, pi_hess, theta)
  term2 <- EstFuncDR(dd, vv, tt, cpt = cpt, rho_est, rho_score, pi_est, pi_score, theta)
  term1_inv <- solve(term1)
  Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
  hh <- h_deriv(theta, npar.veri + npar.dise)
  var_se_sp <- hh %*% Sig %*% t(hh)
  std_se_sp <- sqrt(diag(var_se_sp))
  se_sp <- c(theta[2]/theta[1], theta[3]/(1 - theta[1]))
  return(list(se_sp = se_sp, std_se_sp = std_se_sp, var_se_sp = var_se_sp))
}

