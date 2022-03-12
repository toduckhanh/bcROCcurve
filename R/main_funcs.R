
#' @export
auc_MAR <- function(out_dise_model, out_veri_model, tt, dd, vv){
  rho_est <- out_dise_model$rho_est
  pi_est <- out_veri_model$pi_est
  para.score <- cbind(out_veri_model$pi_score, out_dise_model$rho_score)
  mat.dise <- out_dise_model$mat.dise
  mat.veri <- out_veri_model$mat.veri
  npar.veri <- ncol(out_veri_model$pi_hess)
  npar.dise <- ncol(out_dise_model$rho_hess)
  para.hess <- rbind(cbind(out_veri_model$pi_hess, matrix(0, ncol = npar.dise, nrow = npar.veri)), 
                     cbind(matrix(0, ncol = npar.veri, nrow = npar.dise), out_dise_model$rho_hess))
  out_auc <- aucMAR(T = tt, D = dd, V = vv, rhoEst = rho_est, piEst = pi_est)
  se_auc_fi <- sqrt(asyVarAUC_FI(T = tt, Matrix.dise = mat.dise, rho.hat = rho_est, score = para.score, 
                                 hess = para.hess, npar.pi = npar.veri, auc.est = out_auc[1]))
  se_auc_msi <- sqrt(asyVarAUC_MSI(T = tt, D = dd, V = vv, Matrix.dise = mat.dise, rho.hat = rho_est,
                                   score = para.score, hess = para.hess, npar.pi = npar.veri, auc.est = out_auc[2]))
  se_auc_ipw <- sqrt(asyVarAUC_IPW(T = tt, D = dd, V = vv, Matrix.veri = mat.veri, pi.hat = pi_est, score = para.score, 
                                   hess = para.hess, npar.pi = npar.veri, npar.rho = npar.dise, auc.est = out_auc[3]))
  se_auc_dr <- sqrt(asyVarAUC_DR(T = tt, D = dd, V = vv, Matrix.veri = mat.veri, Matrix.dise = mat.dise, 
                                 rho.hat = rho_est, pi.hat = pi_est, score = para.score, hess = para.hess, 
                                 npar.pi = npar.veri, npar.rho = npar.dise, auc.est = out_auc[4]))
  ci_auc_fi <- out_auc[1] + c(-1,1)*qnorm((1 + 0.95)/2)*se_auc_fi
  ci_auc_msi <- out_auc[2] + c(-1,1)*qnorm((1 + 0.95)/2)*se_auc_msi
  ci_auc_ipw <- out_auc[3] + c(-1,1)*qnorm((1 + 0.95)/2)*se_auc_ipw
  ci_auc_dr <- out_auc[4] + c(-1,1)*qnorm((1 + 0.95)/2)*se_auc_dr
  res_auc <- cbind(out_auc, c(se_auc_fi, se_auc_msi, se_auc_ipw, se_auc_dr),
                   c(ci_auc_fi[1], ci_auc_msi[1], ci_auc_ipw[1], ci_auc_dr[1]), 
                   c(ci_auc_fi[2], ci_auc_msi[2], ci_auc_ipw[2], ci_auc_dr[2]))
  colnames(res_auc) <- c("Est.", "Std. error", "Low.ci", "Up.ci")
  rownames(res_auc) <- c("FI", "MSI", "IPW", "DR")
  printCoefmat(res_auc, has.Pvalue = FALSE, digits = 4, na.print = "--")
  return(res_auc)
}

#' @export
ROC_FI <- function(out_dise_model, tt, plot = FALSE){
  rho_est <- out_dise_model$rho_est
  cpt_full <- c(-Inf, sort(tt), Inf)
  pr.dise <- mean(rho_est)
  Spe_FI <- sapply(cpt_full, function(x) mean((tt <= x)*(1 - rho_est))/(1 - pr.dise))
  Sen_FI <- sapply(cpt_full, function(x) mean((tt > x)*rho_est)/pr.dise)
  if(plot){
    plot(1 - Spe_FI, Sen_FI, type = "l", xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
         main = "FI estimate")
  }
  return(cbind(Sen = Sen_FI, Spe = Spe_FI))
}

#' @export
ROC_MSI <- function(out_dise_model, tt, dd, vv, plot = FALSE){
  rho_est <- out_dise_model$rho_est
  D_MSI <- vv*dd + (1 - vv)*rho_est
  pr.dise <- mean(D_MSI)
  cpt_full <- c(-Inf, sort(tt), Inf)
  Spe <- sapply(cpt_full, function(x) mean((tt <= x)*(1 - D_MSI))/(1 - pr.dise))
  Sen <- sapply(cpt_full, function(x) mean((tt > x)*D_MSI)/pr.dise)
  if(plot){
    plot(1 - Spe, Sen, type = "l", xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
         main = "MSI estimate")
  }
  return(cbind(Sen = Sen, Spe = Spe))
}

#' @export
ROC_IPW <- function(out_veri_model, tt, dd, vv, plot = FALSE){
  pi_est <- out_veri_model$pi_est
  D_IPW <- vv*dd/pi_est
  D0_IPW <- vv*(1 - dd)/pi_est
  pr.dise_IPW <- mean(D_IPW)
  pr.dise0_IPW <- mean(D0_IPW)
  cpt_full <- c(-Inf, sort(tt), Inf)
  Spe <- sapply(cpt_full, function(x) mean((tt <= x)*D0_IPW)/pr.dise0_IPW)
  Sen <- sapply(cpt_full, function(x) mean((tt > x)*D_IPW)/pr.dise_IPW)
  if(plot){
    plot(1 - Spe, Sen, type = "l", xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
         main = "IPW estimate")
  }
  return(cbind(Sen = Sen, Spe = Spe))
}

#' @export
ROC_DR <- function(out_dise_model, out_veri_model, tt, dd, vv, plot = FALSE){
  rho_est <- out_dise_model$rho_est
  pi_est <- out_veri_model$pi_est
  D_DR <- vv*dd/pi_est - (vv - pi_est)*rho_est/pi_est
  pr.dise_DR <- mean(D_DR)
  cpt_full <- c(-Inf, sort(tt), Inf)
  Spe <- sapply(cpt_full, function(x) mean((tt <= x)*(1 - D_DR))/(1 - pr.dise_DR))
  Sen <- sapply(cpt_full, function(x) mean((tt > x)*D_DR)/pr.dise_DR)
  if(plot){
    plot(1 - Spe, Sen, type = "l", xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
         main = "DR estimate")
  }
  return(cbind(Sen = Sen, Spe = Spe))
}

#' @export
ROC_CC <- function(tt, dd, vv, plot = FALSE){
  cpt <- c(-Inf, sort(tt[vv == 1]), Inf)
  Sen <- sapply(cpt, function(x) mean(tt[vv == 1 & dd == 1] > x))
  Spe <- sapply(cpt, function(x) mean(tt[vv == 1 & dd == 0] <= x))
  if(plot){
    plot(1 - Spe, Sen, type = "l", xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
         main = "Complete Case estimate")
  }
  return(cbind(Sen = Sen, Spe = Spe))
}

#' @export
Sen_Spe_YI_FI <- function(out_dise_model, tt, prev){
  rho_est <- out_dise_model$rho_est
  pr.dise <- mean(rho_est)
  # cpt_full <- c(-Inf, sort(tt), Inf)
  # Spe_FI <- sapply(cpt_full, function(x) mean((tt <= x)*(1 - rho_est))/(1 - pr.dise))
  # Sen_FI <- sapply(cpt_full, function(x) mean((tt > x)*rho_est)/pr.dise)
  # cpt_YI <- cpt_full[which.max(Spe_FI + Sen_FI)]
  ff <- function(x, tt, rho_est, pr){
    Spe <- mean((tt <= x)*(1 - rho_est))/(1 - pr)
    Sen <- mean((tt > x)*rho_est)/pr
    return(Spe + Sen)
  }
  cpt_YI <- optimize(ff, interval = range(tt), tt = tt, rho_est = rho_est, pr = pr.dise, maximum = TRUE)$maximum
  Spe_YI <- mean((tt <= cpt_YI)*(1 - rho_est))/(1 - pr.dise)
  Sen_YI <- mean((tt > cpt_YI)*rho_est)/pr.dise
  plr <- Sen_YI/(1 - Spe_YI)
  nlr <- (1 - Sen_YI)/Spe_YI
  ppv <- Sen_YI*prev/(Sen_YI*prev + (1 - Spe_YI)*(1 - prev))
  npv <- Spe_YI*(1 - prev)/((1 - Sen_YI)*prev + Spe_YI*(1 - prev))
  return(c(cpt_YI = cpt_YI, Sen_YI = Sen_YI, Spe_YI = Spe_YI, plr = plr, nlr = nlr, ppv = ppv, npv = npv))
}

#' @export
Sen_Spe_YI_FI_boot <- function(data, name.test, name.dise, name.veri, B, formula.dise, prev){
  data_split <- split(data, f = as.factor(data[, name.dise]))
  data_1 <- data_split[[1]]
  data_2 <- data_split[[2]]
  data_3 <- data_split[[3]]
  n1 <- nrow(data_1)
  n2 <- nrow(data_2)
  n3 <- nrow(data_3)
  bts_fun <- function(data_1, data_2, data_3, n1, n2, n3, formula.dise, name.test, name.veri, prev){
    id1 <- sample(1:n1, n1, replace = TRUE)
    id2 <- sample(1:n2, n2, replace = TRUE)
    id3 <- sample(1:n3, n3, replace = TRUE)
    data_1_b <- data_1[id1,]
    data_2_b <- data_2[id2,]
    data_3_b <- data_3[id3,]
    data_b <- rbind(data_1_b, data_2_b, data_3_b)
    out_dise_model_b <- dise_model(formula.dise, data = data_b, name.veri = name.veri, print = FALSE)
    out <- Sen_Spe_YI_FI(out_dise_model_b, data_b[, name.test], prev = prev)
    return(out)
  }
  out_bts <- sapply(1:B, function(i){
    bts_fun(data_1 = data_1, data_2 = data_2, data_3 = data_3, n1 = n1, n2 = n2, n3 = n3,
            formula.dise = formula.dise, name.test = name.test, name.veri = name.veri, prev = prev)
  })
  return(out_bts)
}

#' @export
Sen_Spe_YI_MSI <- function(out_dise_model, tt, dd, vv, prev){
  rho_est <- out_dise_model$rho_est
  D_MSI <- vv*dd + (1 - vv)*rho_est
  pr.dise <- mean(D_MSI)
  ff <- function(x, tt, D, pr){
    Spe <- mean((tt <= x)*(1 - D))/(1 - pr)
    Sen <- mean((tt > x)*D)/pr
    return(Spe + Sen)
  }
  cpt_YI <- optimize(ff, interval = range(tt), tt = tt, D = D_MSI, pr = pr.dise, maximum = TRUE)$maximum
  Spe_YI <- mean((tt <= cpt_YI)*(1 - D_MSI))/(1 - pr.dise)
  Sen_YI <- mean((tt > cpt_YI)*D_MSI)/pr.dise
  plr <- Sen_YI/(1 - Spe_YI)
  nlr <- (1 - Sen_YI)/Spe_YI
  ppv <- Sen_YI*prev/(Sen_YI*prev + (1 - Spe_YI)*(1 - prev))
  npv <- Spe_YI*(1 - prev)/((1 - Sen_YI)*prev + Spe_YI*(1 - prev))
  return(c(cpt_YI = cpt_YI, Sen_YI = Sen_YI, Spe_YI = Spe_YI, plr = plr, nlr = nlr, ppv = ppv, npv = npv))
}

#' @export
Sen_Spe_YI_MSI_boot <- function(data, name.test, name.dise, name.veri, B, formula.dise, prev){
  data_split <- split(data, f = as.factor(data[, name.dise]))
  data_1 <- data_split[[1]]
  data_2 <- data_split[[2]]
  data_3 <- data_split[[3]]
  n1 <- nrow(data_1)
  n2 <- nrow(data_2)
  n3 <- nrow(data_3)
  bts_fun <- function(data_1, data_2, data_3, n1, n2, n3, formula.dise, name.test, name.veri, prev){
    id1 <- sample(1:n1, n1, replace = TRUE)
    id2 <- sample(1:n2, n2, replace = TRUE)
    id3 <- sample(1:n3, n3, replace = TRUE)
    data_1_b <- data_1[id1,]
    data_2_b <- data_2[id2,]
    data_3_b <- data_3[id3,]
    data_b <- rbind(data_1_b, data_2_b, data_3_b)
    out_dise_model_b <- dise_model(formula.dise, data = data_b, name.veri = name.veri, print = FALSE)
    out <- Sen_Spe_YI_MSI(out_dise_model_b, tt = data_b[, name.test], dd = data_b[, name.dise], 
                          vv = data_b[, name.veri], prev = prev)
    return(out)
  }
  out_bts <- sapply(1:B, function(i){
    bts_fun(data_1 = data_1, data_2 = data_2, data_3 = data_3, n1 = n1, n2 = n2, n3 = n3,
            formula.dise = formula.dise, name.test = name.test, name.veri = name.veri, prev = prev)
  })
  return(out_bts)
}

#' @export
Sen_Spe_YI_IPW <- function(out_veri_model, tt, dd, vv, prev){
  pi_est <- out_veri_model$pi_est
  D_IPW <- vv*dd/pi_est
  D0_IPW <- vv*(1 - dd)/pi_est
  pr.dise_IPW <- mean(D_IPW)
  pr.dise0_IPW <- mean(D0_IPW)
  # cpt_full <- c(-Inf, sort(tt), Inf)
  # Spe <- sapply(cpt_full, function(x) mean((tt <= x)*D0_IPW)/pr.dise0_IPW)
  # Sen <- sapply(cpt_full, function(x) mean((tt > x)*D_IPW)/pr.dise_IPW)
  # YI <- Spe + Sen
  # cpt_YI <- mean(cpt_full[which(YI == max(YI))])
  ff <- function(x, tt, D, D0, pr, pr0){
    Spe <- mean((tt <= x)*D0)/pr0
    Sen <- mean((tt > x)*D)/pr
    return(Spe + Sen)
  }
  cpt_YI <- optimize(ff, interval = range(tt), tt = tt, D = D_IPW, D0 = D0_IPW, pr = pr.dise_IPW, pr0 = pr.dise0_IPW,
                     maximum = TRUE)$maximum
  Spe_YI <- mean((tt <= cpt_YI)*D0_IPW)/pr.dise0_IPW
  Sen_YI <- mean((tt > cpt_YI)*D_IPW)/pr.dise_IPW
  plr <- Sen_YI/(1 - Spe_YI)
  nlr <- (1 - Sen_YI)/Spe_YI
  ppv <- Sen_YI*prev/(Sen_YI*prev + (1 - Spe_YI)*(1 - prev))
  npv <- Spe_YI*(1 - prev)/((1 - Sen_YI)*prev + Spe_YI*(1 - prev))
  return(c(cpt_YI = cpt_YI, Sen_YI = Sen_YI, Spe_YI = Spe_YI, plr = plr, nlr = nlr, ppv = ppv, npv = npv))
}

#' @export
Sen_Spe_YI_IPW_boot <- function(data, name.test, name.dise, name.veri, B, formula.veri, prev){
  data_split <- split(data, f = as.factor(data[, name.dise]))
  data_1 <- data_split[[1]]
  data_2 <- data_split[[2]]
  data_3 <- data_split[[3]]
  n1 <- nrow(data_1)
  n2 <- nrow(data_2)
  n3 <- nrow(data_3)
  bts_fun <- function(data_1, data_2, data_3, n1, n2, n3, formula.veri, name.test, name.veri, prev){
    id1 <- sample(1:n1, n1, replace = TRUE)
    id2 <- sample(1:n2, n2, replace = TRUE)
    id3 <- sample(1:n3, n3, replace = TRUE)
    data_1_b <- data_1[id1,]
    data_2_b <- data_2[id2,]
    data_3_b <- data_3[id3,]
    data_b <- rbind(data_1_b, data_2_b, data_3_b)
    out_veri_model_b <- veri_model(formula.veri, data = data_b, name.veri = name.veri, print = FALSE)
    out <- Sen_Spe_YI_IPW(out_veri_model_b, tt = data_b[, name.test], dd = data_b[, name.dise], 
                          vv = data_b[, name.veri], prev = prev)
    return(out)
  }
  out_bts <- sapply(1:B, function(i){
    bts_fun(data_1 = data_1, data_2 = data_2, data_3 = data_3, n1 = n1, n2 = n2, n3 = n3,
            formula.veri = formula.veri, name.test = name.test, name.veri = name.veri, prev = prev)
  })
  return(out_bts)
}

#' @export
Sen_Spe_YI_DR <- function(out_dise_model, out_veri_model, tt, dd, vv, prev){
  rho_est <- out_dise_model$rho_est
  pi_est <- out_veri_model$pi_est
  D_DR <- vv*dd/pi_est - (vv - pi_est)*rho_est/pi_est
  pr.dise <- mean(D_DR)
  # cpt_full <- c(-Inf, sort(tt), Inf)
  # Spe <- sapply(cpt_full, function(x) mean((tt <= x)*D0_IPW)/pr.dise0_IPW)
  # Sen <- sapply(cpt_full, function(x) mean((tt > x)*D_IPW)/pr.dise_IPW)
  # YI <- Spe + Sen
  # cpt_YI <- mean(cpt_full[which(YI == max(YI))])
  ff <- function(x, tt, D, pr){
    Spe <- mean((tt <= x)*(1 - D))/(1 - pr)
    Sen <- mean((tt > x)*D)/pr
    return(Spe + Sen)
  }
  cpt_YI <- optimize(ff, interval = range(tt), tt = tt, D = D_DR, pr = pr.dise, maximum = TRUE)$maximum
  Spe_YI <- mean((tt <= cpt_YI)*(1 - D_DR))/(1 - pr.dise)
  Sen_YI <- mean((tt > cpt_YI)*D_DR)/pr.dise
  plr <- Sen_YI/(1 - Spe_YI)
  nlr <- (1 - Sen_YI)/Spe_YI
  ppv <- Sen_YI*prev/(Sen_YI*prev + (1 - Spe_YI)*(1 - prev))
  npv <- Spe_YI*(1 - prev)/((1 - Sen_YI)*prev + Spe_YI*(1 - prev))
  return(c(cpt_YI = cpt_YI, Sen_YI = Sen_YI, Spe_YI = Spe_YI, plr = plr, nlr = nlr, ppv = ppv, npv = npv))
}

#' @export
Sen_Spe_YI_DR_boot <- function(data, name.test, name.dise, name.veri, B, formula.dise, formula.veri, prev){
  data_split <- split(data, f = as.factor(data[, name.dise]))
  data_1 <- data_split[[1]]
  data_2 <- data_split[[2]]
  data_3 <- data_split[[3]]
  n1 <- nrow(data_1)
  n2 <- nrow(data_2)
  n3 <- nrow(data_3)
  bts_fun <- function(data_1, data_2, data_3, n1, n2, n3, formula.dise, formula.veri, name.test, name.veri, prev){
    id1 <- sample(1:n1, n1, replace = TRUE)
    id2 <- sample(1:n2, n2, replace = TRUE)
    id3 <- sample(1:n3, n3, replace = TRUE)
    data_1_b <- data_1[id1,]
    data_2_b <- data_2[id2,]
    data_3_b <- data_3[id3,]
    data_b <- rbind(data_1_b, data_2_b, data_3_b)
    out_dise_model_b <- dise_model(formula.dise, data = data_b, name.veri = name.veri, print = FALSE)
    out_veri_model_b <- veri_model(formula.veri, data = data_b, name.veri = name.veri, print = FALSE)
    out <- Sen_Spe_YI_DR(out_dise_model_b, out_veri_model_b, tt = data_b[, name.test], dd = data_b[, name.dise], 
                         vv = data_b[, name.veri], prev = prev)
    return(out)
  }
  out_bts <- sapply(1:B, function(i){
    bts_fun(data_1 = data_1, data_2 = data_2, data_3 = data_3, n1 = n1, n2 = n2, n3 = n3, formula.dise = formula.dise,
            formula.veri = formula.veri, name.test = name.test, name.veri = name.veri, prev = prev)
  })
  return(out_bts)
}

