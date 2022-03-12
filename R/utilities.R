
expit <- function (x) {
  return(as.numeric(1/(1 + exp(-x))))
}


design_Matrix <- function(formula, data){
  ff <- formula
  m <- model.frame(ff, data, na.action = NULL)
  X <- model.matrix(ff, m)
  return(X)
}

#' @export
dise_model <- function(formula, data, name.veri, print = TRUE){
  out_dise <- glm(formula, data = data[data[, name.veri] == 1, ], family = binomial(link = "logit"))
  D <- model.frame(out_dise$terms, data = data)[,1]
  V <- data[, name.veri]
  rho_est <- predict(out_dise, newdata = data, type = "response")
  mat.dise <- design_Matrix(out_dise$terms, data = data)
  rho_score <- V*(D - rho_est)*mat.dise
  rho_hess <- solve(summary(out_dise)$cov.scaled)
  npar.dise <- ncol(mat.dise)
  if(print) print(summary(out_dise))
  return(list(out_dise = out_dise, rho_est = rho_est, rho_score = rho_score, rho_hess = rho_hess, 
              mat.dise = mat.dise, npar.dise = npar.dise))
}

#' @export
veri_model <- function(formula, data, name.veri, print = TRUE){
  out_veri <- glm(formula, data = data, family = binomial(link = "logit"))
  pi_est <- predict(out_veri, type = "response")
  pi_hess <- solve(summary(out_veri)$cov.scaled)
  mat.veri <- design_Matrix(out_veri$terms, data = data)
  npar.veri <- ncol(mat.veri)
  pi_score <- (data[,name.veri] - pi_est)*mat.veri
  if(print) print(summary(out_veri))
  return(list(out_veri = out_veri, pi_est = pi_est, pi_score = pi_score, pi_hess = pi_hess, 
              mat.veri = mat.veri, npar.veri = npar.veri))
}
