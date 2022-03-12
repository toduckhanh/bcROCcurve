#' @useDynLib bcROCcurve
#' @importFrom Rcpp evalCpp

#' @export
aucFull <- function(T, D){
  return(aucFull_C(T, D))
}

#' @export
aucFI <- function(T, rhoEst){
  return(aucFI_C(T, rhoEst))
}

#' @export
aucMSI <- function(T, D, V, rhoEst){
  return(aucMSI_C(T, D, V, rhoEst))
}

#' @export
aucIPW <- function(T, D, V, piEst){
  return(aucIPW_C(T, D, V, piEst))
}

#' @export
aucPDR <- function(T, D, V, rhoEst, piEst){
  return(aucPDR_C(T, D, V, rhoEst, piEst))
}

#' @export
aucMAR <- function(T, D, V, rhoEst, piEst){
  return(aucMAR_C(T, D, V, rhoEst, piEst))
}

