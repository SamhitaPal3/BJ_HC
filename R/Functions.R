## BH-FDR 
#' Title
#'
#' @param p_val 
#' @param alpha 
#'
#' @return multiplicity corrected p-value
#' @export
#'
#' @examples
BH_FDR <- function(p_val, alpha){
  N <- length(p_val)   
  sortp <- sort(p_val) 
  unifp <- (1 : N) / N *alpha
  ind = which(sortp<unifp)
  if (length(ind) ==0) {
    tt = 0
  } else {
    tt = max(ind)
  }  
  return(sortp[tt])
}

