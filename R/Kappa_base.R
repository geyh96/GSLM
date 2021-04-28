#' Calculate the kappa index proposed in the Sun W, Wang J, Fang Y. Consistent selection of tuning parameters via variable selection stability. The Journal of Machine Learning Research, 2013, 14(1): 3419-3440.
#' @param result1 a vector contains the variable selected from the sub dataset 1
#' @param result2 a vector contains the variable selected from the sub dataset 2
#' @return The calculate Kappa index from result1 and result2
#' @export
calculate <- function(result1, result2) {
  cou.p <- length(result1)
  n11 <- 0
  n12 <- 0
  n21 <- 0
  n22 <- 0
  Cal.Kapp <- 0
  for (i in 1:cou.p) {
    if (result1[i] == 1 & result2[i] == 1) {
      n11 <- n11 + 1
    } else if (result1[i] == 1 & result2[i] == 0) {
      n12 <- n12 + 1
    } else if (result1[i] == 0 & result2[i] == 1) {
      n21 <- n21 + 1
    }
    else if (result1[i] == 0 & result2[i] == 0) {
      n22 <- n22 + 1
    }
  }
  if (n22 == cou.p) {
    Cal.Kapp <- -1
  } else if (n11 == cou.p) {
    Cal.Kapp <- -1
  }
  else {
    pra <- (n11 + n22) / cou.p
    pre <- (n11 + n12) * (n11 + n21) / (cou.p^2) + (n12 + n22) * (n21 + n22) / (cou.p^2)
    Cal.Kapp <- (pra - pre) / (1 - pre)
  }
  return(Cal.Kapp)
}


#' calculate the kappa index proposed in the Sun(2013). the kap.lam with largest return value are selected.
#' @param kap.re1 a matrix of b by p . b the the time of simulation and p is the dimension of covariates . kap.re1[i,] contains the variable selected from the sub dataset 1 in the i'th loop
#' @param kap.re2  a matrix of b by p . b the the time of simulation and p is the dimension of covariates . kap.re1[i,] contains the variable selected from the sub dataset 2 in the i'th loop
#' @param kap.lam the candidate of threshold used in the GSLM framework to get exact variable selection.
#' @return The calculate Kappa index corresponding to the kap.lam
#' @export
Kapp.cal <- function(kap.re1, kap.re2, kap.lam) {
  kap.value <- 0
  kap.re1[(abs(kap.re1) < kap.lam)] <- 0
  kap.re1[(abs(kap.re1) >= kap.lam)] <- 1
  kap.re2[(abs(kap.re2) < kap.lam)] <- 0
  kap.re2[(abs(kap.re2) >= kap.lam)] <- 1
  for (k1 in 1:nrow(kap.re1)) {
    kap.value <- kap.value + calculate(kap.re1[k1, ], kap.re2[k1, ])
  }
  return(kap.value)
}