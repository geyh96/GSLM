#' @useDynLib GSLM, .registration = TRUE
#' @importFrom methods as
#' @importFrom MASS ginv
#' @importFrom Matrix solve
#' @importFrom kernlab rbfdot
#' @importFrom kernlab sigest
#' @importFrom testit assert
#' @importFrom stats median
#' @importFrom kernlab kernelMatrix
#' @export
NULL

#' GSLM
#' 
#' @description A package for Structure learning via unstructured kernel-based M-regression in Reproducing Kernel Hilbert Space

#' \tabular{ll}{
#'   Package: \tab GSLM\cr
#'   Type: \tab Package\cr
#'   Version: \tab 1.0.0\cr
#'   Date: \tab 2021-04-28\cr
#'   License: \tab GPL-2\cr
#'   LazyLoad: \tab yes\cr
#' }
#' @docType package
#' @aliases GSLM-package
#' @author Xin He, Yeheng Ge and Xingdong Feng \cr
#' Maintainers: Yeheng Ge <geyh96@163.sufe.edu.cn>;
#' @references
#'  The meanuscripts 
"_PACKAGE"
#> [1] "_PACKAGE"
