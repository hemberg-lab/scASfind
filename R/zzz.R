#' Essential File so module is loaded
#'
#' @importFrom Rcpp loadModule
#' @useDynLib scASfind
#'
#' @include zzz.R
#' @name EliasFanoDB
#' @export
NULL

loadModule("EliasFanoDB", TRUE)
