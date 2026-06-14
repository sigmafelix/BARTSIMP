#' Matérn Covariance Matrix
#'
#' Computes a dense Matérn covariance matrix for spatial coordinates.
#'
#' @param nu Smoothness parameter of the Matérn covariance.
#' @param kappa Inverse range parameter.
#' @param coords Matrix or data frame with one row per spatial location.
#' @param sigma2e Nugget variance added to the diagonal.
#' @param sigma2x Marginal spatial variance.
#'
#' @return A dense covariance matrix with one row and column per coordinate.
#'
#' @examples
#' coords <- data.frame(s1 = c(0, 0.5, 1), s2 = c(0, 0.25, 1))
#' matern.mat(nu = 1, kappa = 2, coords = coords, sigma2e = 0.1, sigma2x = 1)
#'
#' @export
matern.mat <- function(nu, kappa, coords, sigma2e, sigma2x) {
  n <- dim(coords)[1]
  dmat <- dist(coords)
  mcor <- as.matrix(2^(1-nu)*(kappa*dmat)^nu*
                      besselK(dmat*kappa,nu)/gamma(nu))
  diag(mcor) <- 1
  mcov <- sigma2e*diag(n) + sigma2x*mcor
  return(mcov)
}
