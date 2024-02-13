
#' Conditional covariance estimation by kernel-smoothing
#'
#' Computed an estimator of the
#' conditional covariance matrix of Y given X = x
#'
#' @param datax vector or matrix of observations of the conditioning variable
#' @param datay vector or matrix of observations of the conditioned variable,
#' with the same number of observations as `datax`
#'
#' @param x0 point at which the conditional covariance matrix is estimated.
#' @param condmean known or estimated conditional mean of \eqn{(Y_1,..., Y_d)}
#'
#' @examples
#' n = 10
#' X <- 1:n
#' mu_1 <- X^2
#' mu_2 <- 3*X
#' sd_1 <- X
#' sd_2 <- 1 + X
#' cor <- sin(X / 20)
#'
#' # Computation of the covariance matrices
#' cov <- cor * sd_1 * sd_2
#' Sigma = array(dim = c(2,2,n))
#' Sigma[1,1,] <- sd_1^2
#' Sigma[2,2,] <- sd_2^2
#' Sigma[1,2,] <- Sigma[2,1,] <- cov
#'
#' SqrtSigma = apply(X = Sigma, MARGIN = c(3), FUN = chol)
#' SqrtSigma = array(SqrtSigma, dim = c(2,2,n))
#' Z = matrix(rnorm(n = 2*n), ncol = 2)
#'
#' Y1 = mu_1 + Z[,1] * SqrtSigma[1,1,] + Z[,2] * SqrtSigma[1,2,]
#' Y2 = mu_2 + Z[,1] * SqrtSigma[2,1,] + Z[,2] * SqrtSigma[2,2,]
#' datay = cbind(Y1, Y2)
#'
#' # Comparison with the usual (unconditional) covariance
#' condcov = cond.cov(datax = X, datay,
#'                    x0 = 4, h = 1000)
#' print(condcov)
#' print(cov(datay))
#'
#' # We now do smoothing
#' condcov = cond.cov(datax = X, datay,
#'                    x0 = 4, h = 2)
#' print(condcov)
#' print(cov(datay))
#'
#'
# # Test
# datax = c(2,5)
# datay = cbind(c(1,6), c(2,7))
# x0 = 1
# h = 1
# cond.cov(datax = datax, datay = datay, x0 = x0, h = 1, condmean = c(0,0))
# ( (1-0)*(2-0) * exp(-(2 - 1)^2) + (6-0)*(7-0) * exp(-(5 - 1)^2)
# ) / (exp(-(2 - 1)^2) + exp(-(5 - 1)^2))
#'
#' @export
#'
cond.cov <- function(datax, datay, x0, h, condmean = NULL)
{
  # Test of compatibility for the dimensions
  if(NROW(datax) != NROW(datay)){
    stop("The number of observations in datax and datay should be equal.")
  }
  if(NCOL(datax) != length(x0)){
    stop("x0 should have a number of entries equal to the number of columns in datax.")
  }

  # We convert vector (if given) to matrices with one column
  if (is.vector(datax) ){
    datax = matrix(datax, ncol = 1)
  }
  if (is.vector(datay) ){
    datay = matrix(datay, ncol = 1)
  }

  # Centering by the mean
  if (is.null(condmean)){
    result = cond.mean.internal(datax = datax, datay = datay, x0 = x0, h = h)
    condmean = result$condmean
    weights_ = result$weights_
  } else {
    if(NCOL(datay) != length(condmean)){
      stop("condmean should have a length equal to the number of columns in datay.")
    }

    scaled_x = abs( sweep(x = datax, MARGIN = 2, STATS = x0) / h )
    weights_ = exp( - rowSums(scaled_x^2))
    weights_ = weights_ / sum(weights_)
  }

  normalizedy = sweep(x = datay, MARGIN = 2, STATS = condmean)

  condcov = matrix(ncol = NCOL(datay), nrow = NCOL(datay))

  for (i in 1:NCOL(datay)){
    for (j in 1:i){

      condcov[i,j] = sum( weights_ * normalizedy[,i] * normalizedy[,j] ) / sum(weights_)

      condcov[j,i] <- condcov[i,j]
    }
  }

  return(condcov)
}


