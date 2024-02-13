
#' Conditional mean estimation by kernel smoothing
#'
#' @param datax vector or matrix of observed X values (conditioning variables).
#' @param datay vector or matrix of observed Y values (conditioned variables).
#' Must have the same number of observations as \code{datax}
#'
#' @param x0 point at which the U-statistic should be computed.
#' @param h the bandwidth of the kernel smoothing.
#'
#' @return a vector having `NROW(datay)` entries with the conditional mean
#' of the corresponding variables.
#'
#' @examples
#' n = 500
#' datax = rnorm(n)
#' datay = datax + rnorm(n)
#' x0 = 1
#' h = 0.2
#' cond.mean(datax = datax, datay = datay, x0 = x0, h = h)
#'
# # Test
# datax = c(2,5)
# datay = c(1,6)
# x0 = 1
# h = 1
# cond.mean(datax = datax, datay = datay, x0 = x0, h = 1)
# ( 1 * exp(-(2 - 1)^2) + 6 * exp(-(5 - 1)^2) )/ (exp(-(2 - 1)^2) + exp(-(5 - 1)^2))
#'
#'
#' @export
cond.mean <- function(datax, datay, x0, h)
{
  if (is.vector(datax) & is.vector(datay)) {
    result = cond.UStat.univ(
      datax = datax, datay = datay,
      FUN = function(y1){return(y1)}, p=1,
      x0 = x0, h = h, cutoff = Inf)
  }

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

  result = cond.mean.internal(datax, datay, x0, h)$condmean
  return(result)
}

#' Function for computing the conditional mean
#'
#' Internal function, assume that all inputs satisfy the required conditions
#'
#' @param datax matrix of observations of x
#' @param datay matrix of observations of y
#'
#' @return a list whose elements are: \itemize{
#'   \item `condmean`: the estimated conditional mean
#'   \item `weights_`: the weights used (a vector of the same length as `ncol(Y)`)
#' }
#'
#' @examples
#' datax = cbind(1:5, 101:105)
#' datay = cbind(4+rnorm(5), 12 + rnorm(5))
#' x0 = c(2,103)
#' h = 1
#'
#' result = cond.mean.internal(datax, datay, x0, h)
#' weighted.mean(datay[,1], weights_) # Testing
#'
#' @noRd
cond.mean.internal <- function(datax, datay, x0, h)
{
  if (is.matrix(datax)){
    scaled_x = abs( sweep(x = datax, MARGIN = 2, STATS = x0) / h )
    weights_ = exp(- rowSums(scaled_x^2))
  } else {
    scaled_x = abs( (datax - x0) / h )
    weights_ = exp(- scaled_x^2)
  }
  weights_ = weights_ / sum(weights_)
  result = matrix(weights_, nrow = 1) %*% datay

  return(list(condmean = result, weights_ = weights_))
}


#' Conditional variance estimation by kernel smoothing
#'
#' @param datax the vector of observed X values (conditioning variable).
#' @param datay the vector of observed Y values (conditioned variable).
#' Must have the same length as \code{datax}
#' @param x0 the point at which the U-statistic should be computed.
#' @param h the bandwidth of the kernel smoothing.
#' @param condmean value of the conditional mean of y at x=x0
#' (if known or estimated previously).
#'
#' @examples
#' n = 500
#' datax = rnorm(n)
#' datay = datax + rnorm(n)
#' x0 = 1
#' h = 0.2
#' cond.var(datax = datax, datay = datay, x0 = x0, h = h)
#'
#' @export
cond.var <- function(datax, datay, x0, h, condmean = NULL)
{
  if (is.null(condmean)){
    condmean = cond.mean(datax = datax, datay = datay, x0 = x0, h = h)
  }

  result = cond.UStat.univ(
    datax = datax, datay = datay,
    FUN = function(y1){return( (y1 - condmean)^2 ) }, p=1,
    x0 = x0, h = h)

  # Equivalent possibility, but longer:
  # result = cond.UStat(
  #   datax = datax, datay = datay,
  #   FUN = function(y1, y2){return(y1^2 - y1 * y2)}, p=2,
  #   x0 = c(x0,x0), h = h)

  return(result)
}


#' Conditional skewness estimation by kernel smoothing
#'
#' @param datax the vector of observed X values.
#' @param datay the vector of observed Y values. Must have the same length as \code{datax}
#' @param x0 the point at which the U-statistic should be computed.
#' @param h the bandwidth of the kernel smoothing.
#' @param condmean value of the conditional mean (if know or estimated previously).
#'
#' @examples
#' n = 500
#' datax = rnorm(n)
#' datay = datax + rnorm(n)
#' x0 = 1
#' h = 0.3
#' cond.skewness(datax = datax, datay = datay, x0 = x0, h = h)
#'
#' @export
cond.skewness <- function(datax, datay, x0, h, condmean = NULL)
{
  if (is.null(condmean)){
    condmean = cond.mean(datax = datax, datay = datay, x0 = x0, h = h)
  }

  moment3 = cond.UStat.univ(
    datax = datax, datay = datay,
    FUN = function(y1){return( (y1 - condmean)^3 ) }, p=1,
    x0 = x0, h = h)

  # moment3 = cond.UStat(
  #   datax = datax, datay = datay,
  #   FUN = function(y1, y2, y3){return(y1^3 - 3 * y1^2 * y2 + 2 * y1 * y2 * y3)}, p=3,
  #   x0 = c(x0,x0,x0), h = h)

  condvar = cond.var(datax = datax, datay = datay, x0 = x0, h = h, condmean = condmean)

  result = moment3 / condvar^(3/2)

  return(result)
}


#' Conditional kurtosis estimation by kernel smoothing
#'
#' @param datax the vector of observed X values.
#' @param datay the vector of observed Y values. Must have the same length as \code{datax}
#' @param x0 the point at which the U-statistic should be computed.
#' @param h the bandwidth of the kernel smoothing.
#' @param condmean value of the conditional mean (if know or estimated previously).
#'
#' @examples
#' n = 500
#' datax = rnorm(n)
#' datay = datax + rnorm(n)
#' x0 = 1
#' h = 0.2
#' cond.kurtosis(datax = datax, datay = datay, x0 = x0, h = h)
#'
#' @export
cond.kurtosis <- function(datax, datay, x0, h, condmean = NULL)
{
  if (is.null(condmean)){
    condmean = cond.mean(datax = datax, datay = datay, x0 = x0, h = h)
  }

  moment4 = cond.UStat.univ(
    datax = datax, datay = datay,
    FUN = function(y1){return( (y1 - condmean)^4 ) }, p=1,
    x0 = x0, h = h)

  # moment4 = cond.UStat.univ(
  #   datax = datax, datay = datay,
  #   FUN = function(y1, y2, y3, y4){
  #     return(y1^4 - 4 * y1^3 * y2 + 6 * y1^2 * y2 * y3 - 3 * y1 * y2 * y3 * y4)},
  #   p=4,
  #   x0 = c(x0,x0,x0,x0), h = h)

  condvar = cond.var(datax = datax, datay = datay, x0 = x0, h = h)

  result = moment4 / condvar^2

  return(result)
}
