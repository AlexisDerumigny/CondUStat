
#' Conditional mean estimation by kernel smoothing
#'
#' @param datax the vector of observed X values.
#' @param datay the vector of observed Y values. Must have the same length as \code{datax}
#' @param x0 the point at which the U-statistic should be computed.
#' @param h the bandwidth of the kernel smoothing.
#'
#' @examples
#' n = 500
#' datax = rnorm(n)
#' datay = datax + rnorm(n)
#' x0 = 1
#' h = 0.2
#' cond.mean(datax = datax, datay = datay, x0 = x0, h = h)
#'
#' @export
cond.mean <- function(datax, datay, x0, h)
{
  result = cond.UStat(
    datax = datax, datay = datay,
    FUN = function(y1){return(y1)}, p=1,
    x0 = x0, h = h)

  return(result)
}


#' Conditional variance estimation by kernel smoothing
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
#' cond.var(datax = datax, datay = datay, x0 = x0, h = h)
#'
#' @export
cond.var <- function(datax, datay, x0, h, condmean = NULL)
{
  if (is.null(condmean)){
    condmean = cond.mean(datax = datax, datay = datay, x0 = x0, h = h)
  }

  result = cond.UStat(
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

  moment3 = cond.UStat(
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

  moment4 = cond.UStat(
    datax = datax, datay = datay,
    FUN = function(y1){return( (y1 - condmean)^4 ) }, p=1,
    x0 = x0, h = h)

  # moment4 = cond.UStat(
  #   datax = datax, datay = datay,
  #   FUN = function(y1, y2, y3, y4){
  #     return(y1^4 - 4 * y1^3 * y2 + 6 * y1^2 * y2 * y3 - 3 * y1 * y2 * y3 * y4)},
  #   p=4,
  #   x0 = c(x0,x0,x0,x0), h = h)

  condvar = cond.var(datax = datax, datay = datay, x0 = x0, h = h)

  result = moment4 / condvar^2

  return(result)
}
