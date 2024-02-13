

#' Computation of a univariate conditional U-statistics by kernel smoothing
#'
#'
#' @param datax the vector of observed X values.
#' @param datay the vector of observed Y values. Must have the same length as \code{datax}
#' @param FUN the kernel of the U statistic.
#' Should be vectorized and take as parameters \code{y1}, ..., \code{yp}.
#' (ie \code{p} vectors of the same length)
#' @param p the degree of the U-statistic given by FUN.
#' @param x0 the point at which the U-statistic should be computed.
#' @param h the bandwidth of the kernel smoothing.
#'
#' @param cutoff value above which the centered observations (X - x0)/h
#' have no influence.
#'
#' @examples
#' n = 500
#' datax = rnorm(n)
#' datay = datax + rnorm(n)
#' x0 = 1
#' h = 0.2
#' # Estimation of a conditional mean
#' cond.UStat(datax, datay, FUN = function(y1){return(y1)}, p=1, x0, h)
#' # Naive estimation of a conditional variance
#' cond.UStat(
#'   datax = datax, datay = datay,
#'   FUN = function(y1, y2){return(y1^2 - y1 * y2)}, p=2,
#'   x0 = c(x0,x0), h = h)
#'
#' @export
cond.UStat.univ <- function(datax, datay, FUN, p, x0, h, cutoff = 8)
{
  stopifnot(NROW(datax) == NROW(datay),
            length(x0) == p)

  switch(
    p,

    {
      # p = 1
      scaled_x = abs( (datax - x0) / h )
      relevant_x = which( scaled_x < cutoff )
      reduced_x = scaled_x[relevant_x]
      reduced_y = datay[relevant_x]

      weights_ = exp( - reduced_x^2)
      results_FUN = FUN(y1 = reduced_y)
      result = sum( weights_ * results_FUN) / sum(weights_)
    },

    {
      # p = 2
      scaled_x1 = abs( (datax - x0[1]) / h )
      scaled_x2 = abs( (datax - x0[2]) / h )
      
      relevant_x1 = which( scaled_x1 < cutoff )
      reduced_x1 = scaled_x1[relevant_x1]
      reduced_y1 = datay[relevant_x1]
      
      relevant_x2 = which( scaled_x2 < cutoff )
      reduced_x2 = scaled_x2[relevant_x2]
      reduced_y2 = datay[relevant_x2]

      expandX = expand.grid(reduced_x1, reduced_x2)
      expandY = expand.grid(reduced_y1, reduced_y2)

      weights_ = exp( - (expandX$Var1^2 + expandX$Var2^2) )
      results_FUN = FUN(y1 = expandY$Var1, y2 = expandY$Var2)
      result = sum( weights_ * results_FUN) / sum(weights_)
    },

    {
      # p = 3
      scaled_x1 = abs( (datax - x0[1]) / h )
      scaled_x2 = abs( (datax - x0[2]) / h )
      scaled_x3 = abs( (datax - x0[3]) / h )
      relevant_x1 = which( scaled_x1 < cutoff )
      reduced_x1 = scaled_x[relevant_x1]
      reduced_y1 = datay[relevant_x1]
      relevant_x2 = which( scaled_x2 < cutoff )
      reduced_x2 = scaled_x[relevant_x2]
      reduced_y2 = datay[relevant_x2]
      relevant_x3 = which( scaled_x3 < cutoff )
      reduced_x3 = datax[relevant_x3]
      reduced_y3 = datay[relevant_x3]

      expandX = expand.grid(x1 = reduced_x1, x2 = reduced_x2, x3 = reduced_x3)
      expandY = expand.grid(y1 = reduced_y1, y2 = reduced_y2, y3 = reduced_y3)

      weights_ = exp( - (expandX$x1^2 + expandX$x2^2 + expandX$x3^2) )
      results_FUN = FUN(y1 = expandY$y1, y2 = expandY$y2, y3 = expandY$y3)
      result = sum( weights_ * results_FUN) / sum(weights_)
    },

    {
      # p = 4
      scaled_x1 = abs( (datax - x0[1]) / h )
      scaled_x2 = abs( (datax - x0[2]) / h )
      scaled_x3 = abs( (datax - x0[3]) / h )
      scaled_x4 = abs( (datax - x0[4]) / h )
      relevant_x1 = which( scaled_x1 < cutoff )
      reduced_x1 = scaled_x[relevant_x1]
      reduced_y1 = datay[relevant_x1]
      relevant_x2 = which( scaled_x2 < cutoff )
      reduced_x2 = scaled_x[relevant_x2]
      reduced_y2 = datay[relevant_x2]
      relevant_x3 = which( scaled_x3 < cutoff )
      reduced_x3 = scaled_x[relevant_x3]
      reduced_y3 = datay[relevant_x3]
      relevant_x4 = which( scaled_x4 < cutoff )
      reduced_x4 = scaled_x[relevant_x4]
      reduced_y4 = datay[relevant_x4]

      expandX = expand.grid(
        x1 = reduced_x1, x2 = reduced_x2, x3 = reduced_x3, x4 = reduced_x4)
      expandY = expand.grid(
        y1 = reduced_y1, y2 = reduced_y2, y3 = reduced_y3, y4 = reduced_y4)

      weights_ = exp( - (expandX$x1^2 + expandX$x2^2 + expandX$x3^2 + expandX$x4^2) )
      results_FUN = FUN(y1 = expandY$y1, y2 = expandY$y2, y3 = expandY$y3, y4 = expandY$y4)
      result = sum( weights_ * results_FUN) / sum(weights_)
    },

    {
      stop("p>4 not implemented in cond.UStat.univ")
    }
  )

  return(result)
}
