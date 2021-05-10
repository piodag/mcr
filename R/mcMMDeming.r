###############################################################################
##
## mcWDeming.r
##
## Function for computing weighted deming regression for two methods with  proportional errors.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
## Copyright (C) 2020 Giorgio Pioda for the MM-Deming
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Calculate Weighted Deming Regression
#'
#' Calculate weighted deming regression with iterative algorithm suggested by Linnet.
#' This algorithm is avalaible only for positive values. But even in this case there is no guarantee that
#' the algorithm always converges.
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param error.ratio ratio between squared measurement errors of reference- and test method, necessary for Deming regression (Default is 1).
#' @param iter.max maximal number of iterations.
#' @param threshold threshold value.
#' @return a list with elements
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{xw}{average of reference method values.}
#'  \item{iter}{number of iterations.}
#' @references  Linnet K.
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              CLIN. CHEM. 39/3, 424-432 (1993).
#'
#'              Linnet K.
#'              Estimation of the Linear Relationship between the Measurements of two Methods with Proportional Errors.
#'              STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).
mc.mmdemingConstCV <- function(X, Y, error.ratio, iter.max=30, threshold=0.000001)
{
  # Check validity of parameters

  stopifnot(is.numeric(X))
  stopifnot(is.numeric(Y))
  stopifnot(length(X) == length(Y))
  stopifnot(is.numeric(error.ratio))
  stopifnot(error.ratio > 0)
  stopifnot(is.numeric(iter.max))
  stopifnot(round(iter.max) == iter.max)
  stopifnot(iter.max > 0)
  stopifnot(is.numeric(threshold))
  stopifnot(threshold >= 0)

  # This algorithm often doesn't converge if there are negative
  # measurements in data set

#  if (min(X)<0 | min(Y)<0)
#  {
#    return(paste("There are negative values in data set."))
#  }
#  else
#  {
    # 1. Calculate  initial values
    #    (point estimates of unweighted Deming regression)

    n <- length(X)

    if (n >= 100){
      start.n<-20
    } else if(n >= 76){
      start.n<-50
    } else if (n >= 46){
      start.n<-100
    } else if (n >= 36) {
      start.n<-250
    } else {
      #message("there is no convergence warranty below 36 samples")
      start.n<-250
    }

    #mX <- mean(X)
    #mY <- mean(Y)
    #u <- sum((X-mX)^2)
    #q <- sum((Y-mY)^2)
    #p <- sum((X-mX)*(Y-mY))

    ## initial values, b1 as the straight bisector of the two slopes obtained by the robust covariance matrix
    ## and b0 as the derived intercept using its centers. The scale parameter is then obtained with the
    ## mad() of the euclidean calculated residuals.

    #b1 <- ((error.ratio*q-u)+sqrt((u-error.ratio*q)^2+4*error.ratio*p^2))/(2*error.ratio*p)
    #b0 <- mean(Y)-b1*mean(X)

    cov.sest<-rrcov::CovSest(cbind(X,Y),nsamp=start.n,bdp=0.5)
    b1<-mean(c(cov.sest$cov[2,1]/cov.sest$cov[1,1],1/(cov.sest$cov[2,1]/cov.sest$cov[2,2])))
    b0<-cov.sest$center[2]-b1*cov.sest$center[1]

    d <- Y-(b0+b1*X)
    XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
    YHAT <- Y-(d/(1+error.ratio*b1^2))
    euclid.d<-sqrt((X-XHAT)^2+(Y-YHAT)^2)
    scale.lts <- mad(euclid.d)
    #stdeuclid.d<-euclid.d/scale.lts
    #k<-4.685
    #W<-ifelse(abs(stdeuclid.d) <= k, (1-(stdeuclid.d/k)^2)^2, 0)
    #message(paste("start b1:",b1,"start b0:",b0,"scale:",scale.lts,"sum W:",sum(W)))
    error.count<-0


    ## Iterative Algorithm

    i <- 0     # Number of iterations
    warn<-"no warnings"   # Warnings

    repeat
    {
      if (i >= iter.max)
      {
        warning(paste("no convergence after",iter.max,"iterations, last values B1:",B1,"B0:",B0))
        break
      }

      i<-i+1

      # Calculation of weights
      d <- Y-(b0+b1*X)
      XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
      YHAT <- Y-(d/(1+error.ratio*b1^2))
      #W <- ((XHAT+error.ratio*YHAT)/(1+error.ratio))^(-2)

      euclid.d<-sqrt((X-XHAT)^2+(Y-YHAT)^2)
      #euclid.mad<-mad(euclid.d)
      stdeuclid.d<-euclid.d/scale.lts
      k<-4.685
      W<-ifelse(abs(stdeuclid.d) <= k, (1-(stdeuclid.d/k)^2)^2, 0)

      # Calculation of regression coefficients
      XW <- sum(W*X)/sum(W)
      YW <- sum(W*Y)/sum(W)
      U <- sum(W*((X-XW)^2))
      Q <- sum(W*((Y-YW)^2))
      P <- sum(W*(X-XW)*(Y-YW))

      # Point estimates
      B1 <- (error.ratio*Q-U+sqrt((U-error.ratio*Q)^2+4*error.ratio*P^2))/(2*error.ratio*P)

      # If a singularity arises during covariance sfast calculation an alternative Rocke restart is attempted first

      if (!is.finite(B1)){

          if (error.count==0) {
             #message(paste("Rocke prior restart","b1:",b1,"b0:",b0,"sum W:",sum(W)))
             error.count<-error.count+1
             cov.sest<-rrcov::CovSest(cbind(X,Y),method="rocke")
             b1<-mean(c(cov.sest$cov[2,1]/cov.sest$cov[1,1],1/(cov.sest$cov[2,1]/cov.sest$cov[2,2])))
             b0<-cov.sest$center[2]-b1*cov.sest$center[1]

             d <- Y-(b0+b1*X)
             XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
             YHAT <- Y-(d/(1+error.ratio*b1^2))
             euclid.d<-sqrt((X-XHAT)^2+(Y-YHAT)^2)
             scale.lts <- mad(euclid.d)

             # Calculation of weights
             d <- Y-(b0+b1*X)
             XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
             YHAT <- Y-(d/(1+error.ratio*b1^2))
             #W <- ((XHAT+error.ratio*YHAT)/(1+error.ratio))^(-2)

             euclid.d<-sqrt((X-XHAT)^2+(Y-YHAT)^2)
             #euclid.mad<-mad(euclid.d)
             stdeuclid.d<-euclid.d/scale.lts
             k<-4.685
             W<-ifelse(abs(stdeuclid.d) <= k, (1-(stdeuclid.d/k)^2)^2, 0)

             # Calculation of regression coefficients
             XW <- sum(W*X)/sum(W)
             YW <- sum(W*Y)/sum(W)
             U <- sum(W*((X-XW)^2))
             Q <- sum(W*((Y-YW)^2))
             P <- sum(W*(X-XW)*(Y-YW))

             # Point estimates

             B1 <- (error.ratio*Q-U+sqrt((U-error.ratio*Q)^2+4*error.ratio*P^2))/(2*error.ratio*P)
             warning(paste("Rocke post restart","b1:",b1,"b0:",b0,"B1:",B1,"sum W:",sum(W)))

             # If still singular, then the Rocke estimate is returned
             # to avoid bootstrap problems. Perhaps better the previous sfast estimate?

             if (!is.finite(B1)){
                warning(paste("Not even Rocke start helps, passing over the start values","b1:",b1,"b0:",b0))
                B1 <- b1
                B0 <- b0
                break
             }

          } else {
          warning("Catch all break, should never happen")
          B1<-b1
          B0<-b0
          break
          }

      }

      #Calculation of B0 with a finite B1

      B0 <- YW-B1*XW

      if(abs(b1-B1) < threshold & abs(b0-B0) < threshold){
        #if(error.count>0){
        #  message(paste("final B1:",B1,"B0:",B0))
        #}
        break
      }

      # new values
      b1<-B1
      b0<-B0
    } # end of iterative algorithm

    return(list(b1=B1,b0=B0,iter=i,xw=XW, weight=W))
  }
#}

