
SCB_gls_climate = 
function (Z, level, X = NULL, w = NULL, correlation = NULL, corpar = NULL, 
          groups = NULL, V = NULL, alpha = 0.1, N = 1000, mu = NULL, 
          mask = NULL) 
{
  require(nlme)
  x = Z$x
  y = Z$y
  Y = Z$z
  n = dim(Y)[3]
  nloc <- length(x) * length(y)
  if (is.null(X)) {
    X <- matrix(1, n, 1)
    w <- matrix(1, 1, 1)
  }
  p <- ncol(X)
  if (is.null(groups)) {
    groups <- rep(1, n)
  }
  invsqrtm <- function(A) {
    E <- eigen(A)
    U <- E$vectors
    D <- diag(E$values)
    U %*% diag(1/sqrt(E$values)) %*% t(U)
  }
  deR <- array(0, c(length(x), length(y), n))
  vabs <- matrix(0, length(x), length(y))
  norm_est <- matrix(0, length(x), length(y))
  mu_hat <- matrix(0, length(x), length(y))
  if (!is.null(correlation)) 
    correlation = do.call(get(correlation), c(corpar, form = ~1 | 
                                                groups))
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      ytemp <- Y[i, j, ]
      if (sum(is.na(ytemp)) == length(ytemp)) {
        mu_hat[i, j] = NA
        norm_est[i, j] = NA
      }
      else {
        df <- data.frame(cbind(ytemp = ytemp, X, groups = groups))
        df <- df[order(groups), ]
        groups <- sort(groups)
        fo <- paste(names(df)[1], "~", paste(names(df)[-c(1, 
                                                          p + 2)], collapse = " + "), "-1")
        if (is.null(V)) {
          model <- nlme::gls(formula(fo), data = df, 
                             correlation = correlation)
        }
        else {
          model <- MASS::lm.gls(formula(fo), data = df, 
                                W = V[i, j, , ], inverse = TRUE)
        }
        mu_hat[i, j] <- t(w) %*% model$coefficients
        if (!is.null(correlation)) {
          cM <- nlme::corMatrix(model$modelStruct$corStruct, 
                                corr = F)
          if (!is.list(cM)) 
            cM <- list(cM)
          invsqrtmOmega <- Matrix::as.matrix(Matrix::bdiag(cM))
          deR[i, j, ] <- invsqrtmOmega %*% model$residuals
          deR[i, j, ] <- deR[i, j, ]/sd(deR[i, j, ])
        }
        else if (!is.null(V)) {
          deR[i, j, ] <- chol(solve(V[i, j, , ])) %*% 
            model$residuals
          deR[i, j, ] <- deR[i, j, ]/sd(deR[i, j, ])
          model$varBeta = solve(t(X) %*% solve(V[i, j, 
                                                 , ], X))
        }
        else {
          deR[i, j, ] <- model$residuals
          deR[i, j, ] <- deR[i, j, ]/sd(deR[i, j, ])
        }
        vabs[i, j] <- sqrt(t(w) %*% model$varBeta %*% 
                             w)
        norm_est[i, j] <- (mu_hat[i, j] - level)/vabs[i, 
                                                      j]
      }
    }
  }
  if (is.null(mask)) {
    mask = array(1, dim = c(length(x), length(y)))
  }
  else {
    mask[which(!is.na(mask), arr.ind = TRUE)] = 1
  }
  mu_hat <- mu_hat * mask
  norm_est <- norm_est * mask
  deR_mask = sweep(deR, 1:2, mask, FUN = "*" )
  a_MB = quantile(MB_(x = x, y = y, R = deR, N = N), 
                  probs = 1 - alpha, type = 8)
  norm_est[i, j] <- (mu_hat[i, j] - level)/vabs[i, 
                                                j]
  scb_up = mu_hat + a_MB*vabs
  scb_low = mu_hat - a_MB*vabs
  # return index, scb_up, scb_low
  return(list(scb_up = scb_up, scb_low = scb_low, mu_hat = mu_hat, thres = a_MB,
              x = x, y = y))
}

MB_ = function (x, y, R, N = 1000) 
{
  n = dim(R)[3]
  g = matrix(rnorm(n * N), n, N)
  apply(abs(matrix(R, ncol = n) %*% g), 2, max, na.rm = T)/sqrt(n - 
                                                              2)
}




###################################
#Plotting
###################################
library(fields)
library(maps)
data(worldMapEnv)

#Functions
##################################
#Plotting values on the map.
image.map = function(lon, lat, img, mask=NULL, xlab='longitude', ylab='latitude', ...) {
  if(!is.null(mask)) {
    img = img*mask
    xlim = lon[range(which(rowSums(mask, na.rm=TRUE)>0))]
    ylim = lat[range(which(colSums(mask, na.rm=TRUE)>0))]
  } else {
    xlim = range(lon)
    ylim = range(lat)
  }
  image.plot(lon, lat, img, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  map('world', add=TRUE)
}

#Drawing the contour.
drawContour<-function(x,y,z,c,col,lty=1){
  C<-contourLines(x,y,z,levels=c,nlevels=1)
  for(i in 1:length(C))
    lines(C[[i]]$x,C[[i]]$y,pch=20,col=col,lwd=3,lty=lty)
}

#Create single peak test signal.
single.peak = function(sx, sy, x, y, b) {
  mu0 = matrix(0, length(sx), length(sy))
  mu0[x, y] = 1
  mu = image.smooth(mu0, theta = b, dx = sx[2]-sx[1], dy = sy[2]-sy[1])$z
}

#A and B are coincident matrices representing subsets of the plane (non-zero entry = point contained in set).
#This function returns 1 if A is a subset of B and zero otherwise.
subset = function(A,B){
  all(!(A & (!B)))
}

false_pos = function(A_true,A_plus){
  sum(A_plus & !A_true) / sum(!A_true)
}

false_neg = function(A_true,A_minus){
  sum(A_true & !A_minus) / sum(A_true)
}

#This function takes realizations Y of a random field and returns a distribution functions
# for the supremum of the limiting Gaussian field.
MC_gauss = function(Y,N){
  n = dim(Y)[3]
  Y = matrix(Y,ncol=n)
  g = matrix(rnorm(n*N),n,N)
  maxima = apply(abs(Y %*% g), 2, max) / sqrt(n)
  
  function(t) sum(maxima>=t) / length(maxima)
}

#Does the same as MC_gauss but uses a plug-in estimate of the excursion set A_c.
MC_plug = function(Y,N,A){
  n = dim(Y)[3]
  Y = matrix(Y,ncol=n)
  A = as.vector(A)
  A[A==TRUE] = -1
  A[A==FALSE] = 1
  g = matrix(rnorm(n*N),n,N)
  r = apply(Y %*% g * A,2,max) / sqrt(n)
  function(t) sum(r>=t) / length(r)
}

#This function computes an estimate for the variance upper bound of the false area ratios by Monte Carlo integration.
variance_mc = function(x,y,Y,A,a,N){
  n = dim(Y)[3]
  
  rx = range(x)
  ry = range(y)
  
  #Draw uniformly random points in SxS.
  ux1 = runif(N,min=rx[1],max=rx[2])
  uy1 = runif(N,min=ry[1],max=ry[2])
  ux2 = runif(N,min=rx[1],max=rx[2])
  uy2 = runif(N,min=ry[1],max=ry[2])
  
  #Interpolate values of realizations to the locations of points.
  Y1 = matrix(0,n,N)
  for(i in 1:n) Y1[i,] = interp.surface(list(x=x,y=y,z=Y[,,i]),cbind(ux1,uy1))
  Y2 = matrix(0,n,N)
  for(i in 1:n) Y2[i,] = interp.surface(list(x=x,y=y,z=Y[,,i]),cbind(ux2,uy2))
  
  #Compute empirical covariance.
  cv = cov(Y1,Y2)
  
  #Determine if random draws are in A or !A.
  IndA1 = interp.surface(list(x=x,y=y,z=A),cbind(ux1,uy1))
  IndA2 = interp.surface(list(x=x,y=y,z=A),cbind(ux2,uy2))
  
  IndAc1 = interp.surface(list(x=x,y=y,z=!A),cbind(ux1,uy1))
  IndAc2 = interp.surface(list(x=x,y=y,z=!A),cbind(ux2,uy2))
  
  #Function from the variance bound as in the text.
  M = function(cv) cv/(2*pi*sqrt(1-cv^2)) * exp(-1/(1+cv)*a^2)
  
  #Computing volumes of A and !A.
  volA = sum(A) * diff(range(x))/length(x) * diff(range(y))/length(y)
  volAc = sum(!A) * diff(range(x))/length(x) * diff(range(y))/length(y)
  
  #Compute MC integral.
  c(mean(IndA1 %o% IndA2 * M(cv)) * length(A)^2 / sum(A)^2 , mean(IndAc1 %o% IndAc2 * M(cv)) * length(A)^2 / sum(!A)^2)
}