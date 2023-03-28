BernsteinSumNoise = function (N, x = seq(0, 1, length.out = 100), sigma = function(x) {
  1
}, randNumber = rnorm) 
{
  f <- cbind((1 - x)^6, 6 * x * (1 - x)^5, 15 * x^2 * (1 - 
                                                         x)^4, 20 * x^3 * (1 - x)^3, 15 * x^4 * (1 - x)^2, 6 * 
               x^5 * (1 - x), x^6)
  fSqSum <- apply(f^2, 1, sum)
  fNorm <- f/sqrt(fSqSum)
  nBasis <- dim(f)[2]
  return(fNorm %*% matrix(randNumber(nBasis * N), nBasis, N) * 
           sigma(x))
}

GaussDensitySum2DNoise = function (N, x = seq(0, 1, length.out = 50), sigma = function(x) array(outer(x, 
                                                                                                      x, FUN = function(s, t) (s + 1)/(t^2 + 1))/3, c(rep(length(x), 
                                                                                                                                                          2), 1)), randNumber = rnorm, M = 6, bdx = 0.02, Fnorm = NULL) 
{
  if (is.null(Fnorm)) {
    nT <- nT1 <- nT2 <- length(x)
    rx = range(x)
    grd = seq(rx[1] + bdx, rx[2] - bdx, length.out = M)
    ptGrid = as.matrix(expand.grid(grd, grd))
    FF = NULL
    for (k in 1:dim(ptGrid)[1]) {
      F1 = outer(x, x, FUN = Vectorize(function(s, t) {
        dnorm(sqrt(sum((c(s, t) - ptGrid[k, ])^2)), mean = 0, 
              sd = 0.1)
      }, vectorize.args = c("s", "t")))
      FF = abind::abind(FF, F1, along = 3)
    }
    FFsqSum = sqrt(apply(FF^2, MARGIN = c(1, 2), sum))
    Fnorm = matrix(FF/array(rep(FFsqSum, M^2), dim = c(nT, 
                                                       nT, M^2)), nT^2, M^2)
    rm(FF, FFsqSum)
  }
  else {
    dimF = dim(Fnorm)
    M = sqrt(dimF[3])
    nT1 = dimF[1]
    nT2 = dimF[2]
    Fnorm = matrix(Fnorm, nT1 * nT2, M^2)
  }
  sd_eval = array(rep(sigma(x), N), c(rep(length(x), 2), N))/3
  rcoef <- matrix(randNumber(M^2 * N), M^2, N)
  return(sd_eval * array(Fnorm %*% rcoef, dim = c(nT1, nT2, 
                                                  N)))
}


FunctionalDataSample = function (N, x = seq(0, 1, length.out = 100), mu = function(x) {
  rep(0, length(x))
}, noise = SinCosSumNoise, sigma = function(x) {
  rep(1, length(x))
}, sd_ObsNoise = 0, ...) 
{
  m = mu(x)
  if (is.vector(m)) {
    mdim = length(x)
    fac = 1
    m = matrix(m, mdim, N)
  }
  else {
    mdim = dim(m)
    if (is.vector(x)) {
      fac = length(mdim)
    }
    else {
      fac = 1
    }
    m = array(rep(m, N), c(mdim, N))
  }
  return(m + noise(N = N, x = x, sigma = sigma, ...) + array(rnorm(N * 
                                                                     length(x) * fac, mean = 0, sd = sd_ObsNoise), dim = c(mdim, 
                                                                                                                           N)))
}