
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