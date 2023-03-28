
SCB_dense = function(A, alpha_level = 0.05){
  require(SIRF)
  # Get the number of repeats and dimension
  dimA = dim(A)
  N = dimA[length(dimA)]
  D = length(dimA) - 1
  # getting the estimated sd
  sd_A = sqrt(apply(A, 1:D, var))
  # getting the estimated mean
  mean_A <- apply(A, 1:D, mean)
  mean_A_array = array(rep(mean_A, N), dim = c(dimA[1:D], N))
  # Get the threshold
  thres = SIRF::MultiplierBootstrap(sqrt(N/(N - 1))*(A-mean_A_array), alpha = alpha_level)$q # the default is the multiplier-t with Rademacher multipliers, here alpha is a 
  # construct the SCB
  scb_up = mean_A + thres*sd_A/sqrt(N)
  scb_low = mean_A - thres*sd_A/sqrt(N)
  # return index, scb_up, scb_low
  return(list(scb_up = scb_up, scb_low = scb_low, thres = thres))
}