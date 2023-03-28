scb_to_cs = function(scb_up, scb_low, levels, true_mean = NULL,est_mean = NULL, x1, x2, type = "upper", return_contain_only = F, return_plot = F,...)
  {
  # scb_up: the upper simultaneous confidence interval, a matrix with the same dimension with index
  # scb_low: the lower simultaneous confidence interval, a matrix with the same dimension with index
  # levels: a list f scalers for different levels or matrix containing interval sets to construct the confidence sets
  # true_mean: a matrix with the same dimension as index, with elements of the true mean, only for simulation.
  # est_mean: only for plotting: if you do not have the true mean
  # x1: the coordination for the first dimension
  # x2: the coordination for the second dimension
  # type: the type of inverse set we are building if the levels are not a matrix, values: upper, lower, two-sided
  # return_contain_only: whether to return the contain T or F only
  # The type I error rate is determined by SCB type I error rate
  incl_f <- function(A, B) {
    min(B - A) >= 0 #if B includes A
  }
  if(return_plot){
    p_para <- list(...)
    if(is.null(true_mean)){
      pl_together = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low), 
                            levels = levels, x = x1, y = x2, mu_hat = est_mean,
                            together = T, xlab = p_para$xlab,ylab = p_para$ylab)
      pl = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low), 
                   levels = levels, x = x1, y = x2, mu_hat = est_mean, 
                   together = F, xlab = p_para$xlab,ylab = p_para$ylab)
      pl = list(pl_together, pl)
    }else{
      pl_together = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low), 
                            levels = levels, x = x1, y = x2,mu_true = true_mean, 
                            mu_hat = est_mean, together = T,xlab = p_para$xlab,ylab = p_para$ylab)
      pl = plot_cs(SCB = list(scb_up = scb_up, scb_low = scb_low), levels = levels, 
                   x = x1, y = x2,mu_true = true_mean, mu_hat = est_mean, 
                   together = F,
                   xlab = p_para$xlab,ylab = p_para$ylab)
      pl = list(pl_together, pl)
    }
  }else{
    pl = NULL
  }
  contain_v = c()
  in_list = list()
  out_list = list()
  j = 1
  if(is.vector(levels)){
    if(type == "upper"){
      for(c in levels){
        in_set = scb_low >= c
        out_set = scb_up >= c
        if(!return_contain_only){
          in_list[[j]] = in_set
          out_list[[j]]= out_set
        }
        if(!is.null(true_mean)){
          true_set = true_mean >= c
          contain_v = c(contain_v,incl_f(in_set, true_set) & incl_f(true_set, out_set))
        }
        j = j + 1
      }
    }else if(type == "lower"){
      for(c in levels){
        in_set = scb_up <= c
        out_set = scb_low <= c
        if(!return_contain_only){
          in_list[[j]] = in_set
          out_list[[j]]= out_set
        }
        if(!is.null(true_mean)){
          true_set = true_mean <= c
          contain_v = c(contain_v,incl_f(in_set, true_set) & incl_f(true_set, out_set))
        }
        j = j + 1
      }
    }else if(type == "two-sided"){
      out_set_low_list = list()
      out_set_up_list = list()
      for(c in levels){
        out_set_low = scb_low <= c
        out_set_up = scb_up >= c
        if(!return_contain_only){
          out_set_low_list[[j]] = out_set_low
          out_set_up_list[[j]]= out_set_up
        }
        if(!is.null(true_mean)){
          true_set_low = true_mean <= c
          true_set_up = true_mean >= c
          contain_v = c(contain_v,incl_f(true_set_low, out_set_low) & incl_f(true_set_up, out_set_up))
        }
        j = j + 1
      }
      return(list(levels = levels, L_out = out_set_low_list, U_out = out_set_up_list, 
                  contain_individual = contain_v, contain_all = all(contain_v), plot_cs = plot_cs))
    }

  }else{# When we have interval level sets
    l_dim = dim(levels)
    for(i in 1:l_dim[1]){
      c = c(levels[i,])
      in_set = scb_low >= c$low & scb_up <= c$up
      out_set = scb_up >= c$low & scb_low <= c$up
      if(!return_contain_only){
        in_list[[i]] = in_set
        out_list[[i]]= out_set
      }
      if(!is.null(true_mean)){
        true_set = true_mean >= c$low & true_mean <= c$up
        contain_v = c(contain_v,incl_f(in_set, true_set) & incl_f(true_set, out_set))
      }
    }
  }
  return(list(levels = levels, U_in = in_list, U_out = out_list, 
         contain_individual = contain_v, contain_all = all(contain_v), plot_cs = pl))
}