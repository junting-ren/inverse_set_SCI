SCB_linear_outcome = function(df_fit, model, grid_df, n_boot = 1000, alpha = 0.05, grid_df_boot = NULL){
  formula_ = as.formula(model)
  fit = lm(model, df_fit)
  y_hat = predict(fit, grid_df, se.fit = TRUE, level = 1 - alpha) # for constructing the whole 
  res_max_v = rep(0,n_boot)
  if(is.null(grid_df_boot)){
    grid_df_boot = grid_df
    y_hat_level_grid = y_hat
  }else{
    y_hat_level_grid = predict(fit, grid_df_boot, se.fit = TRUE, level = 1 - alpha)# bootstrap true mean
  }
  for(i in 1:n_boot){
    df_boot = df_fit[sample(1:dim(df_fit)[1], replace = T),]
    fit_boot = lm(model, df_boot)
    y_hat_boot = predict(fit_boot, grid_df_boot, level = 1 - alpha, se.fit = TRUE	) 
    residual = abs(y_hat_boot$fit - y_hat_level_grid$fit)/y_hat_boot$se.fit
    res_max_v[i] = max(residual)
  }
  thres = quantile(res_max_v, probs = 1 - alpha)
  sim_CB = data.frame(LowerBound = y_hat$fit - thres*y_hat$se.fit, Mean = y_hat$fit, UpperBound = y_hat$fit + thres*y_hat$se.fit, grid_df)
}

expit = function(x){
  1/(1+exp(-x))
}


SCB_logistic_outcome = function(df_fit, model, grid_df, n_boot = 1000, alpha = 0.05){
  formula_ = as.formula(model)
  fit =  suppressWarnings(glm(model, family = binomial(), data = df_fit)) # Suppress warning forfitted probabilities numerically 0 or 1
  y_hat = predict(fit, grid_df, se.fit = TRUE, level = 1 - alpha) # bootstrap true mean
  res_max_v = rep(0,n_boot)
  for(i in 1:n_boot){
    df_boot = df_fit[sample(1:dim(df_fit)[1], replace = T),]
    fit_boot = suppressWarnings(glm(model, family = binomial(), data = df_boot)) 
    y_hat_boot = predict(fit_boot, grid_df, level = 1 - alpha, se.fit = TRUE	)
    residual = abs(y_hat_boot$fit - y_hat$fit)/y_hat_boot$se.fit
    res_max_v[i] = max(residual)
  }
  thres = quantile(res_max_v, probs = 1 - alpha)
  sim_CB = data.frame(LowerBound = expit(y_hat$fit - thres*y_hat$se.fit), Mean = expit(y_hat$fit),
                      UpperBound = expit(y_hat$fit + thres*y_hat$se.fit), grid_df)
}

SCB_regression_coef = function(df_fit, model, n_boot = 5000, alpha = 0.05, type = "linear"){
  # type: "linear" or "logistic"
  formula_ = as.formula(model)
  if(type == "linear"){
    fit = lm(model, df_fit)
  }else if(type == "logistic"){
    fit = suppressWarnings(glm(model, family = binomial(), data = df_fit)) # Suppress warning for 0 or 1 probability for simulation
  }
  sum_out = summary(fit)
  coef_hat = fit$coefficients # bootstrap true mean
  coef_sd = sum_out$coefficients[,2]
  res_max_v = rep(0,n_boot)
  for(i in 1:n_boot){
    df_boot = df_fit[sample(1:dim(df_fit)[1], replace = T),]
    if(type == "linear"){
      fit_boot = lm(model, df_boot)
    }else if(type == "logistic"){
      fit_boot = suppressWarnings(glm(model, family = binomial(), data = df_boot))
    }
    sum_out_boot = summary(fit_boot)
    coef_hat_boot =fit_boot$coefficients
    coef_sd_boot = sum_out_boot$coefficients[,2]
    residual = abs(coef_hat_boot - coef_hat)/coef_sd_boot
    res_max_v[i] = max(residual)
  }
  thres = quantile(res_max_v, probs = 1 - alpha)
  sim_CB = data.frame(LowerBound = coef_hat - thres*coef_sd, Mean = coef_hat, UpperBound = coef_hat + thres*coef_sd)
}