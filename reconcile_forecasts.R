rec_fc_store <- function(STORE, h_list, results_path,
                         N_samples_IS, N_samples_TD) {
  
  store_path = paste0(results_path, STORE)
  
  A = get_hier_M5()$A
  S = get_hier_M5()$S
  
  n_u = nrow(A)
  n_b = ncol(A)
  n = n_b + n_u
  
  timings <- matrix(nrow=length(h_list), ncol=3)
  
  for (h in h_list) {
    print(paste0("Reconciliation of store ", STORE, " for h = ", h))
    
    h_path = paste0(store_path, "/h=", h)
    
    str_bott = paste0(h_path, "/base_fc_bottom.rds")
    str_upp  = paste0(h_path, "/base_fc_upper.rds")
    fc_bottom = readRDS(str_bott)
    fc_upper = readRDS(str_upp)
    
    rec_str = paste0(h_path, "/rec_fc.rds") 
    
    # If rec file already present, skip
    if (file.exists(rec_str)) {
      next
    } 
    
    # Initialize list for the results
    rec_fc = list(
      gauss      = list(),
      mixed_cond = list(),
      TD_cond    = list()
    )
    
    # Bottom base forecasts (in the form of PMFs)
    bottom_pmf = lapply(fc_bottom, "[[", "pmf")
    # Upper base forecasts (multivariate Gaussian)
    mu_u = unlist(lapply(fc_upper, "[[", "mu"))
    residuals.upper = lapply(fc_upper, "[[", "residuals")
    residuals.upper = t(do.call("rbind", residuals.upper))
    # covariance matrix obtained using shrinkage estimator:
    Sigma_u = bayesRecon::schaferStrimmer_cov(residuals.upper)$shrink_cov
    upper_params = list(mu=mu_u, Sigma=Sigma_u)
    
    ### Gaussian ###
    # Fit Gaussians on the bottom PMFs
    mu_b = unlist(lapply(bottom_pmf, bayesRecon::PMF.get_mean))
    var_b = unlist(lapply(bottom_pmf, bayesRecon::PMF.get_var))
    Sigma_b = diag(var_b)  # bottom are independent
    
    mu    <- c(mu_u, mu_b)
    Sigma <- matrix(0, nrow = n, ncol = n)
    Sigma[1:n_u, 1:n_u] <- Sigma_u
    Sigma[(1+n_u):n, (1+n_u):n] <- Sigma_b
    
    tic("g")
    gauss = reconc_gaussian(S, mu, Sigma)
    toc(log = TRUE, quiet = TRUE)
    
    rec_fc$gauss = list(mu_b    = gauss$bottom_reconciled_mean,
                        Sigma_b = gauss$bottom_reconciled_covariance,
                        mu_u    = A %*% gauss$bottom_reconciled_mean,
                        Sigma_u = A %*% gauss$bottom_reconciled_covariance %*% t(A))

    ### Mixed conditioning ###
    tic("mc")
    mixcond = reconc_MixCond(S, bottom_pmf, upper_params, bottom_in_type = "pmf",
                             num_samples = N_samples_IS, return_type = "pmf")
    toc(log = TRUE, quiet = TRUE)
    
    rec_fc$mixed_cond = list(
      bottom = mixcond$bottom_reconciled$pmf,
      upper = mixcond$upper_reconciled$pmf,
      ESS = mixcond$ESS
    )
      
    ### Top-down conditioning ###  
    tic("td")
    td = reconc_TDcond(S, bottom_pmf, upper_params, bottom_in_type = "pmf", 
                       num_samples = N_samples_TD, return_type = "pmf")
    toc(log = TRUE, quiet = TRUE)
    
    rec_fc$TD_cond = list(
      bottom = td$bottom_reconciled$pmf,
      upper = td$upper_reconciled$pmf
    )
    
    # Save results
    saveRDS(rec_fc, rec_str)

    # Save computational times
    log.lst <- tic.log(format=FALSE)
    timings[h,] <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
    tic.clearlog()
    
  }
  saveRDS(timings, paste0(store_path,"/rec_times.rds"))
  
}