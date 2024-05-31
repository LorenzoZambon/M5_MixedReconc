
compute_scores = function(STORE, h_list = 1:14, alpha_mis = 0.1, 
                          n_samp_gtrunc = 1e4, results_path, jitt = 1e-9) {
  
  A = get_hier_M5()$A
  n_u = nrow(A)
  n_b = ncol(A)
  n = n_u + n_b
  
  store_path = paste0(results_path, STORE)
  scores_str = paste0(store_path, "/scores.rds")
  if (file.exists(scores_str)) {
    return()
  }
  
  print("Computing scores...")
  
  # Initialize lists for the scores
  mase = list("base"        = matrix(nrow = n, ncol = length(h_list)),
              "gauss"       = matrix(nrow = n, ncol = length(h_list)),
              "gauss_trunc" = matrix(nrow = n, ncol = length(h_list)),
              "mixed_cond"  = matrix(nrow = n, ncol = length(h_list)),
              "TD_cond"     = matrix(nrow = n, ncol = length(h_list)))
  mis = mase
  rps = mase
  
  for (j in seq_along(h_list)) {
    
    h = h_list[[j]]
    print(paste0("h = ", h))
    h_path = paste0(store_path, "/h=", h)
    
    fc_bottom = readRDS(paste0(h_path, "/base_fc_bottom.rds"))
    fc_upper  = readRDS(paste0(h_path, "/base_fc_upper.rds"))
    rec_fc    = readRDS(paste0(h_path, "/rec_fc.rds"))
    
    actuals_u = unlist(lapply(fc_upper, "[[", "actual"))
    actuals_b = unlist(lapply(fc_bottom, "[[", "actual"))
    actuals = c(actuals_u, actuals_b)
    
    # compute scaling factors Q for MASE
    train_u = lapply(fc_upper, "[[", "train")
    train_b = lapply(fc_bottom, "[[", "train")
    Q_u = unlist(lapply(train_u, function(x) mean(abs(x[-1] - x[-length(x)]))))
    Q_b = unlist(lapply(train_b, function(x) mean(abs(x[-1] - x[-length(x)]))))
    Q = c(Q_u, Q_b)
    
    ### BASE ###
    # Upper
    mu_u = unlist(lapply(fc_upper, "[[", "mu"))
    sd_u = unlist(lapply(fc_upper, "[[", "sigma")) 
    mase$base[1:n_u,j] = abs(mu_u - actuals_u) / Q_u
    mis$base[1:n_u,j]  = mapply(MIS_gauss, mu_u, sd_u, actuals_u, 
                                MoreArgs = list(alpha=alpha_mis))
    rps$base[1:n_u,j] = scoringRules::crps(actuals_u, "norm", mean=mu_u, sd=sd_u)
    # Bottom
    base_fc_pmf = lapply(fc_bottom, "[[", "pmf")
    mase$base[(n_u+1):n,j] = mapply(PMF.AE, base_fc_pmf, actuals_b) / Q_b
    mis$base[(n_u+1):n,j] = mapply(PMF.IS, base_fc_pmf, actuals_b, 
                                   MoreArgs = list(alpha=alpha_mis))
    rps$base[(n_u+1):n,j] = mapply(PMF.RPS, base_fc_pmf, actuals_b)
    
    ### GAUSSIAN ###
    mu = c(rec_fc$gauss$mu_u, rec_fc$gauss$mu_b)
    sd = c(diag(rec_fc$gauss$Sigma_u)**0.5, diag(rec_fc$gauss$Sigma_b)**0.5)
    sd = sd + jitt  # for numerical stability
    mase$gauss[,j] = abs(mu - actuals) / Q
    mis$gauss[,j] = mapply(MIS_gauss, mu, sd, actuals, 
                           MoreArgs = list(alpha=alpha_mis)) 
    rps$gauss[,j] = scoringRules::crps(actuals, "norm", mean=mu, sd=sd)
    
    ### TRUNCATED GAUSSIAN ### 
    # Bottom
    mu_b = rec_fc$gauss$mu_b
    mu_b_trunc = mu_b
    mu_b_trunc[mu_b_trunc<0] = 0
    sd_b = diag(rec_fc$gauss$Sigma_b)**0.5
    mase$gauss_trunc[(n_u+1):n,j] = abs(mu_b_trunc - actuals_b) / Q_b
    mis$gauss_trunc[(n_u+1):n,j] = mapply(MIS_gauss, mu_b_trunc, sd_b, actuals_b, 
                                          MoreArgs = list(alpha=alpha_mis)) 
    Sigma_b = rec_fc$gauss$Sigma_b
    Sigma_b = Sigma_b + diag(jitt, nrow(Sigma_b))
    samples_b = t(.MVN_sample(n_samp_gtrunc, mu_b, Sigma_b))  # n x n_samp
    samples_b[samples_b<0] = 0
    rps$gauss_trunc[(n_u+1):n,j] = scoringRules::crps_sample(actuals_b, samples_b)
    # Upper
    samples_u = A %*% samples_b
    rm(samples_b)  # for memory
    medians_u = apply(samples_u, 1, median)
    mase$gauss_trunc[1:n_u,j] = abs(medians_u - actuals_u) / Q_u
    mis$gauss_trunc[1:n_u,j] = MIS_samples(samples_u, actuals_u, alpha_mis)
    rps$gauss_trunc[1:n_u,j] = scoringRules::crps_sample(actuals_u, samples_u)
    
    ### MIXED CONDITIONING ###
    mc_pmf = c(rec_fc$mixed_cond$upper, rec_fc$mixed_cond$bottom)
    mase$mixed_cond[,j] = mapply(PMF.AE, mc_pmf, actuals) / Q
    mis$mixed_cond[,j] = mapply(PMF.IS, mc_pmf, actuals, 
                                MoreArgs = list(alpha=alpha_mis))
    rps$mixed_cond[,j] = mapply(PMF.RPS, mc_pmf, actuals)
    
    ### TOP-DOWN CONDITIONING ###
    td_pmf = c(rec_fc$TD_cond$upper, rec_fc$TD_cond$bottom)
    mase$TD_cond[,j] = mapply(PMF.AE, td_pmf, actuals) / Q
    mis$TD_cond[,j] = mapply(PMF.IS, td_pmf, actuals,
                             MoreArgs = list(alpha=alpha_mis))
    rps$TD_cond[,j] = mapply(PMF.RPS, td_pmf, actuals)
    
    scores = list(
      mase = mase,
      mis  = mis,
      rps  = rps
    )
    saveRDS(scores, scores_str)
  }

}