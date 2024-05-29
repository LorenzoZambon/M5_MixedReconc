################################################################################
# Utils for computing base forecasts 

get_upp_ts = function(store.train, store.test, lev, 
                      name, h, len = 1941) {
  
  if (!(lev %in% c("store_id", "cat_id", "dept_id"))) stop("Wrong lev name")
  
  df1 = store.train[store.train[[lev]] == name]
  df2 = store.test[store.test[[lev]] == name]
  st = which(colnames(df1)=="d_1")
  
  serie = cbind(
    df1[,st:(st+len-1)],
    df2[,st:(st+h-1)] )
  data.train = serie[,1:(len+h-1)]
  data.test = serie[,(len+h):(len+h)]
  train.agg = colSums(data.train)
  test.agg = colSums(data.test)
  return(list(train = train.agg, test = test.agg))
}

get_bott_ts = function(store.train, store.test, 
                       item.id, h, len = 1941) {
  
  df1 = store.train[store.train$item_id == item.id]
  df2 = store.test[store.test$item_id == item.id]
  st = which(colnames(df1)=="d_1")
  
  serie = c(as.numeric(df1[1,st:(st+len-1)]), as.numeric(df2[1,st:(st+h-1)]))
  train = serie[1:(len+h-1)]
  test = serie[[len+h]]
  
  return(list(train = train, test = test))
}

# Fit model on upper TS
model_upper = function(data, model_str = "AXX", 
                       occ_str = "none", distr = "dnorm") {
  model = auto.adam(data, model_str, lags = c(7), 
                    occurrence = occ_str, distribution = distr)
  return(model)
}

# Fit model on bottom TS
model_bottom = function(data, model_str = "MNN", 
                        occ_str = "auto", distr = c("dgamma")) {
  model = adam(data, model_str, lags = c(7), 
               occurrence = occ_str, distribution = distr)
  return(model)
}

# Check that the samples are discrete 
.check_discrete_samples <- function(samples) {
  if (!isTRUE(all.equal(unname(samples), as.integer(samples)))) {
    stop("Input error: samples are not all discrete")
  }
}

# Compute pmf from samples
PMF.from_samples = function(v) {
  .check_discrete_samples(v)
  pmf = tabulate(v+1) / length(v)  # the support starts from 0 
  # Tabulate only counts values above 1: if sum(tabulate(v+1)) > length(v),
  # it means that there were negative samples
  if (!isTRUE(all.equal(sum(pmf), 1))) {
    stop("Input error: same samples are negative")
  }
  return(pmf)
}

################################################################################
# Return aggregation and summing matrices (A and S) for the M5-store hierarchy
get_hier_M5 = function() {
  
  n_b = 3049
  n_u = 11
  n = n_u + n_b
  
  A = matrix(0, nrow = n_u, ncol = n_b)
  A[1,] = 1           # store
  A[2,1:565] = 1      # categories
  A[3,566:1612] = 1
  A[4,1613:3049] = 1
  A[5,1:416] = 1      # departments
  A[6,417:565] = 1
  A[7,566:1097] = 1
  A[8,1098:1612] = 1
  A[9,1613:1828] = 1
  A[10,1829:2226] = 1
  A[11,2227:3049] = 1
  
  S = rbind(A, diag(rep(1, n_b)))
  
  return(list(A = A, S = S))
}

################################################################################
# Utils for computing scores

# Computes MIS for Gaussian distribution
# If trunc=TRUE, set all the negative mass to zero
MIS_gauss = function(mu, sd, actual, alpha=0.1, trunc=FALSE) {
  z = qnorm(1-(alpha/2))
  u = mu + z * sds
  l = mu - z * sds
  if (trunc) {  
    l = max(l, 0)  # if negative, set to zero
    u = max(u, 0)
  }  
  return(u - l + (2/alpha)*(l - actual)*(actual < l) + 
           (2/alpha)*(actual - u)*(actual > u) )
}

# Compute the absolute error for a pmf 
PMF.AE = function(pmf, actual) {
  median_pmf = PMF.get_quantile(pmf, 0.5)
  return(abs(median_pmf - actual))
}

# Compute the interval score (specified by alpha) for a pmf 
PMF.IS = function(pmf, actual, alpha = 0.1) {
  u = PMF.get_quantile(pmf, 1-(alpha/2))
  l = PMF.get_quantile(pmf, alpha/2)
  return(u - l + (2/alpha)*(l - actual)*(actual < l) + 
                 (2/alpha)*(actual - u)*(actual > u) )
}

# Compute the RPS for a pmf 
PMF.RPS = function(pmf, actual) {
  cdf <- cumsum(pmf) / sum(pmf)
  M <- length(cdf)
  # if actual is outside the supp of pmf, add ones to the end of the cdf:
  if (actual >= M) {  
    cdf <- c(cdf, rep(1, (actual-M+1)))
    M <- length(cdf)
  }
  cdf_act <- (0:(M-1)) >= actual   # unit step function in actual
  crps_ <- sum((cdf - cdf_act)**2)
  return(crps_)
}

# Compute MIS (specified by alpha) for samples (dim: len(actual) x n_sample)
MIS_samples = function(samples, actual, alpha=0.1) {
  if (nrow(samples)!=length(actual)) stop(paste0(samples, " must be of shape length(", actual, ") x n_samples"))
  u = apply(samples, 1, function(x) quantile(x, 1-(alpha/2)))
  l = apply(samples, 1, function(x) quantile(x, alpha/2))
  return(u - l + (2/alpha)*(l - actual)*(actual < l) + 
                 (2/alpha)*(actual - u)*(actual > u) )
}

################################################################################
# Sample from a multivariate Gaussian distribution with specified mean and cov. matrix
.MVN_sample = function(n_samples, mu, Sigma) {
  n = length(mu)
  if (any(dim(Sigma) != c(n,n))) {
    stop("Dimension of mu and Sigma are not compatible!")
  } 
  .check_cov(Sigma, "Sigma", pd_check = FALSE, symm_check = FALSE)
  
  Z = matrix(stats::rnorm(n*n_samples), ncol = n)
  
  Ch <- tryCatch(base::chol(Sigma),
                 error = function(e) stop(paste0(e,"check the covariance in .MVN_sample, the Cholesky fails.")))
  
  samples = Z %*% Ch + matrix(mu, nrow = n_samples, ncol = n, byrow = TRUE)
  return(samples)
}

################################################################################
# Function for computing the skill score
skill.score <- function(ref, met) {
  s <- (2 * (ref - met) / (ref + met)) * 100
  s[is.na(s)] <- 0  # if both numerator and denominator are 0, set skill score to 0
  return(s)
}

# Computing the skill scores for all methods, scores, and stores.
# Save the mean skill scores and return the entire list of skill scores
compute_SS = function(results_path, STORES, h_list, ref_met = "base", 
                      methods_ = c("gauss", "gauss_trunc", "mixed_cond", "TD_cond"), 
                      metrics = c("mase", "mis", "rps")) {
  
  # Load scores
  scores = list()
  for (STORE in STORES) {
    scores[[STORE]] = readRDS(paste0(results_path, STORE, "/scores.rds"))
  }
  
  # Compute skill scores
  SS_u = list()
  SS_b = list()
  for (m in metrics) {
    SS_u[[m]] = list()
    SS_b[[m]] = list()
    for (met in methods_) {
      ss_b = c()
      ss_u = c()
      for (STORE in STORES) {
        ss = skill.score(scores[[STORE]][[m]][[ref_met]],   # matrix 3060 x 14
                         scores[[STORE]][[m]][[met]])
        ss_u = c(ss_u, ss[1:11,])
        ss_b = c(ss_b, ss[12:3060,])
      }
      SS_u[[m]][[met]] = ss_u
      SS_b[[m]][[met]] = ss_b
    }
  }
  
  # Compute mean skill scores
  mean_SS_u = list()
  mean_SS_b = list()
  for (m in metrics) {
    mean_SS_u[[m]] = list()
    mean_SS_b[[m]] = list()
    for (met in methods_) {
      mean_SS_u[[m]][[met]] = mean(SS_u[[m]][[met]])
      mean_SS_b[[m]][[met]] = mean(SS_b[[m]][[met]])
    }
  }
  
  # Save mean skill scores
  mean_SS = list(
    upper  = mean_SS_u,
    bottom = mean_SS_b
  )
  saveRDS(mean_SS, paste0(results_path, "mean_skill_scores.rds"))
  
  # Return all skill scores
  SS = list(
    upper  = SS_u,
    bottom = SS_b
  )
  return(SS)
}
















