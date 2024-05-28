base_fc_store = function(STORE, h_list, n_samples_bott, 
                         data_path, results_path, n_cpu) {
  
  # Create folder to save results
  store_path = paste0(results_path, STORE)
  if (!dir.exists(store_path)) dir.create(store_path)
  
  # Load train and test
  store.train = m5::m5_get_raw_evaluation(data_path)[[1]]
  store.train = store.train[store.train$store_id == STORE]
  store.test  = m5::m5_get_raw_evaluation(data_path)[[2]]
  store.test  = store.test[store.test$store_id == STORE]
  len = 1941  # length of the training set
  
  ### Base forecasts upper time series
  print("Computing upper base forecasts...")
  for (h in h_list) {
    
    # Create folder for h
    h_path = paste0(store_path, "/h=", h)
    if (!dir.exists(h_path)) {
      dir.create(h_path)
    }
    
    # If upper base fc have already been computed, skip
    str_upp = paste0(h_path, "/base_fc_upper.rds")
    if (file.exists(str_upp)) {
      print(paste0("Upper base forecasts for h = ", h, " already computed"))
      next
    } 
    
    fc = list()
    
    # Store
    ts.agg = get_upp_ts(store.train, store.test, "store_id", 
                        STORE, h = h, len = len)
    train.agg = ts.agg$train
    test.agg = ts.agg$test
    
    model = model_upper(train.agg)       
    fc.model = forecast(model, h = 1) 
    fc[[STORE]] = list(
      mu = fc.model$mean,
      sigma = model$scale,
      residuals = model$residuals,
      actual = test.agg,
      train = train.agg
    )
    
    # Category
    for (cat.id in CAT) {
      ts.agg = get_upp_ts(store.train, store.test, "cat_id", 
                          cat.id, h = h, len = len)
      train.agg = ts.agg$train
      test.agg = ts.agg$test
      
      model = model_upper(train.agg)       
      fc.model = forecast(model, h = 1) 
      fc[[cat.id]] = list(
        mu = fc.model$mean,
        sigma = model$scale,
        residuals = model$residuals,
        actual = test.agg,
        train = train.agg
      )
    }
    
    # Department
    for (dept.id in DEPT) {
      ts.agg = get_upp_ts(store.train, store.test, "dept_id", 
                          dept.id, h = h, len = len)
      train.agg = ts.agg$train
      test.agg = ts.agg$test
      
      model = model_upper(train.agg)       
      fc.model = forecast(model, h = 1) 
      fc[[dept.id]] = list(
        mu = fc.model$mean,
        sigma = model$scale,
        residuals = model$residuals,
        actual = test.agg,
        train = train.agg
      )
    }
    
    # Save results
    saveRDS(fc, str_upp)
  }
  
  ### Base forecasts bottom time series
  n_cores = min(n_cpu, length(h_list))  # for parallel computation
  print(paste0("Computing bottom base forecasts using ", n_cores, " cores..."))

  cl <- makeCluster(n_cores) 
  registerDoSNOW(cl)
  foreach(h = h_list, .packages = c("data.table", "smooth"),
          .export=c("get_bott_ts", "model_bottom", "PMF.from_samples",
                    ".check_discrete_samples")) %dopar% {
            
            # If bottom base fc have already been computed, skip
            h_path = paste0(store_path, "/h=", h)
            str_bott = paste0(h_path, "/base_fc_bottom.rds")
            if (file.exists(str_bott)) return(NULL)
            
            fc = list()

            for (item.id in unique(store.train$item_id)) {
              
              bts = get_bott_ts(store.train, store.test, item.id, h)
              train = bts$train
              test = bts$test

              model = model_bottom(train)
              fc.model = forecast(model, h = 1, interval="simulated",
                                  scenarios=TRUE, nsim = n_samples_bott)

              # round to integer
              samples = round(fc.model$scenarios[1,])
              samples[samples<0] = 0    # set negative to zero
              # get pmf from samples
              pmf = PMF.from_samples(samples)

              fc[[item.id]] = list(
                pmf = pmf,
                actual = test,
                train = train
              )
            }
            saveRDS(fc, str_bott)
          }
  stopCluster(cl)
}

