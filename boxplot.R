# Computes the skill scores and produces the boxplots, which are saved in BP_path
produce_boxplots = function(STORES, h_list, BP_met, results_path) {
  
  BP_path = paste0(results_path, "boxplots/")
  
  # Compute skill scores
  SS = compute_SS(results_path, STORES, h_list)
  
  ### Boxplots ###
  custom_colors <- c("#a8a8e4", 
                     "#a9c7e4",
                     "#aae4df")
  
  # Boxplots of MASE skill scores
  pdf("mase.pdf")
  par(mfrow = c(2, 1))
  boxplot(SS$upper$mase[BP_met], main = "MASE upper time series", 
          col = custom_colors, ylim = c(-80,80))
  abline(h=0,lty=3)
  boxplot(SS$bottom$mase[BP_met], main = "MASE bottom time series", 
          col = custom_colors, ylim = c(-200,200))
  abline(h=0,lty=3)
  dev.off()
  
  # Boxplots of MIS skill scores
  pdf("mis.pdf")
  par(mfrow = c(2, 1))
  boxplot(SS$upper$mis[BP_met], main = "MIS upper time series", 
          col = custom_colors, ylim = c(-150,150))
  abline(h=0,lty=3)
  boxplot(SS$bottom$mis[BP_met], main = "MIS bottom time series", 
          col = custom_colors, ylim = c(-200,200))
  abline(h=0,lty=3)
  dev.off()
  
  # Boxplots of RPS skill scores
  pdf("rps.pdf")
  par(mfrow = c(2,1))
  boxplot(SS$upper$rps[BP_met], main = "RPS upper time series", 
          col = custom_colors, ylim = c(-80,80))
  abline(h=0,lty=3)
  boxplot(SS$bottom$rps[BP_met], main = "RPS bottom time series", 
          col = custom_colors, ylim = c(-200,200))
  abline(h=0,lty=3)
  dev.off()
  
}










