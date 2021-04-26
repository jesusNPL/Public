
pred1 = "loo_eff_reflectance_canopy"
pred2 = "loo_eff_functional_diversity_canopy"
#flujo3 = "loo_eff_spectra_diversity_canopy"
#flujo4 = "loo_eff_veg_indices_canopy"

ruta = "C:/Users/jpintole/Documents/Sandra/Fluxes/Ecosystem_predictors/out/"


predictores <- c("reflectance", "FD", "SD", "VI")
flux <- c("NEE", "GPP", "Reco",	"Rsoil",	"SOC",	"ET", "WUE")
covs <- c("reflec", "FD")

demonSandra_loo <- function(RDS1, RDS2, RDS3, RDS4, path, fluxNames, covars) { 
  rds1 <-  readRDS(paste0(path, RDS1, ".rds"))
  rds2 <-  readRDS(paste0(path, RDS2, ".rds"))
  #rds3 <-  readRDS(paste0(path, RDS3, ".rds"))
  #rds4 <-  readRDS(paste0(path, RDS4, ".rds"))
  
  ll1 <- list()  #loo por flujo si agrego X predic then x ll 
  ll2 <- list()
  #ll3 <- list()
  #ll4 <- list()
  
  for(j in 1:length(rds1)) { 
    
    ll1[[j]] <- rds1[[j]]$LOO
    ll2[[j]] <- rds2[[j]]$LOO
    #ll3[[j]] <- rds3[[j]]$LOO
    #ll4[[j]] <- rds4[[j]]$LOO
    
  }
  
  nee <- list(ll1[[1]], ll2[[1]]) #pred
  names(nee) <- covars
  gpp <- list(ll1[[2]], ll2[[2]]) 
  names(gpp) <- covars
  
  NEE <- loo_compare(x = nee)
  NEE <- data.frame(print(NEE, simplify = FALSE, digits = 3))
  NEE$Flux <- "NEE"
  NEE$Covariates <- covars
  
  GPP <- loo_compare(x = gpp)
  GPP <- data.frame(print(GPP, simplify = FALSE, digits = 3))
  GPP$Flux <- "GPP" 
  GPP$Covariates <- covars
  
  results <- rbind(NEE, GPP)
  return(results)
}

demonSandra_loo(RDS1 = pred1, RDS2 = pred2, path = ruta, fluxNames = flux[1:2], covars = covs)


