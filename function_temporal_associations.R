# Functions for multispatial CCM and S-map to obtain significatn temporal associations

# Detrend data

detrend_data <- function(series_all){
  series <- series_all$abond
  data <- data.frame(x=1:length(series), y=series)
  fit <- lm(y~x, data=data)
  fit_summary <- summary(fit)
  p_value <- fit_summary$coefficient[, 4]['x']
  
  if(!is.na(p_value)){
    if (p_value < 0.05){
      series_all$abond <- as.numeric(series - fit$fitted.values)
    }
  }

  return(series_all)
}

# Multispatial CCM

multisp_CCM_intern <- function(x, niter, i){
  
  x_multi_ccm <- droplevels(na.omit(x[,c(1,2,to_test[,i])]))
  
  if(nrow(x_multi_ccm)>500){
    while(nrow(x_multi_ccm)>500){
      lev <- levels(as.factor(x_multi_ccm$code_square))
      x_multi_ccm <- droplevels(x_multi_ccm[as.factor(x_multi_ccm$code_square) != lev[length(lev)],])
    }
  }
  
  if(length(levels(as.factor(x_multi_ccm$code_square)))>1 & nrow(x_multi_ccm)>=50){
    
    x_multi_ccm <- ddply(x_multi_ccm, .(code_square), .fun=add_row, .parallel = F)
    Accm <- x_multi_ccm[-nrow(x_multi_ccm),3]
    Bccm <- x_multi_ccm[-nrow(x_multi_ccm),4]
    
    maxE <- 5
    
    Emat <- matrix(nrow=maxE-1, ncol=2)
    
    colnames(Emat) <- c("A", "B")
    
    for(E in 2:maxE) {
      Emat[E-1,"A"] <- SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
      Emat[E-1,"B"] <- SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
    }
    
    E_A <- which.max(na.omit(Emat[,1]))+1
    E_B <- which.max(na.omit(Emat[,2]))+1
    
    if(length(E_A)>0 & length(E_B)>0){
      signal_A_out <- SSR_check_signal(A=Accm, E=E_A, tau=1, predsteplist=1:10)
      signal_B_out <- SSR_check_signal(A=Bccm, E=E_B, tau=1, predsteplist=1:10)
      
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_A <- CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_boot_B <- CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_A,CCM_boot_B)
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        
        CCM_boot_A <- CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_A,CCM_boot_A)
        CCM_significance_test[2] <- 1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_B <- CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_B,CCM_boot_B)
        CCM_significance_test[1] <- 1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        CCM_significance_test <- c(1,1)
      }
    }
  }else{
    CCM_significance_test <- c(1,1)
  }
  
  result <- data.frame(a_cause_b=CCM_significance_test[1], # a = species 1, b = species 2
                       b_cause_a=CCM_significance_test[2])
  names(result) <- c(paste0(names(x_multi_ccm)[3],"_cause_",names(x_multi_ccm)[4]),
                     paste0(names(x_multi_ccm)[4],"_cause_",names(x_multi_ccm)[3]))
  return(result)
}

multisp_CCM_intern2 <- function(x, niter, i){
  tryCatch(multisp_CCM_intern(x, niter, i),
           error=function(e) NA)
}

multisp_CCM <- function(x, niter=100){
  
  x <- droplevels(x[,c("code_sp","code_square","year","abond")])
  x <- dcast(x, code_square+year ~ code_sp, value.var="abond")
  
  if(ncol(x)>3){
    to_test <- combn(3:ncol(x),2)
    
    library(parallel)
    
    no_cores <- detectCores() - 1
    
    cl <- makeCluster(no_cores, type = "FORK")
    clusterExport(cl, c("multisp_CCM_intern","to_test","x","niter"), envir=environment())
    
    result_tot <- parLapply(cl, c(1:ncol(to_test)),
                            fun=function(y){
                              multisp_CCM_intern2(x = x,
                                                 niter = niter,
                                                 i = y)
                            })
    
    stopCluster(cl)
    
    
    res_multiCCM <- cbind(unlist(result_tot))
    
    res_multiCCM <- data.frame(species = row.names(res_multiCCM), rho = res_multiCCM[,1])
    
  }else{
    res_multiCCM <- data.frame(species = NA, rho = NA)
  }
  
  
  return(res_multiCCM)
}

# S-map

smap_fun <- function(block, res_ccm){
  
  #block <- droplevels(df_press4[df_press4$Species=="Alauda arvensis" & df_press4$country=="Austria",c("Abd_std","temp_std","urb_std","hico_std","forest_std")])
  block <- block[,c("Abd_std","temp_std","urb_std","hico_std","forest_std")]
  sub_res_ccm <- res_ccm[1,c("temp_cause_species","urb_cause_species",
                        "hico_cause_species","forest_cause_species")]
  num_ccm_sig <- which(sub_res_ccm < 0.05)#sub_res_ccm[1,sub_res_ccm < 0.05]
  sub_res_ccm <- sub_res_ccm[num_ccm_sig]
  to_keep <- order(unlist(sub_res_ccm))
  if(res_ccm[1,"E_A"]<=length(to_keep)){
    to_keep <- to_keep[1:res_ccm[1,"E_A"]]
  }
  sub_res_ccm <- sub_res_ccm[to_keep]
  
  names(block) <- sub("_.*","",names(block))
  block_sub <- block[,c("Abd",sub("_.*","",names(sub_res_ccm)))]
  block_sub <- data.frame(apply(block_sub,2,Zscore))
  
  # Determine the best theta
  theta.examined <- seq(0, 10, by = 0.1)
  th.test <-
    pforeach(
      i      = theta.examined,
      .c     = rbind,
      .cores = 7
    )({
      th.test0 <- block_lnlp(block_sub, method = "s-map", tp = 1,columns=c(1:dim(block_sub)[2]),
                             target_column = 1,theta  = i,
                             silent = T, num_neighbors = 0)
    })
  
  best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']
  
  # Perform multivariate S-map to quantify interaction strengths
  smapc.res <- block_lnlp(block_sub, method = "s-map", tp = 1, columns=c(1:dim(block_sub)[2]),
                          target_column = 1, theta  = best.th,
                           num_neighbors = 0, silent = T,
                          save_smap_coefficients = T)
  smapc.tmp <- data.frame(smapc.res[[1]]$smap_coefficients)
  
  colnames(smapc.tmp) <- c(colnames(block_sub), "Constant")
  
  rho <- smapc.res[[1]]$stats$rho
  
  n_pred <- smapc.res[[1]]$stats$num_pred
  t <- rho*sqrt(n_pred-2)/sqrt(1-rho^2)
  if (t >= 0){
    pvalue <- 1 - pt(t, df = n_pred-2)    
  } else {
    pvalue <- pt(t, df = n_pred-2)
  }
  
  theta <- round(best.th, 2)
  
  return(list(coefficients = smapc.tmp, rho = rho, pvalue = pvalue, theta = theta, res_ccm=res_ccm))
}

smap_fun_signif <- function(data_ccm, res_ccm){
  data_ccm <- droplevels(data_ccm)
  res_ccm <- droplevels(res_ccm[res_ccm$Species==levels(as.factor(data_ccm$Species)),])
  res_smap2 <- res_ccm[,c("temp_cause_species","urb_cause_species",
                          "hico_cause_species","forest_cause_species","E_A")]
  if((res_ccm$temp_cause_species < 0.05 |
     res_ccm$urb_cause_species < 0.05 |
     res_ccm$hico_cause_species < 0.05 |
     res_ccm$forest_cause_species < 0.05) & sum(abs(diff(na.omit(data_ccm$urb_std))))!=0){
    res_smap <- smap_fun(data_ccm,res_smap2)
  }else{
    res_smap <- list(coefficients = NA, rho = NA, pvalue = NA, theta = NA, res_ccm =NA)
  }

  return(res_smap)
}


# Function for test

# Multispatial CCM

multisp_CCM_test <- function(x, niter){
  
  x <- droplevels(x)
  
  if(length(levels(as.factor(x$country)))>1 & nrow(x)>=50){
    
    x_multi_ccm <- ddply(x, .(country), .fun=add_row, .parallel = F)
    Accm <- x_multi_ccm$Index[-nrow(x_multi_ccm)]
    Bccm <- x_multi_ccm[-nrow(x_multi_ccm),"value"]
      
    maxE <- 4
    
    Emat <- matrix(nrow=maxE-1, ncol=2)
      
    colnames(Emat) <- c("A", "B")
      
    for(E in 2:maxE) {
      Emat[E-1,"A"] <- SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
      Emat[E-1,"B"] <- SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
    }
      
    E_A <- which.max(na.omit(Emat[,1]))+1
    E_B <- which.max(na.omit(Emat[,2]))+1
      
    if(length(E_A)>0 & length(E_B)>0){
      signal_A_out <- SSR_check_signal(A=Accm, E=E_A, tau=1, predsteplist=1:10)
      signal_B_out <- SSR_check_signal(A=Bccm, E=E_B, tau=1, predsteplist=1:10)
        
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
       summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
          
        CCM_boot_A <- CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_boot_B <- CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_A,CCM_boot_B)
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
          
        CCM_boot_A <- CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_A,CCM_boot_A)
        CCM_significance_test[2] <- 1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_B <- CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_B,CCM_boot_B)
        CCM_significance_test[1] <- 1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        CCM_significance_test <- c(1,1)
      }
    }
      
    result <- data.frame(a_cause_b=CCM_significance_test[1], # a = species, b= pressure
                           b_cause_a=CCM_significance_test[2])
    names(result) <- c(paste0("species_cause_pressure"),
                         paste0("pressure_cause_species"))

    
  }else{result <- data.frame(species_cause_pressure=NA,pressure_cause_species=NA)}
  
return(result)
}
