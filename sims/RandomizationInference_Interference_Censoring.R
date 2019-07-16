SurvRonlyTestStats <- function(y0_, noncens_, z_, g_, 
                               neighs, aft_dist_="lognormal") {
  # y0_= uniformity trial outcomes under H0
  
  # y_= survival object based on y0_ under H0 and failing indicators (noncens_)
  y_ <- Surv(y0_, event=noncens_)
  
  teststatR <- NULL
  resnames <- NULL
  
  fit_diff <- survdiff(y_ ~ z_)
  teststatR <- c(teststatR, fit_diff$chisq)
  resnames <- c(resnames, "pv.km_logr.R") # KM

  fit_surv <- survreg(y_ ~ z_*g_+neighs,
                      control=survreg.control(maxiter=1e6),
                      score=TRUE,
                      dist=aft_dist_)
  yhats <- predict(fit_surv)
  
  fit_surv_ts <- diff(fit_surv$loglik) # LRaft
  teststatR <- c(teststatR,fit_surv_ts)
  resnames <- c(resnames, "pv.loglik_survreg")
  
  # BFPa model
  g_z0 <- (1-z_)*g_
  fit_surv <- survreg(y_ ~ z_+g_z0+neighs,
                      control=survreg.control(maxiter=1e6),
                      score=TRUE,
                      dist=aft_dist_)
  yhats <- predict(fit_surv)
  
  fit_surv_ts <- 2*diff(fit_surv$loglik) # LRaft_BFPa
  teststatR <- c(teststatR,fit_surv_ts)
  resnames <- c(resnames, "pv.loglik_survreg_BFPa")
  
  names(teststatR) <- resnames
  return(teststatR)
}


# Modular function to impute failure times for censored observations ----------
ImputeFailureTimes <- function(y_o, # min of failure times and censoring times
                               noncens_, # failure indicator
                               impute_num, # number of imputations
                               n_ # sample size
                               ) {
  # survival object
  y_ <- Surv(y_o, event=noncens_)
  # (1) KM estimator
  sfit <- survfit(y_~1L)
  surv_cdf <- data.table(cbind(time=sfit$time,cdfhat=1.0-sfit$surv))
  surv_cdf_max <- max(surv_cdf$cdfhat)
  time_max <- max(surv_cdf$time)

  # (2) Generate multiple imputations (random draws)
  tildeY_MIlist <- lapply(1:impute_num, function(u) {
    # sample failure time from common CDF
    unlist(sapply(1:n_, function(i) {
      tildeY <- y_o[i]
      if (noncens_[i]==0L) {
        vi_lb <- 0.0
        cdf_idx <- findInterval(x=y_o[i],vec=surv_cdf$time)
        if (cdf_idx>0) {
          vi_lb <- surv_cdf$cdfhat[cdf_idx]  
        }
        v_i <- runif(1, min=vi_lb, max=1.0)
        if (v_i < surv_cdf_max) {
          tildeY <- surv_cdf$time[findInterval(x=v_i,vec=surv_cdf$cdfhat)+1]
        } else {
          tildeY <- time_max
        }
      }
      return(tildeY)
    }))
  })
  return(tildeY_MIlist)
}

# Modular function to impute censoring times for rerandomizations -------------
ImputeCensoringTimes <- function(tmpdat, # data frame with 'time','status','z'
                                 omega_mc_z, # MC samples of Z
                                 n_ # sample size
                                 ) {
  
  # for imputing all n x MC censoring times at once 
  # (faster than imputing for each subset by each rerandomization)
  n_mc <- n_*ncol(omega_mc_z)
  
  # survival object
  cens_surv <- Surv(time=tmpdat$time, event=tmpdat$status)
  # conditional KM estimators for Z=0,1 (IPZ)
  cens_km <- lapply(0:1, function(zz) {
    sfit <- survfit(Surv(time=time, event=status)~1L,data=tmpdat[tmpdat$z==zz,])
    surv_cdf <- data.table(cbind(time=sfit$time,cdfhat=1.0-sfit$surv))
    surv_cdf_max <- max(surv_cdf$cdfhat)
    time_max <- max(surv_cdf$time)
    cens_ast_km_MI <- vapply(1:n_mc, function(i) {
      cens_i <- time_max
      ## no truncation for censoring times
      v_i <- runif(1, min=0.0, max=1.0)
      if (v_i < surv_cdf_max) {
        cens_i <- surv_cdf$time[findInterval(x=v_i,vec=surv_cdf$cdfhat)+1]
      }
      return(cens_i)
    }, FUN.VALUE = numeric(1))
    cens_ast_km_MI <- matrix(cens_ast_km_MI,nrow=n_)
    return(cens_ast_km_MI)
  })
  cens_ast_km_MI <- (1L-omega_mc_z)*cens_km[[1]]+omega_mc_z*cens_km[[2]]
  return(cens_ast_km_MI)
}


# Directly in R -----------------------------------------
SurvRonly_pvalues <- function(obs_Ys, obs_Zs, obs_g,
                              mc_Zs, mc_Gs,
                              noncens,a_mtx_rowsums,dt_H0,
                              aft_dist="lognormal",
                              model="add") {

  m <- sum(obs_Zs); n <- length(obs_Zs)
  
  # estimators of censoring time distribution ---------------------------------
  tmpdat <- data.frame(cbind("time"=obs_Ys,"status"=1L-noncens,"z"=obs_Zs))
  C_ast_km <- ImputeCensoringTimes(tmpdat=tmpdat,omega_mc_z=mc_Zs,n_=n)
  rm(tmpdat)
  
  results_R <- list()
  ptm <- proc.time()[3]
  for (dti in 1:nrow(dt_H0)) {
    delta<-dt_H0[dti,1];tau<-dt_H0[dti,2]
    pvalues <- NULL

    # Uniformity trial potential outcomes under hypothesized (delta,tau) ------
    if (model!="BFP") {
      y0unif <- obs_Ys / exp(delta*obs_Zs + tau*obs_g)
    } else {
      y0unif <- obs_Ys / exp(delta+log(
        1+(1-obs_Zs)*(exp(-delta)-1))*exp(-tau*tau*obs_g))
    }
    
    # uniformity failure times under hypothesized (delta,tau) -----------------
    ## imputations based on KM estimator (one per re-randomization)
    tilde_y0unif_ast_km <- ImputeFailureTimes(y_o=y0unif, noncens_=noncens,
                                              impute_num=ncol(mc_Zs), n_=n)
    ## transform to failure times for each re-randomization
    mc_y0d_boot_list <- lapply(seq.int(ncol(mc_Zs)), function(sim_mc) {
      mc_z <- mc_Zs[,sim_mc]
      mc_g <- mc_Gs[,sim_mc]
      if (model!="BFP") {
        mc_caF <- exp(delta*mc_z+tau*mc_g)
      } else {
        mc_caF <- exp(delta+log(1+(1-mc_z)*(exp(-delta)-1))*exp(-tau*tau*mc_g))
      }
      
      # Impute uniformity failure times
      tilde_y0unif_ast <- list()
      tilde_y0unif_ast[[1]] <- tilde_y0unif_ast_km[[sim_mc]]
      names(tilde_y0unif_ast)[1] <- "km"
      
      # Failure times for a particular re-randomization
      tildeY_ast <- lapply(tilde_y0unif_ast, function(tilde_y0unif_ast_) {
        tilde_y0unif_ast_*mc_caF
      })
      
      # Impute censoring times for a particular re-randomization
      C_ast <- list()
      C_ast[[1]] <- C_ast_km[,sim_mc]
      names(C_ast)[1] <- "kmz"
      
      # Combinations of imputed failure times and censoring times
      mc_y0d_boot <- NULL
      for (yy in 1:length(tildeY_ast)) {
        for (cc in 1:length(C_ast)) {
          ## uniformity outcomes and failure indicators for tests
          D_ast <- tildeY_ast[[yy]] <= C_ast[[cc]]
          Y_ast <- pmin(tildeY_ast[[yy]],C_ast[[cc]])
          y0_ast <- Y_ast / mc_caF
          mc_y0d_boot <- c(mc_y0d_boot,list(data.table(
            mc_id=sim_mc+1L, # increase index since first column will be obs Z
            y0_meth=names(tildeY_ast)[yy],
            C_meth=names(C_ast)[cc],
            i=1:n, # needed to maintain order of observations
            y0_ast=y0_ast,
            D_ast=D_ast)))
        }
      }
      mc_y0d_boot.dt <- rbindlist(mc_y0d_boot)
      return(mc_y0d_boot.dt)
    })
    mc_y0d_boot_dt <- rbindlist(mc_y0d_boot_list)
    setkey(mc_y0d_boot_dt); rm(mc_y0d_boot_list)
    
    # unique bootstrap methods
    boot_meths <- unique(mc_y0d_boot_dt[, list(y0_meth,C_meth)])
    boot_meths_nofixedD <- 0 # 0 to include tests holding D fixed
    
    pvalues_allmeth <- NULL
    omega_mc_z <- as.matrix(cbind(obs_Zs,mc_Zs))
    omega_mc_g <- as.matrix(cbind(obs_g,mc_Gs))
    
    for (booti in boot_meths_nofixedD:nrow(boot_meths)) {
      # Test statistics sampling distribution
      tmp <- sapply(seq.int(ncol(omega_mc_z)), function(sim_mc) {
        mc_z <- omega_mc_z[,sim_mc]
        mc_g <- omega_mc_g[,sim_mc]
        if (sim_mc==1 | booti==0) {
          mc_y0 <- y0unif
          mc_d <- noncens
        } else {
          mc_y0d_boot <- mc_y0d_boot_dt[
            data.table(mc_id=sim_mc,boot_meths[booti])]
          mc_y0 <- mc_y0d_boot[,y0_ast]
          mc_d <- mc_y0d_boot[,D_ast]
        }
        SurvRonlyTestStats(
          y0_=mc_y0, noncens_=mc_d,
          z_=mc_z, g_=mc_g,
          neighs=a_mtx_rowsums)
      })
      pvalues <- apply(tmp, 1, function(ts_) {
        mean(ts_[-1]>=ts_[1]-.Machine$double.eps*1e1,na.rm = TRUE)
      })
      
      if (booti==0) {
        boot_meth <- data.table(cbind("fixedD","fixedD"))
        names(boot_meth) <- names(boot_meths)
      } else {
        boot_meth <- boot_meths[booti]
      }
      pvalues_allmeth <- c(pvalues_allmeth, list(c(boot_meth, pvalues)))
      rm(pvalues)
    }
    pvalues <- rbindlist(pvalues_allmeth)
    pvalues <- cbind(delta=delta,tau=tau,pvalues)
    
    results_R <- c(results_R, list(pvalues))
    cat(dti, "/", nrow(dt_H0), "H0s | mins =",
        round((proc.time()[3]-ptm)/60), '\n')
  }
  res <- rbindlist(lapply(results_R,function(x) data.table(rbind(x))),
                   use.names=TRUE,fill=TRUE)
  setkey(res)
  return(res)
}