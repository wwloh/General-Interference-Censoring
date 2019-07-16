# 'Randomization inference with general interference and censoring' -----------
# R code for simulations in Section 5 -----------------------------------------

# Load required libraries (and install if needed)
rm(list=ls())
libraries_check <- c("data.table","igraph","survival","stats","msm","mvtnorm")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

n <- 128 # sample size
mc_m <- 2e4 # number of re-randomizations for calculating each p-value
nsimobs <- 10 # number of simulated datasets

# all possible simulation settings
sim_settings <- expand.grid(model_sim=c("add"),
                            network_sim=c("LH","PA"),
                            m_sim=c(124,96,64,32),
                            failrate_sim=0:1)
sim_settings <- sim_settings[rep(1:nrow(sim_settings),each=250),]
# # for assessing power/coverage
power <- FALSE
if (power) {
  sim_settings <- sim_settings[sim_settings$m_sim==96 & sim_settings$failrate_sim==1,]
  sim_settings <- sim_settings[rep(1:nrow(sim_settings),each=nsimobs),]
  nsimobs <- 1
}

# initialize for parallel MC jobs ---------------------------------------------
args <- 13
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
}
(seed <- as.integer(args[1]))
filename <- paste0("sims-withcens-",ifelse(power,"power-",""),seed,".Rdata")


(m <- sim_settings[seed,"m_sim"]) # number of individuals assigned to Z=1
(model <- sim_settings[seed,"model_sim"]) # assumed causal model
(network <- sim_settings[seed,"network_sim"]) # interference structure
(failrate <- sim_settings[seed,"failrate_sim"]) # censoring rate

# adjacency matrix
OneInterfMatrix <- function(n, mean_neighbors) {
  sapply(1:n, function(i) {
    interf_siz <- rpois(n = 1, lambda = mean_neighbors)
    interf_idx <- sample(x = (1:n)[-i], size = interf_siz, replace = FALSE)
    interf_set <- rep(0L, n)
    interf_set[interf_idx] <- 1L
    return(interf_set)
  })
}
## preferential attachment network
OneInterfMatrix_PA <- function(n, mean_neighbors) {
  g <- sample_pa(n=n, m=mean_neighbors/2, directed=FALSE) #half due to symmetry
  adjmat <- as.matrix(as_adjacency_matrix(g))
  return(adjmat)
}
if (network=="LH") {
  a_mtx <- OneInterfMatrix(n,16)
} else {
  a_mtx <- OneInterfMatrix_PA(n,16)  
}

all(diag(a_mtx) == 0L) # check diagonal entries are all zero
a_mtx_rowsums <- rowSums(a_mtx) # interference set sizes
a_mtx_rowsums.nonzero <- sapply(rowSums(a_mtx), max, 1)
a_mtx_edgelist <- as.matrix(get.data.frame(graph.adjacency(a_mtx)))

# generate re-randomizations
omega_mc <- sapply(1:mc_m, function(i) {
  z <- rep(0L, n)
  z[sample.int(n, m, replace = FALSE)] <- 1L
  return(z)
})
unique(colSums(omega_mc))==m # check m is fixed under complete randomization

# calculate values of t and g for each re-randomization z
omega_mc_T <- crossprod(t(a_mtx), omega_mc)
omega_mc_G <- apply(omega_mc_T, 2, function(t) t/a_mtx_rowsums.nonzero)
rm(a_mtx_rowsums.nonzero)

# Step 0 ----------------------------------------------------------------------
delta_true <- 0.7
tau_true <- 2.8; if (model=="BFP") {tau_true <- 0.1}

# # true uniformity trial failure times
mu_true <- 4.5
sigma_true <- .25
# # correlated uniformity failure times between individuals
a_mtx_std <- a_mtx/a_mtx_rowsums # Aij/A_i
all(rowSums(a_mtx_std)==1) # check all row sums equal 1
a_mtx_std <- a_mtx_std* # include noise s.t. corr. matrix is p.d.
  runif(nrow(a_mtx_std)*ncol(a_mtx_std),min=.9,max=1)
corr_mat <- a_mtx_std+t(a_mtx_std)
diag(corr_mat) <- 1 # set diagonal entries equal 1
det(corr_mat) # check positive-definite
y_unif_log <- rmvnorm(1,mean=rep(mu_true,n),sigma=corr_mat*sigma_true^2)[1,]
y_unif <- exp(y_unif_log)

# # true null for checking Type I error
dt.H0 <- cbind(delta_true,tau_true,deparse.level=0) 

if (power) {
  # Values of delta and tau for hypothesis tests
  d.H0vals <- seq(from=0.4,to=1.0,by=0.05)
  t.H0vals <- seq(from=0.8,to=4.0,by=0.2)
  dt.H0 <- expand.grid(d.H0vals,t.H0vals)
  dt.H0 <- as.matrix(dt.H0)
  colnames(dt.H0) <- NULL
}

# Load helper R functions
source("RandomizationInference_Interference_Censoring.R")

ptm <- proc.time()[3]
results <- list()
for (sim_obs in 1:nsimobs) {
  # Step 1 --------------------------------------------------------------------
  # # observed treatments Z,G
  Z <- omega_mc[,sim_obs]
  G <- omega_mc_G[,sim_obs]
  
  # # determine failure times under assigned treatment
  Ytilde <- y_unif*exp(delta_true*Z+tau_true*G)
  if (model=="BFP") {
    GT <- omega_mc_T[,sim_obs]
    Ytilde <- y_unif*exp(delta_true+log(
      1+(1-Z)*(exp(-delta_true)-1))*exp(-tau_true*tau_true*GT))
  }
  
  # censoring times
  cens_admi <- exp(mu_true+2*sigma_true+tau_true)
  # # Z=1
  omega <- sqrt(1-.25^2)
  # cens_Z1_drop <- exp(rnorm(n,mean=mu_true+tau_true*G,sd=omega))
  cens_Z1_drop.log <- rmvnorm(1,mean=mu_true+tau_true*G,
                              sigma=corr_mat*omega^2)[1,]
  cens_Z1_drop <- exp(cens_Z1_drop.log)
  cens_Z1 <- pmin(cens_Z1_drop,cens_admi)
  # # Z=0
  censfact <- (failrate==0)*.6 + (failrate==1)*1
  cens_Z0 <- cens_admi*censfact

  C_ <- cens_Z1*Z + cens_Z0*(1-Z)
  # observed times
  data <- data.frame(time=pmin(Ytilde,C_), status=Ytilde<=C_, Z=Z)
  
  # Step 2 --------------------------------------------------------------------
  # # Carry out hypothesis tests
  if (model=="BFP") {
    G <- GT
    omega_mc_G <- omega_mc_T
  }
  res <- SurvRonly_pvalues(
    obs_Ys=data$time, obs_Zs=data$Z, obs_g=G,
    mc_Zs=omega_mc, mc_Gs=omega_mc_G,
    noncens=data$status,
    a_mtx_rowsums=a_mtx_rowsums,
    dt_H0=dt.H0,
    model=model)
  
  # save simulations results
  setkey(res)
  results[[sim_obs]] <- cbind(true_d=delta_true, true_t=tau_true, 
                              sim=paste(rownames(sim_settings)[seed],sim_obs,sep="."), 
                              res,
                              propD1_Z1=sum(data$status*Z)/m,
                              propD1_Z0=sum(data$status*(1-Z))/(n-m),
                              m=m,model=model,network=network,failrate=failrate)
  rm(res)
  cat("#### sim number", sim_obs, "/", nsimobs,
      "| mins taken:", round((proc.time()[3]-ptm)/60),
      ", remaining:", round((proc.time()[3]-ptm)/sim_obs*(nsimobs-sim_obs)/60),
      "\n")
}
save(results,file=filename)
q()


# Plots -----------------------------------------------------------------------
source("RandomizationInference_Interference-plotting.R")

# setwd("withcens-corr/")
# setwd("withcens-corr-Crho/")
setwd("withcens-corr-rhoU-Crho/")
myfiles <- list.files()
myfiles <- myfiles[grep("Rdata",myfiles)]

results_all <- NULL
for (filename in myfiles) {
  load(filename)
  results <- rbindlist(results,use.names=TRUE,fill=TRUE)
  setkey(results)
  results_all <- rbind(results_all,results)
  rm(results)
}
results <- results_all; rm(results_all)
setkey(results)

# unique number of sims
results[,lapply(.SD, function(x) length(x)/2),
        by=list(model,network,failrate,m),.SDcols="sim"]

# average proportion of failures in each Z group
avD1 <- results[,lapply(.SD, mean),by=list(model,network,failrate,m),
                .SDcols=c("propD1_Z1","propD1_Z0")]
setkey(avD1)
summary(avD1)

teststats <- c("pv.loglik_survreg_BFPa","pv.km_logr.R","pv.loglik_survreg")
names(teststats) <- c("LRaft_BFP","LogR","LRaft")

# # ECDF plots ----------------------------------------------------------------
sim_settings <- unique(results[,list(model,network,failrate,m)])
setkey(sim_settings)

for (model_sim in unique(sim_settings[,model])) {
  for (network_sim in unique(sim_settings[,network])) {
    for (failrate_sim in unique(sim_settings[,failrate])) {
      pdf(paste0("sims-withcens-corr-",
                 model_sim,"-",network_sim,"-",failrate_sim,".pdf"),
          height=6,width=5)
      par(mfcol=c(2,2))
      for (m_sim in c(32,96,64,124)) {
        # model_sim="add";network_sim="LH";failrate_sim=1;m_sim=124
        results_tmp <- results[
          model==model_sim & network==network_sim &
            m==m_sim & failrate==failrate_sim]
        if (nrow(results_tmp)>0) {
          setkey(results_tmp)
          avD1_tmp <- avD1[
            model==model_sim & network==network_sim &
              m==m_sim & failrate==failrate_sim,
            list(propD1_Z1,propD1_Z0)]
          
          sim_ <- c(#ifelse(failrate_sim==0,0.6,1),
            m_sim,round(avD1_tmp,2))
          names(sim_) <- c(#"k",
            "m","p1","p0")
          sim_ <- paste(names(sim_),sim_,sep="=",collapse=",")
          
          ts_plot<-1:3;if(model_sim=="add"){ts_plot<-2:3}
          # OneECDF(results_plot=results_tmp[y0_meth=="fixedD" & C_meth=="fixedD"],
          #         teststats=teststats[ts_plot],ecdf.main=sim_,
          #         legend_loc="bottomright")
          OneECDF(results_plot=results_tmp[y0_meth=="km" & C_meth=="kmz"],
                  teststats=teststats[ts_plot],ecdf.main=sim_,
                  legend_loc="bottomright")
        } else {
          plot.new()
          # plot.new()
        }
      }
      dev.off()
    }
  }
}

## single setting plots
model_sim="add";network_sim="LH";failrate_sim=0;m_sim=124
results_tmp <- results[
  model==model_sim & network==network_sim &
    m==m_sim & failrate==failrate_sim]
setkey(results_tmp)
ts_plot<-2:3
OneECDF(filename="plot-pvaluesECDF-complete_rand-sims-cens-n_128-m_124-d_07-t_28-5-.pdf",
        results_plot=results_tmp[y0_meth=="fixedD" & C_meth=="fixedD"],
        teststats=teststats[ts_plot],
        legend_loc="bottomright")
OneECDF(filename="plot-pvaluesECDF-complete_rand-sims-cens-n_128-m_124-d_07-t_28-1-.pdf",
        results_plot=results_tmp[y0_meth=="km" & C_meth=="kmz"],
        teststats=teststats[ts_plot],
        legend_loc="bottomright")

# empirical coverage of confidence sets ---------------------------------------
delta_true <- unique(results$true_d)
tau_true <- unique(results$true_t)
(teststats <- names(results)[grepl("pv.",names(results))])
names(teststats) <- c("LogR","LRaft","LRaft_BFPa")
(teststats <- teststats[1:2])

results_sims_ci_list <- results[, lapply(.SD, function(pv) {
  sum(pv>=0.05-.Machine$double.eps*1e1, na.rm=T)/sum(!is.na(pv))
  }),
  by=list(y0_meth,C_meth,delta,tau), # average over sims
  .SDcols=teststats]
setkey(results_sims_ci_list)
# nominal coverage at true parameter values
results_sims_ci_list[abs(delta-delta_true)<1e-2 & abs(tau-tau_true)<1e-2]

(boot_meths <- unique(results[,list(y0_meth,C_meth)]))
for (bb in 1:nrow(boot_meths)) {
  results_sims_all <- results_sims_ci_list[
    y0_meth==boot_meths[bb,y0_meth] & C_meth==boot_meths[bb,C_meth]]
  setkey(results_sims_all)
  
  filename <- paste0("plot-coverage-n_256m96-",bb,".pdf")
  
  OneCoveragePlot(
    filename=filename,delta_true=delta_true,tau_true=tau_true,
    results_plot=results_sims_all,
    teststats=rev(teststats))
}

# power ECDF plots
dts <- expand.grid(delta_true+c(-1,0,1)*.1, tau_true+c(-1,0,1)*.4)
for (dt in 1:nrow(dts)) {
  delta0 <- dts[dt,1]
  tau0 <- dts[dt,2]  
  filename <- paste0("plot-power-n_256m96-",dt,".pdf")
  results_tmp <- results[y0_meth=="km" & C_meth=="kmz" & 
                           abs(delta-delta0)<1e-3 & abs(tau-tau0)<1e-3]
  setkey(results_tmp)
  OneECDF(filename=filename,
          results_plot=results_tmp,
          teststats=teststats,
          ecdf.main=bquote(delta[0] == .(delta0) ~ "," ~ tau[0] == .(tau0)),
          legend_loc="bottomright")
}