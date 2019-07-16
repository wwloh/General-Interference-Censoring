# 'Randomization inference with general interference and censoring' -----------
# R code for plots -----------------------------------------

# Load required libraries (and install if needed)
libraries_check <- c("Hmisc", "ggplot2","gridExtra","directlabels")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# # ECDF plots ----------------------------------------------------------------
OneECDF <- function(filename=NULL,
                    results_plot,
                    teststats,
                    ecdf.main=NULL,
                    legend_loc="topleft") {
  if (!is.null(filename)) {
    pdf(filename,width=3.5,height=4)
  }
  if (is.null(ecdf.main)) {
    par(mar=c(4,4,2,1)+0.1) # when there is no title  
  }
  
  Ecdf(seq(0,1,length.out=nrow(results_plot)),
       lty=0,subtitles=FALSE,main=ecdf.main,cex.main=1,
       xlab="p",ylab=bquote("Proportion" <= "p"))
  lines(0:1,0:1,col="grey80")
  
  setkey(results_plot)
  for (pp in 1:length(teststats)) {
    Ecdf(as.numeric(unlist(results_plot[,teststats[pp],with=FALSE])),
         lwd=2,subtitles=FALSE,add=TRUE,
         lty=length(teststats)+1-pp,col=length(teststats)+1-pp)
  }
  
  legend(legend_loc, legend=names(teststats),#xpd = TRUE,horiz=TRUE,
         cex=(length(teststats))^{-.5},bty="n",lwd=2,
         lty=length(teststats):1, col=length(teststats):1)
  if (!is.null(filename)) {
    dev.off()
  }
}

# # Coverage plots for each (delta,tau) tested --------------------------------
OneCoveragePlot <- function(filename=NULL,
                            delta_true,tau_true,
                            results_plot,
                            teststats) {
  # # One plot for each test statistic
  plots_cis <- list()
  for (pp in 1:length(teststats)) {
    sims_cis <- results_plot[, c("delta","tau",teststats[pp]), with=FALSE]
    setnames(sims_cis, ncol(sims_cis), "Coverage")
    setkey(sims_cis)
    
    onep <- ggplot(sims_cis, aes(x=delta,y=tau),
                   xlim = sims_cis[, range(delta)+c(-1,1)*.1],
                   ylim = sims_cis[, range(tau)+c(-1,1)*.1])+
      xlab(expression(delta))+
      ylab(expression(tau))+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle= element_text(hjust = 0.5))
    if (0 >= sims_cis[, min(delta)]) {
      onep <- onep+geom_vline(xintercept=0,colour="grey80",size=.25)
    }
    if (0 >= sims_cis[, min(tau)]) {
      onep <- onep+geom_hline(yintercept=0,colour="grey80",size=.25)
    }
    
    onep <- onep +
      labs(title=names(teststats)[pp]) +
      geom_vline(xintercept=delta_true,colour="red",size=.25)
    if (sims_cis[, tau_true <= max(tau)]) {
      onep <- onep + geom_hline(yintercept=tau_true,colour="red",size=.25)
    }
    
    if (sims_cis[,max(Coverage)>.05]) {
      onep <- onep +
        geom_point(aes(fill = Coverage>=(0.95-1e-2)),
                   shape=22,stroke=.5, size=1e3/nrow(results_plot), 
                   color="grey80")+
        scale_fill_manual(values=c("white", "black"), guide=F)+
        geom_contour(aes(z=Coverage, color = ..level..),
                     breaks=c(.1,.5,.8), size=.3)+
        scale_colour_gradient(na.value="white", low="black", high="black", guide=F)
      onep <- direct.label(onep,method=list("top.pieces",cex=.6))
    } else {
      onep <- onep + geom_point(color = "grey80",shape=0,
                                size=1e3/nrow(results_plot))
    }
    plots_cis[[pp]] <- onep
  }
  
  px <- length(plots_cis); py=1
  if (!is.null(filename)) {
    pdf(filename,width=3.5*px,height=4*py)
    grid.arrange(grobs=plots_cis,nrow=py)
    dev.off()
  } else {
    grid.arrange(grobs=plots_cis,nrow=py)
  }
}


# # Plots for one confidence region -------------------------------------------
OneConfidenceRegion <- function(filename=NULL,
                                delta_true,tau_true,
                                results_plot,
                                teststats,
                                alpha) {
  # # One plot for each test statistic
  plots_cis <- list()
  for (pp in 1:length(teststats)) {
    sims_cis <- results_plot[, c("delta","tau",teststats[pp]), with=FALSE]
    setnames(sims_cis, ncol(sims_cis), "pvalue")
    setkey(sims_cis)
    
    onep <- ggplot(sims_cis, aes(x=delta,y=tau),
                   xlim = sims_cis[, range(delta)+c(-1,1)*.1],
                   ylim = sims_cis[, range(tau)+c(-1,1)*.1])+
      xlab(expression(delta))+
      ylab(expression(tau))+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle= element_text(hjust = 0.5))
    if (0 >= sims_cis[, min(delta)]) {
      onep <- onep+geom_vline(xintercept=0,colour="grey80",size=.25)
    }
    if (0 >= sims_cis[, min(tau)]) {
      onep <- onep+geom_hline(yintercept=0,colour="grey80",size=.25)
    }
    
    onep <- onep +
      labs(title=names(teststats)[pp]) +
      geom_vline(xintercept=delta_true,colour="grey20",size=.25) +
      geom_hline(yintercept=tau_true,colour="grey20",size=.25)
    
    if (sims_cis[,max(pvalue)>alpha]) {
      onep <- onep +
        geom_point(aes(fill = pvalue>=alpha),
                   shape=22,stroke=.5, size=2.5, color="grey80")+
        scale_fill_manual(values=c("white", "black"), guide=F)+
        geom_contour(aes(z=alpha, color = ..level..),
                     breaks=c(.05,.5,.8), size=.3)+
        scale_colour_gradient(na.value="white", low="black", high="black", guide=F)
      onep <- direct.label(onep,method="top.pieces")
    } else {
      onep <- onep + geom_point(color = "grey80",shape=0,size=1.5)
    }
    plots_cis[[pp]] <- onep
  }
  
  px <- length(plots_cis); py=1
  if (!is.null(filename)) {
    pdf(filename,width=3.5*px,height=4*py)
    grid.arrange(grobs=plots_cis,nrow=py)
    dev.off()
  } else {
    grid.arrange(grobs=plots_cis,nrow=py)
  }
}
