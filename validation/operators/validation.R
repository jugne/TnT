library(ggplot2)
library(grid)
library(gridExtra)
library(coda)

# clear workspace
rm(list = ls())
# Set the directory to the directory of the file
# this.dir <- dirname(parent.frame(2)$ofile)
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)


doCompare <- function(dsFileName, mcmcFileName) {
  dfs <- read.table(dsFileName, header=T)
  df1 <- read.table( mcmcFileName, header=T)
  for (per in c(0.1,0.5,1.0)){
    df<-df1[(1:ceiling(per*dim(df1)[1])),]
    
    # Remove 10% for burnin
    N <- dim(df)[1]
    df <- df[-(1:ceiling(0.1*N)),]
    
    # Tree age/length plot
    pdf(gsub(".log", ".pdf", paste0(per,"_",mcmcFileName)), width = 10, height = 10) 
    
    
    par(mfcol=c(2,1))
    
    maxLength <- max(quantile(df$geneTree.treeLength, probs=0.99),
                     quantile(dfs$geneTree.treeLength, probs=0.99),
                     quantile(df$geneTree.height, probs=0.99),
                     quantile(dfs$geneTree.height, probs=0.99))
    
    maxDensity <- max(density(dfs$geneTree.height)$y,
                      density(df$geneTree.height)$y,
                      density(dfs$geneTree.treeLength)$y,
                      density(df$geneTree.treeLength)$y)
    
    plot(density(df$geneTree.height), 'l', col='red', lwd=2, lty=2,
         xlim=c(0,maxLength), ylim=c(0,maxDensity),
         xlab="Statistic", ylab="Density",
         main="")
    lines(density(dfs$geneTree.height), col='blue', lwd=2, lty=2)
    lines(density(df$geneTree.treeLength), col='red', lwd=2)
    lines(density(dfs$geneTree.treeLength), col='blue', lwd=2)
    legend("topright", inset=0.05,
           c("MCMC", "Direct simulation","Gene tree length", "Gene tree height"),
           lty=c(1,1,1,2), lwd=2, col=c("red","blue","black","black"))
    
    # counts <- table(dfs$PolytomyCount)
    # counts <- counts/sum(counts)
    # sim <- data.frame(counts)
    # sim$method = "Direct Simulation"
    # 
    # counts2 <- table(df$PolytomyCount)
    # counts2 = counts2/sum(counts2)
    # inf <- data.frame(counts2)
    # inf$method = "MCMC"
    # rbind(sim, inf)
    
    
    # Kolmogorov-Smirnov test for network statistics
    
    H <- c()
    L <- c()
    l <- length(df$geneTree.height)
    
    # For comparison to run faster, calculations need to be done at fewer ponts.
    # It can be achieved by taking a fraction smaller than 0.01 below.
    nPoints <- round(l*0.0001)
    by <- as.integer(l/nPoints)
    print(l)
    print(nPoints)
    print(by)
    
    for (i in seq(l, 0, by=-by)){
      H <- append(H, log(ks.test(df$geneTree.height[1:(l-i)], dfs$geneTree.height)$statistic))
      L <- append(L, log(ks.test(df$geneTree.treeLength[1:(l-i)], dfs$geneTree.treeLength)$statistic))
    }
    plot(seq(1, length(L)), L, type='l', col='red',  lwd=2,
         xlab=paste("Number of logged iterations (minus burnin) / ", by),
         ylab="Kolmogorov-Smirnov stat. (log scale)", 
         main="Kolmogorov-Smirnov test for network statistics",
         ylim=c(min(L, H), max(L, H)))
    
    lines(seq(1, length(H)), H, col='blue', lwd=2, lty=2)
    legend("topright", inset=0.05,
           c("Gene tree length", "Gene tree height"),
           lty=c(2,2), lwd=2, col=c("red","blue"))
    dev.off()
  }
  
  
}

# doCompare("simulate_geneTree_SPR.log",
#           "test_geneTree_SPR.geneTree.log")
# 
# doCompare("simulate_geneTree_SPR_2spBranches.log",
#           "test_geneTree_SPR_2spBranches.geneTree.log")

doCompare("simulate_geneTree_SPR.log",
          "test_geneTree_SPR2.geneTree.log")



