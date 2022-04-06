library(ggplot2)
library(greekLetters)

# parameters
lambda_ <- 4.0
mu_ <- 1.1
psi_ <- 4.4 # sampling proportion psi/(psi+mu)=0.8
r <- c(0.0, 0.3, 0.5, 1.0)
rho_ <- c(0.5, 0.0, 1.0, 1.0)
origin <- 1.0
seed = c(614010, 183931, 793, 2836)
minSamples <- 1
maxSamples <- 2000
nSims <- 400000

fileName <- c("noR","noRho","allParams_1", "allParams_2") 

wd<-"/Users/jugne/Documents/Source/TnT/validation/mapper/validation/"
setwd(wd)

for (i in 1:4){
  # run simulations and mappings
  x <- paste('java -jar',paste0(wd,"hiddenEventsMapperTestFigures_06042022.jar"), 
             lambda_, mu_, psi_, r[i], rho_[i], origin, seed[i], minSamples, maxSamples, nSims, fileName[i])
  system(x)
  
  tt<- read.table(paste0(wd,fileName[i],".hiddenEventsMap.log"), header = T)
  tt<-tt[,2:3]
  
  
  maxCount <- max(tt$hiddenEventsSim, tt$hiddenEventsMap)+1
  hs <- hist(tt$hiddenEventsSim, plot=F,
             breaks=seq(-0.5,maxCount+0.5,by=1))
  h <- hist(tt$hiddenEventsMap, plot=F,
            breaks=seq(-0.5,maxCount+0.5,by=1))
  
  dt_h <- data.frame(h$mids, h$density)
  dt_hs <- data.frame(hs$mids, hs$density)
  
  colnames(dt_h) <- c("mids", "density")
  colnames(dt_hs) <- c("mids", "density")
  
  p1 <- ggplot() + geom_point(data=dt_h, aes(x=mids, y=density, colour="Mapped Hidden Events"))+ geom_line(data=dt_h, aes(x=mids, y=density, colour="Mapped Hidden Events"))
  p1 <- p1 + geom_point(data=dt_hs, aes(x=mids, y=density, colour="True Hidden Events"))+ geom_line(data=dt_hs, aes(x=mids, y=density, colour="True Hidden Events"))
  p1 <- p1  + theme(plot.title = element_text(hjust = 0.5, size=10), legend.position = c(0.75, 0.65)) + labs(colour = "", x="Statistic", y="Density")
  p1 <- p1 + ggtitle(paste0(greeks("lambda"), "=",lambda_,", ",greeks("mu"), "=",mu_,", ",greeks("psi"),"=",psi_,", ","r=",r[i],", ",greeks("rho"),"=",rho_[i],", origin=",origin))
  p1
  
  ggsave(paste0("/Users/jugne/Documents/Source/TnT/validation/mapper/validation/",fileName[i],".pdf"), p1)
}




