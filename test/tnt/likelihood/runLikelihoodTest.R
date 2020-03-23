library(ggplot2)
# library(colorblindr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
seed = 1656718
set.seed(seed)

for (i in 1:10) {
  setwd(this.dir)
  directory_name = paste0("test_",i)
  Ne = runif(1, 1,10)
  tau = runif(1, 1,10)
  s = sample(1:10000000, 1)
  nGeneSamples = sample(1:10, 1)
  
  mainDir <- this.dir;
  subDir <- directory_name;
  
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }
  
  system(paste0("java -jar ./LikelihoodTest.jar -seed",
               s," -Ne ",Ne," -tau ",tau," -nGeneSamples",nGeneSamples,
               " -calcEvery 10 -runs 10000000"))
  data = read.table("likelihood_test.txt", header=TRUE)
  data$expected_score_exact_log = log10(data$expected_score_exact)
  
  p<- ggplot(data, aes(x=n_run,y=data$expected_score_exact_log)) + geom_point() +
    xlab("Number of simulations") + ylab("Expected empirical score statistic (log scale)")
  
  ggsave(plot=p,paste0(directory_name,".pdf", sep=""))
  
}

data = read.table("likelihood_test.txt", header=TRUE);

data$expected_score_exact_log = log10(data$expected_score_exact)

p<- ggplot(data, aes(x=n_run,y=data$expected_score_exact_log)) + geom_point() +
  xlab("Number of simulations") + ylab("Expected empirical score statistic (log scale)")

ggsave(plot=p,paste0("test_1.pdf", sep=""))

# ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/", fnames[[fn]], "_nrevents.unfit.pdf", sep=""), width=3, height=3)
