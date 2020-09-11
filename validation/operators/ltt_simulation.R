library(ggplot2)
library(Hmisc)
Ne = 1.0;
lambda=1.0;
psi = 0.1;
rho=0.1;
tau=1.0;
mu=0.5;
c_1=abs(sqrt((lambda-mu-psi)^2+4*lambda*psi))
c_2=-(lambda-mu-2*lambda*rho-psi)/c_1
origin=2.0
proposedTree = NULL
geneTree = NULL

# Simple variant: It is now written only for the case of one transmission sample and one sampled ancestor.
# Does not include observed transmission, only non observed.
simulate<-function(transmissionTreeTimes, nSamplesPerTransmissionBranch){
  x = c()
  u =0.0

  for (i in 1:length(transmissionTreeTimes)){
    currentT <<- transmissionTreeTimes[i]
    if (i>1){
      nSamplesPerTransmissionBranch[i] = nSamplesPerTransmissionBranch[i]+nSamplesPerTransmissionBranch[i-1]
    }
    for (j in 1:length(nSamplesPerTransmissionBranch[i])){
      x=append(x, rep(transmissionTreeTimes[i], nSamplesPerTransmissionBranch[i]))
    }

    minCoalTime = Inf;
    minNonObsTransmissionTime = Inf;
    while(nSamplesPerTransmissionBranch[i]>1){
      if (i<length(transmissionTreeTimes) && transmissionTreeTimes[i+1]<currentT)
        break
        minCoalTime = rexp(n=1, rate=nSamplesPerTransmissionBranch[i]*(nSamplesPerTransmissionBranch[i]-1)*0.5/Ne);
        u <<- runif(1)
        if (i<length(transmissionTreeTimes))
          upper = transmissionTreeTimes[i+1]
        else
          upper = 1000#origin
        if (E(currentT) >= 0 && E(upper)<=0 || E(currentT) <= 0 && E(upper) >= 0)
          minNonObsTransmissionTime = uniroot(E, lower=currentT, upper=upper, tol=10^(-16))$root-currentT;

        if (i<length(transmissionTreeTimes) && transmissionTreeTimes[i+1] <= min(minCoalTime, minNonObsTransmissionTime)+currentT )
          break
        if (minCoalTime < minNonObsTransmissionTime && (i+1>length(transmissionTreeTimes) || currentT+minCoalTime < transmissionTreeTimes[i+1])){
          nSamplesPerTransmissionBranch[i]=nSamplesPerTransmissionBranch[i]-1
          currentT <<- currentT+minCoalTime;
          x = append(x, rep(currentT,nSamplesPerTransmissionBranch[i]))
        } else if(minNonObsTransmissionTime < Inf &&  (i+1>length(transmissionTreeTimes) || currentT+minNonObsTransmissionTime < transmissionTreeTimes[i+1])){
          duplicateT = currentT;
          while (duplicateT<currentT+tau && nSamplesPerTransmissionBranch[i]>1){
            coalTime = rexp(n=1, rate=nSamplesPerTransmissionBranch[i]*(nSamplesPerTransmissionBranch[i]-1)*0.5/Ne);
            if (duplicateT+coalTime < currentT+tau){
              nSamplesPerTransmissionBranch[i]=nSamplesPerTransmissionBranch[i]-1
              duplicateT = duplicateT+coalTime;
            } else
              break;
          }
          currentT <<- currentT+minNonObsTransmissionTime;
            x = append(x, rep(currentT,nSamplesPerTransmissionBranch[i]))
      } else
        next
    }
  }
  return(x)
}

# get time of the next hidden event
E <-function(t){
  E_0 = (-log(((c_2-1)*exp(-c_1*currentT)-c_2-1)^2)-c_1*currentT+(lambda+mu+psi)*currentT)/2*lambda
  E_1 = (-log(((c_2-1)*exp(-c_1*t)-c_2-1)^2)-c_1*t+(lambda+mu+psi)*t)/2*lambda

  result = E_1-E_0+log(u)
  return(result)
}

df=data.frame("Time"=numeric(), "Lin"=numeric(), "N"=numeric())
n_ = 1000
for (s in 1:n_ ){
  x=simulate(c(0.0, 0.01), c(1, 17))
  df_s=data.frame(table(x))
  colnames(df_s) = c("Time", "Lin")
  df_s$N=s
  df_s$Time = as.numeric(as.character(df_s$Time))
  df = rbind(df, df_s)
}

maxTime = max(df$Time)

ddf=data.frame("Time"=numeric(n_*100), "Lin"=numeric(n_*100), "N"=numeric(n_*100))
j=1
for (s in 1:n_){
  df_n = df[which(df$N==s),]
  ddf$Time[j]=head(df_n,1)$Time#rbind(ddf, head(df_n,1))
  ddf$Lin[j]=head(df_n,1)$Lin
  ddf$N[j]=head(df_n,1)$N
  i=1
  t_n = df_n$Time[i]
  t_n_max = tail(df_n, 1)$Time
  for (t in seq(0, maxTime, length.out=100)){
    if (t==0)
      next
    j=j+1
    if (t <= t_n_max){
      if(t >= t_n & t < df_n$Time[i+1]){
        ddf$Time[j]=t
        ddf$Lin[j]=df_n[which(df_n$Time==t_n), ]$Lin
        ddf$N[j]=s
      }
        # ddf = rbind(ddf, data.frame("Time"=t, "Lin"=df_n[which(df_n$Time==t_n), ]$Lin, "N"=s))
      else{
        i = i + 1
        t_n = df_n$Time[i]
        ddf$Time[j]=t
        ddf$Lin[j]=df_n[which(df_n$Time==t_n), ]$Lin
        ddf$N[j]=s
        # ddf = rbind(ddf, data.frame("Time"=t, "Lin"= df_n[which(df_n$Time==t_n), ]$Lin, "N"=s))
      }
    }
    else{
      ddf$Time[j]=t
      ddf$Lin[j]=df_n[which(df_n$Time==t_n_max), ]$Lin
      ddf$N[j]=s
      # ddf=rbind(ddf, data.frame("Time"=t, "Lin"=df_n[which(df_n$Time==t_n_max),]$Lin, "N"=s))
    }
  }
}
save(df, ddf, file = "simulated.RData")

p <- ggplot(ddf, aes(Time, Lin)) +
      geom_line(size=0.2, aes(group=N))+
      stat_summary(fun="mean", geom="line", size = 1.5) +
      stat_summary(fun=function(z) { quantile(z,probs=0.025) }, geom="line", size = 0.5, color="red") +
      stat_summary(fun=function(z) { quantile(z,probs=0.975) }, geom="line", size = 0.5, color="red")

ggsave("simulation.pdf", p)
# #  data.frame(table(x))
# likelihood<-function(geneTreeTable){
#   logP=0.0
# 
#   for (int in 1:(nrow(geneTreeTable)-1)){
#     # intervals
#     logP = logP - geneTreeTable[int,2]*(geneTreeTable[int,2]-1)/(2.0*Ne)*(geneTreeTable[int+1,1] - geneTreeTable[int,1])
#     logP = logP - (1-g_ij (geneTreeTable[int,2], geneTreeTable[int,2]))*lambda*integral_p0(geneTreeTable[int,1], geneTreeTable[int+1,1])
#     
#     if (int>1){
#       if (geneTreeTable[int-1,2]-geneTreeTable[int,2]==1){
#         # bifurcation
#         bif = 1.0/Ne + g_ij(geneTreeTable[int-1,2], geneTreeTable[int,2])*lambda*p_0(geneTreeTable[int,1])
#         logP = logP + log(bif)
#       } else if (geneTreeTable[int-1,2]-geneTreeTable[int,2]>1){
#         # multifurcation
#         mult=g_ij(geneTreeTable[int-1,2], geneTreeTable[int,2])*lambda*p_0(geneTreeTable[int,1])
#         logP=logP + log(mult)
#       }
#     }
#   }
#   return(logP)
# }
# 
# 
# 
# operator<-function(geneTreeTable){
#   logP=0
#   # choose timepoint, skip past sampling times (first two) and root
#   if (length(3:(which(geneTreeTable[,2] == 1))[2])==1){
#     int_s=3
#   } else{
#     int_s = sample(3:(which(geneTreeTable[,2] == 1))[2],1)
#   }
#   # if the root was already reached before, skip
#   if (geneTreeTable[int_s-1,2]==1)
#     return(-Inf)
#   
#   if (runif(1) < 0.5){
#     #decrease
#     diff_down = geneTreeTable[int_s,2]-1
#     if (diff_down==0)
#       return(-Inf)
#     by = 1 #sample(1:diff_down, 1)
#     for (int_d in int_s:99){
#       # if (by > 1 && geneTreeTable[int_d,1]> origin)
#       #   break
#       if(geneTreeTable[int_d,2]==1)
#         break
#       val = geneTreeTable[int_d,2] - by
#       if (val< 1){
#         val=1
#       }
#       geneTreeTable[int_d,2] = val
#     }
#     # diff_up = geneTreeTable[int_s-1,2]-geneTreeTable[int_s,2]
#     # if(diff_up<0){
#     #   print(geneTreeTable[int_s-1,2])
#     # }
#     # logP=log(diff_up/diff_down)
#     proposedTree <<-geneTreeTable
#     return(logP)  
#   } else {
#     # increase
#     diff_up = geneTreeTable[int_s-1,2]-geneTreeTable[int_s,2]
#     if (diff_up<=0)
#       return(-Inf)
#     by = 1 #sample(1:diff_up, 1)
#     for (int_d in int_s:99){
#       val = geneTreeTable[int_d,2] + by
#       geneTreeTable[int_d,2] = val
#       if (val == 1+by)
#         break
#     }
#     # diff_down = geneTreeTable[int_s,2]-1
#     # if(diff_up<0)
#     #   print("yoyo")
#     # logP=log(diff_down/diff_up)
#     proposedTree <<-geneTreeTable
#     return(logP) 
#   }
# }
# 
# set.seed(42)
# n_iter=100000
# mcmc = data.frame("Time"=numeric(), "Lin"=numeric(), "N"=numeric())
# geneTree = head(ddf,100)
# old_likelihood = likelihood(geneTree)
# for (s in 1:n_iter){
#   logHR = operator(geneTree)
#   if (logHR==-Inf){
#     geneTree$N=s
#     colnames(geneTree) = c("Time", "Lin", "N")
#   } else{
#     new_likelihood = likelihood(proposedTree)
#     logAlpha = new_likelihood - old_likelihood + logHR
#     if (logAlpha >=0 || runif(1)<exp(logAlpha)){
#       # accept
#       geneTree <- proposedTree
#       geneTree$N=s
#       colnames(geneTree) = c("Time", "Lin", "N")
#       old_likelihood = new_likelihood
#       mcmc = rbind(mcmc, geneTree)
#     } else{
#       # reject
#       geneTree$N=s
#       colnames(geneTree) = c("Time", "Lin", "N")
#     }
#   } 
#   # if(s%%1000==0){
#     
#   # }
#   
# }
# 
# save(mcmc, file = "mcmc.RData")
# 
# ggplot(mcmc, aes(Time, Lin)) +
#       geom_line(size=0.2, aes(group=N), color="grey")+
#       stat_summary(fun="mean", geom="line", size = 1.5) +
#       stat_summary(fun=function(z) { quantile(z,probs=0.025) }, geom="line", size = 0.5, color="red") +
#       stat_summary(fun=function(z) { quantile(z,probs=0.975) }, geom="line", size = 0.5, color="red")
# 

integral_p0 <- function(t0, t1){
  ans= (1.0/(2.0*lambda)) * ((t1-t0)*(mu+psi+lambda-c_1)+2.0*log(((c_2-1)*exp(-c_1*t0)-c_2-1) / ((c_2-1)*exp(-c_1*t1)-c_2-1)))
  return(ans)
}

p_0 <- function(t){
  p0=(lambda+mu+psi+c_1*( (exp(-c_1*t)*(1-c_2)-(1+c_2))/(exp(-c_1*t)*(1-c_2)+(1+c_2))))/(2.0*lambda)
  return(p0)
}

g_ij <- function(i, j){
  ans = 0.0
  for (k in j:i){
    ans = ans + (2.0*k-1)*((-1.0)^(k-j))*f_1(j, k-1)*f_2(i, k)*exp(-(k*(k-1)*tau*0.5)/Ne) / (factorial(j)*factorial(k-j)*f_1(i, k))
  }
  return(ans)
}

f_1 <- function(x, y){
  ans = 1.0
  if(y==0)
    return(ans)
  for(f1 in 0:(y-1)){
    ans = ans * (x+f1)
  }
  
  return(ans)
}

f_2 <- function(x, y){
  ans=1.0
  for(f2 in 0:(y-1)){
    ans = ans * (x-f2)
  }
  return(ans)
}
