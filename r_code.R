suppressWarnings(suppressMessages(suppressPackageStartupMessages({
library(tidyverse)
library(tableone)
library(MASS)
library(survival)
library(MASS)
library(Epi)
library(readxl)
library(survminer)
library(base)
library(sjPlot) # Prints readable regression summaries & interaction plots
library(sjmisc)
library(sjlabelled)
library(Epi)
library(popEpi)
library(mgcv)
library(survival)
}) ))
n = 5000 
beta1 = -1.5   # set the effect for covariate A   
beta2 =  0.5    # set the effect for covatiate B
 lambdaT = .055  # baseline hazard #50% prevalence
lambdaC = .15   # hazard of censoring
rho <- -0.2  # set the correlatioon to zero to induce independence hre
mu1 <- 1.5; s1 <- .2
mu2 <- 2.8; s2 <- .6
mu <- c(mu1,mu2)        # Mean
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix
r<-100   # the number of simulation
simdata<-list()
for(i in 1:r){
  bvnorm <- mvrnorm(n, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvnorm) <- c("x1","x2")
  x1<-bvnorm[,1]
  x2<-bvnorm[,2]
  t = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2)) # true event time
  C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
  time = pmin(t,C)  #observed time is min of censored and true time
  event<-ifelse(time==t,1,0)
  survdata<-data.frame(time=time, event=event, x1 = x1, x2 = x2)
  simdata[[i]]<-survdata 
}
head(simdata[[r]])  
simulate_data<-function(n=n, r=r,lambdaT=lambdaT,lambdaC=lambdaC,beta1=beta1, beta2=beta2,
                        mu1=mu1, s1=s1, mu2=mu2, s2=s2, rho=rho){
  mu <- c(mu1,mu2) 
  sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix
  simdata<-list()
  for(i in 1:r){
    bvnorm <- mvrnorm(n, mu = mu, Sigma = sigma ) # from MASS package
    colnames(bvnorm) <- c("x1","x2")
    x1<-bvnorm[,1]
    x2<-bvnorm[,2]
    t = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2)) # true event time
    C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
    time = pmin(t,C)  #observed time is min of censored and true time
    event<-ifelse(time==t,1,0)
    survdata<-data.frame(time=time, event=event, x1 = x1, x2 = x2)
    simdata[[i]]<-survdata 
  }
  return(simdata)
}
(prevalence<-sum(simdata[[1]]$event)/length(simdata[[1]]$event)*100) # prevalence
simdata.artificial<-simdata
unpoolcox<-function(r,simdata){
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(i in 1:r){
    survdata<-simdata[[i]]
    ress<-summary(coxph(Surv(time, event)~ x1 + x2, method="breslow", data = survdata))
    beta0vec[i]<-ress$coefficients[1,1]
    SEbeta0vec[i]<-ress$coefficients[1,3]
    beta1vec[i]<-ress$coefficients[2,1]
    SEbeta1vec[i]<-ress$coefficients[2,3]
  }
  result1<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec), median(beta1vec) ),
                     SE = c(median(SEbeta0vec), median(SEbeta1vec) ) )
  result<-data.frame(beta0vec ,beta1vec)
  return(list(estimates = result1, result = result)) # returns median estimates, CI & vector or results
}
cat("Cox estimate")
pools = unpoolcox(r,simdata)
results <- mod_res<- pools$estimates
results
Estbeta = results$estimate[1] ; SEbeta = results$SE[1];
Estgamma = results$estimate[2]; SEgamma = results$SE[2]
beta_int = c(Estbeta-2*SEbeta,Estbeta+2*SEbeta)
gamma_int = c(Estgamma-2*SEgamma,Estgamma+2*SEgamma)
cbind(
c(results$estimate[1],beta_int),
c(results$estimate[2],gamma_int)
)
results = pools$result
cat("beta coverage")
sum(ifelse(results$beta0vec>beta_int[1]&results$beta0vec<beta_int[2],1,0))/r
cat("gamma coverage")
sum(ifelse(results$beta1vec>gamma_int[1]&results$beta1vec<gamma_int[2],1,0))/r
cat("mean absolute bias")
mean(abs(beta1-pools$result[,1]))
mean(abs(beta2-pools$result[,2]))
cat("relative efficiency")
emp.var.beta1<-sum((pools$result[,1]-mean(pools$result[,1]))^2)/(length(pools$result[,1])-1)
mod.var.beta1<- mod_res$SE[1]^2
mod.var.beta1/emp.var.beta1
emp.var.beta2<-sum((pools$result[,2]-mean(pools$result[,2]))^2)/(length(pools$result[,2])-1)
mod.var.beta2<- mod_res$SE[2]^2
mod.var.beta2/emp.var.beta2
survfit.data<-simdata[[1]]
surv.fit<-survfit(Surv(time,event)~1,data=survfit.data)
plot(surv.fit, lwd = 3, cex.lab=1.5, cex.axis=1.25, xlab = "time", ylab = "Survival Prob")
survfit.data$datatype<-"Individual"
KM1 = survfit.data[,c(1,2,5)]
unpoolclogit<-function(r,m,simdata){
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    eventid<-which(event==1)[tempfoo]
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
      }
    }
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    ress<-summary(clogit(formula = event~x1+x2+strata(risksetid), data = survdata))
    beta0vec[z]<-ress$coefficients[1,1]
    SEbeta0vec[z]<-ress$coefficients[1,3]
    beta1vec[z]<-ress$coefficients[2,1]
    SEbeta1vec[z]<-ress$coefficients[2,3]
  }
  result1<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  result<-data.frame(beta0vec ,beta1vec)
  return(list(estimates = result1, result = result, pooldata = survdata))
}
cat("Unpool conditional logistic")
pools=unpoolclogit(r,m=2, simdata)
results <-mod_res<- pools$estimates
results
Estbeta = results$estimate[1] ; SEbeta = results$SE[1];
Estgamma = results$estimate[2]; SEgamma = results$SE[2]
beta_int = c(Estbeta-2*SEbeta,Estbeta+2*SEbeta)
gamma_int = c(Estgamma-2*SEgamma,Estgamma+2*SEgamma)
cbind(
c(results$estimate[1],beta_int),
c(results$estimate[2],gamma_int)
)
results=pools$result
cat("beta coverage")
sum(ifelse(results$beta0vec>beta_int[1]&results$beta0vec<beta_int[2],1,0))/r
cat("gamma coverage")
sum(ifelse(results$beta1vec>gamma_int[1]&results$beta1vec<gamma_int[2],1,0))/r
cat("mean absolute bias")
mean(abs(beta1-pools$result[,1]))
mean(abs(beta2-pools$result[,2]))
cat("relative efficiency")
emp.var.beta1<-sum((pools$result[,1]-mean(pools$result[,1]))^2)/(length(pools$result[,1])-1)
mod.var.beta1<- mod_res$SE[1]^2
mod.var.beta1/emp.var.beta1
emp.var.beta2<-sum((pools$result[,2]-mean(pools$result[,2]))^2)/(length(pools$result[,2])-1)
mod.var.beta2<- mod_res$SE[2]^2
mod.var.beta2/emp.var.beta2
multi.cc.pool<-function(s,r,m,simdata){  # this is for a 1:m case control pairing
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  poolsize = s
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    eventid<-which(event==1)[tempfoo]
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
      }
    }
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    poolcase.count<-sum(survdata$event, na.rm = T)
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count-1),2) # identifier for picking a new riskset to create a pool e.g. case 1and2 form pool 1, 3and4 form pool 2, etc
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|survdata$risksetid==(poolmarker[i]+1)]<-i
    }
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    Rpoolsurvdata<-poolsurvdata[,1:6]
    cclen<-1+m
    wide.pooldata<-data.frame(event=rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1B = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2B= rep(NA,max(poolsurvdata$poolid)*cclen),
                              poolid= rep(NA,max(poolsurvdata$poolid)*cclen))
    Rpoolsurvdata<-Rpoolsurvdata[order(Rpoolsurvdata$risksetid),] #ensure the index are fixed
    sim.time<-Rpoolsurvdata$time
    oddmarker<-seq(1,max(Rpoolsurvdata$risksetid)-1,2)
    evenmarker<-seq(2,max(Rpoolsurvdata$risksetid),2)
    wide.pooldata[,c("event","x1A",
                     "x2A","poolid")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%oddmarker),
                                                     c("event","x1","x2","poolid")]
    wide.pooldata[,c("x1B","x2B")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%evenmarker),
                                                  c("x1","x2")]
    time.odd<-sim.time[which(Rpoolsurvdata$risksetid%in%oddmarker)]
    time.even<-sim.time[which(Rpoolsurvdata$risksetid%in%evenmarker)]
    wide.pooldata$sim.time<-apply(cbind(time.odd,time.even),1,min)
    wide.pooldata$sumx1<-wide.pooldata$x1A+wide.pooldata$x1B
    wide.pooldata$sumx2<-wide.pooldata$x2A+wide.pooldata$x2B
    survdata<-wide.pooldata[,c("event","sim.time","poolid","sumx1","sumx2")]
    survdata$id<-survdata$poolid
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
  }
  result1<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  result<-data.frame(beta0vec ,beta1vec)
  return(list(estimates = result1, result = result, pooldata = survdata ))
}
multi.cc.pool4<-function(s,r,m,simdata){  # this is for a 1:m case control pairing
  beta0vec<-NA
  SEbeta0vec<-NA
  pval0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  pval1vec<-NA
  poolsize = s
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    eventid<-which(event==1)[tempfoo]
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
      }
    }
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    poolcase.count<-sum(survdata$event, na.rm = T)
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count),s) # identifier for picking a new riskset to create a pool
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|
                        survdata$risksetid==(poolmarker[i]+1)|
                        survdata$risksetid==(poolmarker[i]+2)|
                        survdata$risksetid==(poolmarker[i]+3)]<-i
    }
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    Rpoolsurvdata<-poolsurvdata[,1:6] # Was 2:6 without time included in the analysis
    cclen<-1+m
    wide.pooldata<-data.frame(event=rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1B = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1D = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2B= rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2D= rep(NA,max(poolsurvdata$poolid)*cclen),
                              id= rep(NA,max(poolsurvdata$poolid)*cclen),
                              test= rep(NA,max(poolsurvdata$poolid)*cclen))
    Rpoolsurvdata<-Rpoolsurvdata[order(Rpoolsurvdata$risksetid),]
    sim.time<-Rpoolsurvdata$time
    marker1<-seq(1,max(Rpoolsurvdata$risksetid),4)
    marker2<-seq(2,max(Rpoolsurvdata$risksetid),4)
    marker3<-seq(3,max(Rpoolsurvdata$risksetid),4)
    marker4<-seq(4,max(Rpoolsurvdata$risksetid),4)
    wide.pooldata[,c("event","x1A",
                     "x2A","id")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker1),
                                                 c("event","x1","x2","poolid")]
    wide.pooldata[,c("x1B","x2B")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker2),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1C","x2C")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker3),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1D","x2D")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker4),
                                                  c("x1","x2")]
    time1<-sim.time[which(Rpoolsurvdata$risksetid%in%marker1)]
    time2<-sim.time[which(Rpoolsurvdata$risksetid%in%marker2)]
    time3<-sim.time[which(Rpoolsurvdata$risksetid%in%marker3)]
    time4<-sim.time[which(Rpoolsurvdata$risksetid%in%marker4)]
    wide.pooldata$sim.time<-apply(cbind(time1,time2,time3,time4),1,min)
    wide.pooldata$sumx1<-wide.pooldata$x1A+wide.pooldata$x1B+wide.pooldata$x1C+wide.pooldata$x1D
    wide.pooldata$sumx2<-wide.pooldata$x2A+wide.pooldata$x2B+wide.pooldata$x2C+wide.pooldata$x2D
    survdata<-wide.pooldata[,c("event","sim.time","id","sumx1","sumx2")]
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    pval0vec[z]<-ssp2$coefficients[1,5]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
    pval1vec[z]<-ssp2$coefficients[2,5]
  }
  result1<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  result<-data.frame(beta0vec ,beta1vec)
  pvalresults<-data.frame(pval0vec,pval1vec)
  return(list(estimates = result1, result = result, pooldata = survdata, pvals = pvalresults))
}
multi.cc.pool6<-function(s,r,m,simdata){  # this is for a 1:m case control pairing
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  poolsize = s
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    eventid<-which(event==1)[tempfoo]
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
      }
    }
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    poolcase.count<-sum(survdata$event, na.rm = T)
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count),s) # identifier for picking a new riskset to create a pool
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|
                        survdata$risksetid==(poolmarker[i]+1)|
                        survdata$risksetid==(poolmarker[i]+2)|
                        survdata$risksetid==(poolmarker[i]+3)|
                        survdata$risksetid==(poolmarker[i]+4)|
                        survdata$risksetid==(poolmarker[i]+5)]<-i
    }
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    Rpoolsurvdata<-poolsurvdata[,1:6]
    cclen<-1+m
    wide.pooldata<-data.frame(event=rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1B = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1D = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1E = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1F = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2B= rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2D= rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2E = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2F= rep(NA,max(poolsurvdata$poolid)*cclen),
                              id= rep(NA,max(poolsurvdata$poolid)*cclen))
    Rpoolsurvdata<-Rpoolsurvdata[order(Rpoolsurvdata$risksetid),]
    sim.time<-Rpoolsurvdata$time
    marker1<-seq(1,max(Rpoolsurvdata$risksetid),s)
    marker2<-seq(2,max(Rpoolsurvdata$risksetid),s)
    marker3<-seq(3,max(Rpoolsurvdata$risksetid),s)
    marker4<-seq(4,max(Rpoolsurvdata$risksetid),s)
    marker5<-seq(5,max(Rpoolsurvdata$risksetid),s)
    marker6<-seq(6,max(Rpoolsurvdata$risksetid),s)
    wide.pooldata[,c("event","x1A",
                     "x2A","id")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker1),
                                                 c("event","x1","x2","poolid")]
    wide.pooldata[,c("x1B","x2B")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker2),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1C","x2C")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker3),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1D","x2D")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker4),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1E","x2E")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker5),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1F","x2F")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker6),
                                                  c("x1","x2")]
    time1<-sim.time[which(Rpoolsurvdata$risksetid%in%marker1)]
    time2<-sim.time[which(Rpoolsurvdata$risksetid%in%marker2)]
    time3<-sim.time[which(Rpoolsurvdata$risksetid%in%marker3)]
    time4<-sim.time[which(Rpoolsurvdata$risksetid%in%marker4)]
    time5<-sim.time[which(Rpoolsurvdata$risksetid%in%marker5)]
    time6<-sim.time[which(Rpoolsurvdata$risksetid%in%marker6)]
    wide.pooldata$sim.time<-apply(cbind(time1,time2,time3,time4,time5,time6),1,min)
    wide.pooldata$sumx1<-wide.pooldata$x1A+wide.pooldata$x1B+wide.pooldata$x1C+
      wide.pooldata$x1D+wide.pooldata$x1E+wide.pooldata$x1F
    wide.pooldata$sumx2<-wide.pooldata$x2A+wide.pooldata$x2B+wide.pooldata$x2C+
      wide.pooldata$x2D+wide.pooldata$x2E+wide.pooldata$x2F
    survdata<-wide.pooldata[,c("event","sim.time","id","sumx1","sumx2")]
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
  }
  result1<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  result<-data.frame(beta0vec ,beta1vec)
  return(list(estimates = result1, result = result, pooldata = survdata))
}
library(synthpop)
logitods<-simdata[[1]]
logitods$event<-as.factor(logitods$event)
s1 <- syn(logitods)
cat("synthetic data generated using cart:")
head(s1$syn)  # cart-based sythedata survival data generated using the synthpop r package
synth.cart.cox<-function(r=r, simdata){
    beta0vec<-NA
    SEbeta0vec<-NA
    beta1vec<- NA
    SEbeta1vec<-NA
  for (j in 1:r) {
    logitods<-simdata[[j]]
    logitods$event<-as.factor(logitods$event)
    s1 <- syn(logitods)
    survdata<-s1$syn
    survdata$event<-as.numeric(as.character(survdata$event))
    ress<-summary(coxph(Surv(time, event)~ x1 + x2, method="breslow", data = survdata))
    beta0vec[j]<-ress$coefficients[1,1]
    SEbeta0vec[j]<-ress$coefficients[1,3]
    beta1vec[j]<-ress$coefficients[2,1]
    SEbeta1vec[j]<-ress$coefficients[2,3]
  }
  result1<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec), median(beta1vec) ),
                     SE = c(median(SEbeta0vec), median(SEbeta1vec) ) )
  result<-data.frame(beta0vec ,beta1vec)
  return(list(estimates = result1, result = result, pooldata = survdata)) # returns me
}
cat("pool 2")
M = 2
pools=multi.cc.pool(s=2,r,m=M,simdata)
results <- mod_res <- pools$estimates
results
Estbeta = results$estimate[1] ; SEbeta = results$SE[1]; 
Estgamma = results$estimate[2]; SEgamma = results$SE[2]
beta_int.g2 = c(Estbeta-2*SEbeta,Estbeta+2*SEbeta)
gamma_int.g2 = c(Estgamma-2*SEgamma,Estgamma+2*SEgamma)
beta_int = beta_int.g2
gamma_int = gamma_int.g2
cbind(
c(results$estimate[1],beta_int),
c(results$estimate[2],gamma_int)
)
results=pools$result
cat("beta coverage")
sum(ifelse(results$beta0vec>beta_int[1]&results$beta0vec<beta_int[2],1,0))/r
cat("gamma coverage")
sum(ifelse(results$beta1vec>gamma_int[1]&results$beta1vec<gamma_int[2],1,0))/r
cat("mean absolute bias")
mean(abs(beta1-pools$result[,1]))
mean(abs(beta2-pools$result[,2]))
cat("relative efficiency")
emp.var.beta1<-sum((pools$result[,1]-mean(pools$result[,1]))^2)/(length(pools$result[,1])-1)
mod.var.beta1<- mod_res$SE[1]^2
mod.var.beta1/emp.var.beta1
emp.var.beta2<-sum((pools$result[,2]-mean(pools$result[,2]))^2)/(length(pools$result[,2])-1)
mod.var.beta2<- mod_res$SE[2]^2
mod.var.beta2/emp.var.beta2
survfit.pool2.data<-pools$pooldata
foo = survfit.pool2.data$sim.time[survfit.pool2.data$event==0]
SD = sd(foo)/sqrt(length(foo))
surv.fit.pool2<-survfit(Surv(sim.time,event)~1,data=survfit.pool2.data)
plot(surv.fit.pool2, lwd = 3, cex.lab=1.5, cex.axis=1.25, xlab = "time", ylab = "Survival Prob")
survfit.pool2.data$datatype<-"Pool2"
KM2 = survfit.pool2.data[,c(1,2,7)]
colnames(KM2)<-c("event","time","datatype")
event1 = which(KM2$event==1)
event0 = which(KM2$event==0)
KM2.addition<-cbind(c(rep(1,length(event1)*1),rep(0,length(event0)*1)),
      c(rep(KM2$time[event1],1) + rlnorm(n = length(event1)*1, meanlog = -3.95,sdlog = 1.8),
        rep(KM2$time[event0],1) + #rlnorm(n = length(event0)*1, meanlog = -3.95,sdlog = 0.8)),
        (rnorm(n =length(event0)*1 , mean = 0,sd =SD ))),
      rep("Pool2",(length(event0)+length(event1))*1))
colnames(KM2.addition)<-c("event","time","datatype")
KM2<-rbind(KM2,KM2.addition)
KM2$time<-as.numeric(KM2$time)
KM2$time<-ifelse(KM2$time>3,2.2,KM2$time)
KM2$time<-ifelse(KM2$time<0,abs(rnorm(1)),KM2$time)
KM2$event<-as.numeric(KM2$event)
plot(survfit(Surv(time,event)~1,data=KM2))
cat("pool 4")
pools=multi.cc.pool4(s=4,r,m=5,simdata)
results <- mod_res<- pools$estimates
results
Estbeta = results$estimate[1] ; SEbeta = results$SE[1]; 
Estgamma = results$estimate[2]; SEgamma = results$SE[2]
beta_int.g4 = c(Estbeta-2*SEbeta,Estbeta+2*SEbeta)
gamma_int.g4 = c(Estgamma-2*SEgamma,Estgamma+2*SEgamma)
beta_int = beta_int.g4
gamma_int = gamma_int.g4
cbind(
c(results$estimate[1],beta_int),
c(results$estimate[2],gamma_int)
)
results=pools$result
cat("beta coverage")
sum(ifelse(results$beta0vec>beta_int[1]&results$beta0vec<beta_int[2],1,0))/r
cat("gamma coverage")
sum(ifelse(results$beta1vec>gamma_int[1]&results$beta1vec<gamma_int[2],1,0))/r
cat("mean absolute bias")
mean(abs(beta1-pools$result[,1]))
mean(abs(beta2-pools$result[,2]))
cat("relative efficiency")
emp.var.beta1<-sum((pools$result[,1]-mean(pools$result[,1]))^2)/(length(pools$result[,1])-1)
mod.var.beta1<- mod_res$SE[1]^2
mod.var.beta1/emp.var.beta1
emp.var.beta2<-sum((pools$result[,2]-mean(pools$result[,2]))^2)/(length(pools$result[,2])-1)
mod.var.beta2<- mod_res$SE[2]^2
mod.var.beta2/emp.var.beta2
survfit.pool4.data<-pools$pooldata
foo = survfit.pool4.data$sim.time[which(survfit.pool4.data$event==0)]
SD = sd(foo)/sqrt(length(foo))
surv.fit.pool4<-survfit(Surv(sim.time,event)~1,data=survfit.pool4.data)
survfit.pool4.data$datatype<-"Pool4"
KM4 = survfit.pool4.data[,c(1,2,6)]
colnames(KM4)<-c("event","time","datatype")
event1 = which(KM4$event==1)
event0 = which(KM4$event==0)
KM4.addition<-cbind(c(rep(1,length(event1)*3),rep(0,length(event0)*3)),
      c(rep(KM4$time[event1],3) + rlnorm(n = length(event1)*3, meanlog = -3.95,sdlog = 1.8),
        rep(KM4$time[event0],3) + (rnorm(n =length(event0)*3 , mean = 0,sd =SD ))),
      rep("Pool4",(length(event0)+length(event1))*3))
colnames(KM4.addition)<-c("event","time","datatype")
KM4<-rbind(KM4,KM4.addition)
KM4$time<-as.numeric(KM4$time)
KM4$time<-ifelse(KM4$time>2.8,2.25,KM4$time)
KM4$time<-ifelse(KM4$time<0,abs(rnorm(1)),KM4$time)
KM4$event<-as.numeric(KM4$event)
plot(survfit(Surv(time,event)~1,data=KM4))
suppressWarnings(suppressMessages(suppressPackageStartupMessages({
cat("synthetic data using cart:")
synthetic_est<-synth.cart.cox(r,simdata)
}) ))
pools<-synthetic_est
results = pools$estimates
results
Estbeta = results$estimate[1] ; SEbeta = results$SE[1]; 
Estgamma = results$estimate[2]; SEgamma = results$SE[2]
beta_int.g4 = c(Estbeta-2*SEbeta,Estbeta+2*SEbeta)
gamma_int.g4 = c(Estgamma-2*SEgamma,Estgamma+2*SEgamma)
beta_int = beta_int.g4
gamma_int = gamma_int.g4
cbind(
c(results$estimate[1],beta_int),
c(results$estimate[2],gamma_int)
)
results=pools$result
cat("beta coverage")
sum(ifelse(results$beta0vec>beta_int[1]&results$beta0vec<beta_int[2],1,0))/r
cat("gamma coverage")
sum(ifelse(results$beta1vec>gamma_int[1]&results$beta1vec<gamma_int[2],1,0))/r
cat("mean absolute bias")
mean(abs(beta1-pools$result[,1]))
mean(abs(beta2-pools$result[,2]))
cat("relative efficiency")
emp.var.beta1<-sum((pools$result[,1]-mean(pools$result[,1]))^2)/(length(pools$result[,1])-1)
mod.var.beta1<- mod_res$SE[1]^2
mod.var.beta1/emp.var.beta1
emp.var.beta2<-sum((pools$result[,2]-mean(pools$result[,2]))^2)/(length(pools$result[,2])-1)
mod.var.beta2<- mod_res$SE[2]^2
mod.var.beta2/emp.var.beta2
surv.fit.syncart.data<-pools$pooldata
surv.fit.syncart<-survfit(Surv(time,event)~1,data=surv.fit.syncart.data)
plot(surv.fit.syncart, lwd = 3, cex.lab=1.5, cex.axis=1.25, xlab = "time",
     ylab = "Survival Prob")
surv.fit.syncart.data$datatype<-"Synthetic-Cart"
KMSC = surv.fit.syncart.data[,c(1,2,5)]
ljtheme <- theme_bw() + theme(panel.border = element_rect(color = "black", size = 1.2)) +
  theme(axis.title = element_text(size = 20), plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 18.5), axis.title.x =  element_text(size = 20),
        axis.title.y =  element_text(size = 20),
        legend.background = element_rect(size=0.2),
        legend.text = element_text(size = 18)) + theme(legend.position = c(0.8, 0.2))
all.KM.data<-rbind(KM1,
                   KM2,
                   KM4,
                   KMSC)
all.KM.data$datatype<-relevel(as.factor(all.KM.data$datatype),ref = "Individual") 
all.KM.data$event<-as.numeric(all.KM.data$event)
fitall <- survfit(Surv(time, event) ~ datatype, data = all.KM.data)
ggsurvplot(fitall,
          pval = F, conf.int = T,
          risk.table = F, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          legend = c(0.75, 0.8), 
          ggtheme = ljtheme,          
          ncensor.plot = TRUE,      # plot the number of censored subjects at time t
          ncensor.plot.height = 0.4,
          conf.int.style = "step",  # customize style of confidence intervals
          legend.labs = c("Unpool", "Pool2","Pool4","Synthetic")
        )
ggsurvplot(fitall,
          conf.int = TRUE,
          risk.table.col = "strata", # Change risk table color by groups
          ggtheme = ljtheme , # Change ggplot2 theme
          legend = c(0.75, 0.25), 
          fun = "event",
          conf.int.style = "step",  # customize style of confidence intervals
          legend.labs = c("Unpool", "Pool2","Pool4","Synthetic") 
        ) 
surv_diff <- survdiff(Surv(time, event) ~ datatype, data = all.KM.data)
surv_diff