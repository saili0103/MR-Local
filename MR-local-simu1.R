
library(mr.raps)
library(RobustIV)
library(MendelianRandomization)
Niter = 500

p <- 2000
n<-10^5
sd.gamD<-runif(p,0.8,1)/sqrt(n)
sd.gamY<-runif(p,0.8,1)/sqrt(n)
mean(gam.D^2 / sd.gamD^2) #IV strength
plot(density(gam.D/sd.gamD))
source('funQ2-IVW.r') #the main functions
tau0.list=seq(1.2,1.8,length.out=4) #difference choice of tau_0
e0=1
set.seed(e0)
if(e0==1){#majority rule
  gam.D<-rnorm(p,0, sqrt(0.1/p))
  pi.Y<-rnorm(p, 0, sqrt(0.05/p)) + 2.5*gam.D
  Vdy<- sample(1:length(gam.D), floor(0.5*length(gam.D)), replace = F)
  mean(gam.D^2 / sd.gamD^2) #IV strength
  pi.Y[Vdy]<-0
}else if(e0==2){#balanced pleiotropy
  gam.D<-rnorm(p,0, sqrt(0.1/p))
    pi.Y<-rnorm(p, 0, sqrt(0.1/p)) 
    Vdy<-1:length(gam.D)
}else if(e0==3){
  gam.D<-rnorm(p,0, sqrt(0.1/p))
  pi.Y<-rnorm(p,0, sqrt(0.05/p)) + 3*gam.D
  Vdy<- sample(1:length(gam.D), floor(0.28*length(gam.D)), replace = F)
  pi.Y[Vdy]<-0
}else if(e0==4){ #plurality rule
  gam.D<-rnorm(p,0, sqrt(0.1/p))
  pi.Y<-rnorm(p, 0, sqrt(0.5/p)) + 2.5*gam.D
  Vdy<- sample(1:length(gam.D), floor(0.5*length(gam.D)), replace = F)
  pi.Y[Vdy]<-0
}else if(e0==5){ #balaced pleiotropy
  gam.D<-rnorm(p,0, sqrt(0.1/p))
  pi.Y<-rnorm(p, 0, sqrt(0.05/p)) 
  Vdy<-1:length(gam.D)
}
re.all<-NULL
  for(betady in c(0,0.1)){
    gam.Y <- gam.D * betady + pi.Y
    for (tau0 in tau0.list) {
      prop.summ0 <- matrix(NA, ncol = 5, nrow = Niter)
      prop.summ <- matrix(NA, ncol = 5, nrow = Niter)
      tsht.summ <- matrix(NA, ncol = 5, nrow = Niter)
      mode.summ <- matrix(NA, ncol = 5, nrow = Niter)
      raps.summ <- matrix(NA, ncol = 5, nrow = Niter)
      for (it in 1:Niter) {
        #data generation
        gam.D.hat <-sapply(1:p, function(k) rnorm(1, gam.D[k], sd.gamD[k]))
        gam.Y.hat <-sapply(1:p, function(k) rnorm(1, gam.Y[k], sd.gamY[k]))
        #plot(density((gam.Y.hat/gam.D.hat)[abs(gam.D.hat/sd.gamD)>=tau0]))
        prop.re0 <- MR.local(beta.out = gam.Y.hat, beta.exp = gam.D.hat, tau.s=tau0,    #MR-Local
                            se.out = sd.gamY, se.exp = sd.gamD, tau0 = tau0, skew.ind = F)
        cov.prop0 <- (prop.re0$divw.ci[1] <= betady) & (prop.re0$divw.ci[2] >= betady)
        prop.summ0[it, ] <- c(abs(prop.re0$beta.divw - betady), cov.prop0,prop.re0$sd.divw,
                             length(intersect(prop.re0$clust.sel, Vdy)) / length(prop.re0$clust.sel),
                             prop.re0$balanced.pleiotropy)
        #cat(it,prop.summ0[it,],'\n')
        prop.re <- MR.local(beta.out = gam.Y.hat, beta.exp = gam.D.hat,     #MR-Local+
                    se.out = sd.gamY, se.exp = sd.gamD, tau0 = tau0, tau.s=tau0)
        cov.prop <- (prop.re$divw.ci[1] <= betady) & (prop.re$divw.ci[2] >= betady)
        prop.summ[it, ] <- c(abs(prop.re$beta.divw - betady), cov.prop,prop.re$sd.divw,
                             length(intersect(prop.re$clust.sel, Vdy)) / length(prop.re$clust.sel),
                             prop.re$balanced.pleiotropy)
        cat(it,prop.summ[it,],'\n')
        if (tau0 > tau0.list[1]) {next}
        #MR-raps
        raps.re<-mr.raps(b_exp=gam.D.hat,b_out=gam.Y.hat,se_exp=sd.gamD, se_out=sd.gamY,
                         over.dispersion=T)
        cov.raps<-(abs(raps.re$beta.hat - betady)<=qnorm(0.975)*raps.re$beta.se)
        raps.summ[it, ] <- c(abs(raps.re$beta.hat - betady), cov.raps,raps.re$beta.se, 
                             length(Vdy)/length(gam.D.hat),length(gam.D.hat))
        #cat(it,raps.summ[it,],'\n')
        #TSHT
        Sd.hat<- which(abs(gam.D.hat)>sqrt(2*log(p))*sd.gamD)
        s=length(Sd.hat)
        re.tsht <-RobustIV:::TSHT.VHat(n=s, ITT_Y=gam.Y.hat[Sd.hat], ITT_D=gam.D.hat[Sd.hat], V.Gamma=diag((sd.gamY[Sd.hat])^2*s),
                            V.gamma=diag((sd.gamD[Sd.hat])^2*s), C=0, voting ='MP') #MR-Median
        tsht.ivw <- ivw.fun(gam.D.hat[Sd.hat[re.tsht$VHat]], gam.Y.hat[Sd.hat[re.tsht$VHat]],
                        (sd.gamD[Sd.hat[re.tsht$VHat]])^2,(sd.gamY[Sd.hat[re.tsht$VHat]])^2)
        cov.tsht <- (tsht.ivw$ci.ivw[1] <= betady) & (tsht.ivw$ci.ivw[2] >= betady)
        tsht.summ[it, ] <-
          c(abs(tsht.ivw$beta.ivw - betady), cov.tsht, tsht.ivw$sd.ivw, 
            length(intersect(Vdy, Sd.hat[re.tsht$VHat]))/length(re.tsht$VHat),length(re.tsht$VHat))
    
        # #MR-MBE
        re.mode <-MendelianRandomization::mr_mbe(mr_input(
          bx = gam.D.hat[Sd.hat], bxse = sd.gamD[Sd.hat],
          by = gam.Y.hat[Sd.hat], byse = sd.gamY[Sd.hat]), iterations = 100, weighting='unweighted')
        cov.mode <- (re.mode@CILower <= betady) & (re.mode@CIUpper >= betady)
        mode.summ[it, ] <- c(abs(re.mode@Estimate - betady),cov.mode, re.mode@StdError,
                             length(intersect(Vdy, Sd.hat)) / length(Sd.hat), length(Sd.hat))
        #cat(it,mode.summ[it,],'\n')
        
      }

      if (tau0==tau0.list[1]) {
        result.mat <- rbind(
          colMeans(prop.summ0),
          colMeans(prop.summ),
          colMeans(tsht.summ),
          colMeans(mode.summ),
          colMeans(raps.summ))
        result.df <-
          data.frame(e=rep(e0, 5), betady = rep(betady, 5), tau0 = rep(tau0, 5), 
                     method = c('prop0','prop','TSHT', 'mode','raps'))
        result.df <- cbind(result.df, result.mat)
      }else{
        result.mat <- rbind(
          colMeans(prop.summ0),
          colMeans(prop.summ))
        result.df <- data.frame(e=rep(e0,2), betady = rep(betady,2), tau0 = rep(tau0, 2), method = c('prop0','prop'))
        result.df <- cbind(result.df, result.mat)
      }
      re.all<-rbind(re.all, result.df)
    }
  }


write.table(re.all, file=paste('MRlocal-numerical-setting',e0,'.txt',sep=""), row.names=F)

