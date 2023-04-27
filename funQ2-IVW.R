#beta.out = gam.Y.hat;beta.exp = gam.D.hat;se.out = sd.gamY;se.exp = sd.gamD
library(MASS)
library(MendelianRandomization) #use mr_divw 
MR.local<-function(beta.out, beta.exp, se.out, se.exp, tau0=NULL,alpha=0.95, skew.ind=T, tau.s=0){
  C.beta=2
  p<-length(beta.out)
  if(is.null(tau0)){ #default value for tau_0
    tau0=qnorm(1-1/p)
  }
  ivw.re<-VA.IVW.est(beta.out=beta.out, beta.exp=beta.exp, 
                     se.out=se.out, se.exp=se.exp, tau0=tau0,, tau.s=tau.s,skew.ind=skew.ind)
  clust.sel=ivw.re$clust.sel
  sd.divw0<-ivw.re$sd.divw
  if(ivw.re$bp==T){ #sd for the balanced pleiotropy case
    sd.divw=sd.divw0
  }else{ #bootstrap
    bs.vec<-NULL
    for(i in 1:50){
      gamD.bs<-rnorm(p, beta.exp, se.exp)
      gamY.bs<-rnorm(p, beta.out, se.out)
      x.grid<-seq(-C.beta, C.beta,length.out=max(C.beta*200,p)+1)
      s=length(x.grid)
      B.bs<-which(sapply(1:s, function(k) min(abs(x.grid[k]-ivw.re$grid.B))<=2/s))
      clust.size<-rep(0,length(B.bs))
      clust.bs<-list()
      for(j in 1:length(B.bs)){
        x.cur<-x.grid[B.bs[j]]
        clust.bs[[j]]<- which(abs(gamY.bs/gamD.bs-x.cur) <=tau0*sqrt(se.out^2/gamD.bs^2+min(x.cur^2,C.beta^2)*se.exp^2/gamD.bs^2)) 
        clust.size[j]<-length(clust.bs[[j]])
      }
      max.size=max(clust.size)
      j.bs<-which(clust.size==max.size)
      clust.bs=unique(as.vector(unlist(sapply(j.bs, function(k) clust.bs[[k]]))))
      divw.bs<- mr_divw(mr_input(bx=gamD.bs[clust.bs],by=gamY.bs[clust.bs], bxse=se.exp[clust.bs], byse=se.out[clust.bs]))
      bs.vec<-c(bs.vec, divw.bs@Estimate)
    }
   sd.divw=sqrt(mean((bs.vec-ivw.re$beta.divw)^2))
   # sd.divw=sd(bs.vec)
  }

  divw.ci=c(ivw.re$beta.divw-qnorm(0.975)*sd.divw, ivw.re$beta.divw+qnorm(0.975)*sd.divw)

  return(list(beta.divw=ivw.re$beta.divw, divw.ci=divw.ci, 
              sd.divw0=sd.divw0, sd.divw=sd.divw, 
              clust.sel=clust.sel, balanced.pleiotropy=ivw.re$bp))
}
sd.skewness<-function(tau0){
  z<-rnorm(5000)
  sub1<-which(abs(z)<=tau0)
  z.truc<-z[sub1]
  sd(z.truc^3)
}

VA.IVW.est<-function(beta.out, beta.exp, se.out, se.exp, tau0,tau.s, C.beta=2, diagnostics=F, skew.ind){
  p<-length(beta.out)
  Sd.hat<-which(abs(beta.exp/se.exp)>=tau.s)
  sd.YS<-se.out[Sd.hat]; sd.DS<-se.exp[Sd.hat]
  gam.DS<-beta.exp[Sd.hat]; gam.YS<-beta.out[Sd.hat]
  clust.hat=list() #list for C(b)
  rho=list() #list for rho
  betaS<-gam.YS/gam.DS
  #plot(density(betaS))
  #the grid values of b
  x.grid<-seq(-C.beta, C.beta,length.out=max(C.beta*200,length(betaS))+1)
  s=length(x.grid)
  for(j in 1:s){
    clust.temp<- which(abs(betaS-x.grid[j]) <=tau0*sqrt(sd.YS^2/gam.DS^2+min(x.grid[j]^2,C.beta^2)*sd.DS^2/gam.DS^2)) 
    rho[[j]]<-sqrt(sd.YS^2/gam.DS^2+x.grid[j]^2*sd.DS^2/gam.DS^2)[clust.temp]
    clust.hat[[j]]<-clust.temp
  }
  #uncertainty measure
  r1<-1-2*dnorm(tau0)*tau0/(pnorm(tau0)-pnorm(-tau0))
  r2<-(3-2*dnorm(tau0)*tau0*(3+tau0^2)/(pnorm(tau0)-pnorm(-tau0)))-r1^2 #variance of Qj
  beta.summ<-matrix(0,nrow=s,ncol=4)#b_j, Q.hat, size of cluster
  for(j in 1:s){
      if(length(clust.hat[[j]])<=sqrt(s) |length(clust.hat[[j]])<=s/log(s)){ #too small cliques
        next}
      var.beta.j<-mean((betaS[clust.hat[[j]]]-x.grid[j])^2/rho[[j]]^2)    #Q-statistic
      skew.j=mean((betaS[clust.hat[[j]]]-x.grid[j])^3/rho[[j]]^3)*sqrt(length(clust.hat[[j]])) # a skewness measure
      beta.summ[j,]<-c(x.grid[j], var.beta.j-r1, length(clust.hat[[j]]),skew.j)
  }
  B.hat<-which(beta.summ[,2]<=sqrt(log(s))*sqrt(r2/beta.summ[,3]) & beta.summ[,3]>=s/log(s))#threshold for Qj
  if(skew.ind){
    sd.sk<-sd.skewness(tau0)
    B.hat<-intersect(B.hat, which(abs(beta.summ[,4])/sd.sk<=sqrt(log(s))))
  }
  if(length(B.hat)==0){
    clust.sel=1:p
    grid.B=0
  #  divw.est<-mr_divw(mr_input(bx=beta.exp,by=beta.out, bxse=se.exp, byse=se.out))
  #  Sd.hat<- which(abs(beta.exp/se.exp)>sqrt(2*log(p)))
    bp=T#balanced pleiotropy
  }else{
    max.size=max(beta.summ[B.hat,3])
    j.hat<-B.hat[beta.summ[B.hat,3]==max.size]
    clust.sel=Sd.hat[unique(as.vector(unlist(sapply(j.hat, function(k) clust.hat[[k]]))))]
    grid.B=x.grid[B.hat]
    bp=F
  }
  #debiased IVW
  divw.est<-mr_divw(mr_input(bx=beta.exp[clust.sel],by=beta.out[clust.sel], bxse=se.exp[clust.sel], byse=se.out[clust.sel]))
  if(diagnostics){#diagnostics
    h=0.9*min(sd(betaS[abs(betaS)<=C.beta]), mad(betaS[abs(betaS)<=C.beta]))/length(betaS)^(1/5)
    fhat=density(betaS,bw=h, from=-C.beta, to=C.beta, n=max(C.beta*200,length(betaS))+1)
    plot(fhat)
    points(fhat$x[B.hat],fhat$y[B.hat],col='red')
    beta.summ[j.hat,]
    #j.star<-which.min(abs(x.grid-betady))
    #beta.summ[j.star,]
  }

  return(list(beta.divw=divw.est@Estimate, sd.divw=divw.est@StdError, grid.B=grid.B,
              clust.sel=clust.sel, bp=bp))
}


ivw.fun<-function(gam.D,gam.Y, sig.gamD,sig.gamY){
  beta.vec=gam.Y/gam.D
  beta.ivw=sum(beta.vec*gam.D^2/sig.gamY)/sum(gam.D^2/sig.gamY)
  rho0.sq<-(sig.gamY/gam.D^2)
  sd.ivw=1/sqrt(sum(1/rho0.sq))
  ci.ivw<-c(beta.ivw-qnorm(0.975)*sd.ivw, beta.ivw+qnorm(0.975)*sd.ivw)
  return(list(beta.ivw=beta.ivw, sd.ivw=sd.ivw, ci.ivw=ci.ivw))
}

divw.fun<-function(gam.D, gam.Y, sig.gamD,sig.gamY){
  betaS=gam.Y/gam.D
  beta.divw <- sum(betaS*(gam.D^2/sig.gamY))/sum((gam.D^2/sig.gamY-sig.gamD/sig.gamY))
  tausq<-sum((gam.Y-beta.divw*gam.D)^2/sig.gamY-1-  beta.divw^2*sig.gamD/sig.gamY)/sum(1/sig.gamY)
  tausq=max(tausq,0)
  var.num<- sum(gam.D^2/sig.gamY*(1+tausq/sig.gamY)+beta.divw^2*sig.gamD*(gam.D^2+sig.gamD)/sig.gamY^2)
  sd.divw <- sqrt(var.num)/sum((gam.D^2/sig.gamY-sig.gamD/sig.gamY))
  return(list(beta.divw=beta.divw, sd.divw=sd.divw))
}
