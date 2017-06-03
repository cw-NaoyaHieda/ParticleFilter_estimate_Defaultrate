old_d<-getwd()
setwd("C:/Users/naoya/Desktop/Graduage_study/2017ParticleFilter_estimate_Defaultrate/script")
library(rgl)
library(mvtnorm)
library(reshape2)
library(ggplot2)
library(doSNOW)
library(foreach)
library(copula)
source('DR_density.R', encoding = 'UTF-8')



AR_sim<-function(time=100,rho=0.08,PD=0.035){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  mu<-c(sig_env(rho),sig_env(PD))
  tau<-c(0.99,0.99)
  sigma<-matrix(c(0.0005,0,0,0.0004),ncol=2)
  x_0<-matrix(c(sig_env(rho),sig_env(PD)),ncol=2)
  x<-x_0
  
  for(i in 1:(time-1)){
    tmp_x<-mu+tau*(x[i,]-mu)+rmvnorm(1,c(0,0),sigma)
    x<-rbind(x,tmp_x)
  }
  x<-sig(x)
  x<-data.frame(x)
  colnames(x)<-c("rho","PD")
  
  
  DR <- {}
  for(i in 1:time){
    #SIR_DRはDR_density.Rの関数
    DR<-c(DR,SIR_DR(L=1,rho=x[i,1],pd=x[i,2])$q)
  }
  
  
  
  print(head(x))
  x <- cbind(x,DR)
  plot(x$rho,type="l",ylab=expression(rho))
  plot(x$PD,type="l",ylab=expression(PD))
  plot(x$DR,type="l",ylab=expression(DR))
  
  out<-data.frame(x)
}

AR_sim2<-function(time=100,rho=0.08,PD=0.035){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  mu<-c(sig_env(rho))
  tau<-c(0.99)
  sigma<-0.05
  x_0<-sig_env(rho)
  x<-x_0
  
  for(i in 1:(time-1)){
    tmp_x<-mu+tau*(x[i]-mu)+rnorm(1,0,sigma)
    x<-rbind(x,tmp_x)
  }
  
  x<-sig(x)
  row.names(x)<-NULL
  data.frame(rho=x)
  
  print(head(x))
  out<-data.frame(x)
}

AR_sim3<-function(time=100,rho=0.08,PD=0.035){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  mu<-c(sig_env(PD))
  tau<-c(0.99)
  sigma<-0.05
  x_0<-sig_env(PD)
  x<-x_0
  
  for(i in 1:(time-1)){
    tmp_x<-mu+tau*(x[i]-mu)+rnorm(1,0,sigma)
    x<-rbind(x,tmp_x)
  }
  
  x<-sig(x)
  row.names(x)<-NULL
  data.frame(rho=x)
  
  print(head(x))
  
  out<-data.frame(x)
}

local_trend_sim<-function(time=100,rho=0.08,PD=0.035,sigma_1=0.001,sigma_2=0.001,sigma_3=0.0005){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  
  beta_0<-c(0,0)
  mu_0<-c(sig_env(rho),sig_env(PD))
  sigma_1<-diag(c(sigma_1,sigma_1))
  sigma_2<-diag(c(sigma_2,sigma_2))
  sigma_3<-diag(c(sigma_3,sigma_3))
  beta<-beta_0+rmvnorm(1,c(0,0),sigma_3^2)
  mu<-beta_0+mu_0+rmvnorm(1,c(0,0),sigma_2^2)
  x<-mu+rmvnorm(1,c(0,0),sigma_1^2)
  
  for(i in 2:time){
    beta<-rbind(beta,beta[i-1,]+rmvnorm(1,c(0,0),sigma_3^2))
    mu<-rbind(mu,mu[i-1,]+beta[i-1,]+rmvnorm(1,c(0,0),sigma_2^2))
    x<-rbind(x,mu[i,]+rmvnorm(1,c(0,0),sigma_1^2))
  }
  
  
  x<-sig(x)
  x<-data.frame(x)
  colnames(x)<-c("rho","PD")
  print(head(x))
  DR<-c(PD)
  
  plot(x$rho,type="l",ylab=expression(rho))
  plot(x$PD,type="l",ylab=expression(PD))
  
  DR <- {}
  for(i in 1:time){
    #SIR_DRはDR_density.Rの関数
    DR<-c(DR,SIR_DR(L=1,rho=x[i,1],pd=x[i,2])$q)
  }
  
  plot(DR,type="l",ylab=expression(DR))
  out<-data.frame(x,DR)
  out
}

local_trend_sim2<-function(time=100,rho=0.08,PD=0.035,sigma_1=0.001,sigma_2=0.001,sigma_3=0.0005){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  
  beta_0<-0
  mu_0<-c(sig_env(rho))
  beta<-beta_0+rnorm(1,0,sigma_3)
  mu<-beta_0+mu_0+rnorm(1,0,sigma_2)
  x<-mu+rnorm(1,0,sigma_1)
  
  for(i in 2:time){
    beta<-rbind(beta,beta[i-1]+rnorm(1,0,sigma_3))
    mu<-rbind(mu,mu[i-1]+beta[i-1,]+rnorm(1,0,sigma_2))
    x<-rbind(x,mu[i]+rnorm(1,0,sigma_1))
  }
  
  
  x<-sig(x)
  rownames(x)<-NULL
  x<-data.frame(x)
  colnames(x)<-c("rho")
  print(head(x))
  DR<-c(PD)
  
  plot(x$rho,type="l")
  
  out<-data.frame(x,DR,beta,mu)
}

local_trend_sim3<-function(time=100,rho=0.08,PD=0.035,sigma_1=0.001,sigma_2=0.001,sigma_3=0.0005){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  rho_0<-rho
  beta_0<-0
  mu_0<-c(sig_env(PD))
  beta<-beta_0+rnorm(1,0,sigma_3)
  mu<-beta_0+mu_0+rnorm(1,0,sigma_2)
  x<-mu+rnorm(1,0,sigma_1)
  
  for(i in 2:time){
    beta<-rbind(beta,beta[i-1]+rnorm(1,0,sigma_3))
    mu<-rbind(mu,mu[i-1]+beta[i-1,]+rnorm(1,0,sigma_2))
    x<-rbind(x,mu[i]+rnorm(1,0,sigma_1))
  }
  
  
  x<-sig(x)
  rownames(x)<-NULL
  x<-data.frame(x)
  colnames(x)<-c("PD")
  print(head(x))
  DR<-c(PD)
  
  plot(x$PD,type="l",ylab=expression(PD))
  
  for(i in 2:time){
    density<-sapply(1:9999/10000,function(y) g_DR.fn(rho=rho_0,PD=x$PD[i],DR=y))
    density_range<-which(density>0)
    max_denstiy<-max(density)
    
    check<-TRUE
    while(check){
      y<-runif(1,density_range[1]/10000,density_range[length(density_range)]/10000)
      if(g_DR.fn(rho=rho_0,PD=x$PD[i],DR=y) >max_denstiy*runif(1,0,1)){
        DR<-c(DR,y)
        check<-FALSE
      }
    }
  }
  plot(DR,type="l")
  
  
  rownames(beta)<-NULL
  rownames(mu)<-NULL
  colnames(mu)<-NULL
  out<-data.frame(x,DR,beta,mu)
}

copula_approach<-function(time=100,rho=0.08,PD=0.035,sigma_1=0.001,sigma_2=0.001,sigma_3=0.0005){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  
  sigma_1<-diag(c(sigma_1,sigma_1))
  sigma_2<-diag(c(sigma_2,sigma_2))
  sigma_3<-diag(c(sigma_3,sigma_3))
  
  beta_0<-c(0,0)
  mu_0<-c(sig_env(rho),sig_env(PD))
  beta<-beta_0+rmvnorm(1,c(0,0),sigma_3^2)
  mu<-beta_0+mu_0+rmvnorm(1,c(0,0),sigma_2^2)
  x<-mu+rmvnorm(1,c(0,0),sigma_1^2)
  
  for(i in 2:time){
    beta<-rbind(beta,beta[i-1,]+rmvnorm(1,c(0,0),sigma_3^2))
    mu<-rbind(mu,mu[i-1,]+beta[i-1,]+rmvnorm(1,c(0,0),sigma_2^2))
    x<-rbind(x,mu[i,]+rmvnorm(1,c(0,0),sigma_1^2))
  }
  
  
  x<-sig(x)
  x<-data.frame(x)
  colnames(x)<-c("rho","PD")
  print(head(x))
  DR<-c(PD)
  
  plot(x$rho,type="l")
  plot(x$PD,type="l")
  
  
  for(i in 2:time){
    
    n_1_G<-pnorm((sqrt(1-x$rho[i-1])*qnorm(DR[i-1])-qnorm(x$PD[i-1]))/sqrt(x$rho[i-1]))
    
    G_DR_n<-function(DR){
      pnorm((sqrt(1-x$rho[i])*qnorm(DR)-qnorm(x$PD[i]))/sqrt(x$rho[i]))
    }
    
    g_DR_g_n_1<-function(DR){
      dmvnorm(c(qnorm(G_DR_n(DR)),qnorm(n_1_G)),mean=c(0,0),sigma=matrix(c(0.5,0,0,0.5),ncol=2))/n_1_G
    }
    
    
    density<-sapply(1:9999/10000,g_DR_g_n_1)
    density_range<-which(density>0)
    max_density<-max(density[density_range])
    
    check<-TRUE
    while(check){
      
      y<-runif(1,density_range[1]/10000,density_range[length(density_range)]/10000)
      
      
      if(g_DR_g_n_1(y)/
         (max_density)>runif(1,0,1)){
        DR<-c(DR,y)
        print(DR[i])
        check=FALSE
      }
    }
  }
  plot(DR,type="l")
  out=data.frame(x,DR,beta_rho=beta[,1],beta_pd=beta[,2],mu_rho=mu[,1],mu_pd=mu[,2])
}

copula_approach2<-function(time=100,rho=0.08,PD=0.035,sigma_1=0.001,sigma_2=0.001,sigma_3=0.0005){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  
  beta_0<-0
  mu_0<-sig_env(rho)
  beta<-beta_0+rnorm(1,0,sigma_3)
  mu<-beta_0+mu_0+rnorm(1,0,sigma_2)
  x<-mu+rnorm(1,0,sigma_1)
  
  for(i in 2:time){
    beta<-c(beta,beta[i-1]+rnorm(1,0,sigma_3))
    mu<-c(mu,mu[i-1]+beta[i-1,]+rnorm(1,0,sigma_2))
    x<-c(x,mu[i]+rnorm(1,0,sigma_1))
  }
  
  
  x<-sig(x)
  x<-data.frame(x)
  rownames(x)<-NULL
  colnames(x)<-c("rho")
  print(head(x))
  DR<-c(PD)
  
  plot(x$rho,type="l")
  
  
  for(i in 2:time){
    
    n_1_G<-pnorm((sqrt(1-x$rho[i-1])*qnorm(DR[i-1])-qnorm(PD))/sqrt(x$rho[i-1]))
    G_DR_n<-function(DR){
      pnorm((sqrt(1-x$rho[i])*qnorm(DR)-qnorm(PD))/sqrt(x$rho[i]))
    }
    
    g_DR_g_n_1<-function(DR){
      dmvnorm(c(qnorm(G_DR_n(DR)),qnorm(n_1_G)),mean=c(0,0),sigma=matrix(c(0.5,0,0,0.5),ncol=2))/n_1_G
    }
    
    max_density<-max(sapply(1:9999/10000,g_DR_g_n_1))
    check<-TRUE
    while(check){
      y<-runif(1,0,1)
      if(g_DR_g_n_1(y)/
         (max_density)>runif(1,0,1)){
        DR<-c(DR,y)
        print(DR[i])
        check=FALSE
      }
    }
  }
  plot(DR,type="l")
  out=data.frame(x,DR)
  
}

copula_approach3<-function(time=100,rho=0.08,PD=0.035,sigma_1=0.001,sigma_2=0.001,sigma_3=0.0005){
  sig<-function(x){(tanh(x)+1)/2}
  sig_env<-function(y){(1/2)*log(y/(1-y))}
  
  beta_0<-0
  mu_0<-sig_env(PD)
  beta<-beta_0+rnorm(1,0,sigma_3)
  mu<-beta_0+mu_0+rnorm(1,0,sigma_2)
  x<-mu+rnorm(1,0,sigma_1)
  
  for(i in 2:time){
    beta<-c(beta,beta[i-1]+rnorm(1,0,sigma_3))
    mu<-c(mu,mu[i-1]+beta[i-1]+rnorm(1,0,sigma_2))
    x<-c(x,mu[i]+rnorm(1,0,sigma_1))
  }
  
  
  x<-sig(x)
  x<-data.frame(x)
  rownames(x)<-NULL
  colnames(x)<-c("PD")
  print(head(x))
  DR<-c(PD)
  
  plot(x$PD,type="l")
  
  
  for(i in 2:time){
    
    n_1_G<-pnorm((sqrt(1-rho)*qnorm(DR[i-1])-qnorm(x$PD[i-1]))/sqrt(rho))
    G_DR_n<-function(DR){
      pnorm((sqrt(1-rho)*qnorm(DR)-qnorm(x$PD[i]))/sqrt(rho))
    }
    
    g_DR_g_n_1<-function(DR){
      dmvnorm(c(qnorm(G_DR_n(DR)),qnorm(n_1_G)),mean=c(0,0),sigma=matrix(c(0.5,0,0,0.5),ncol=2))/n_1_G
    }
    
    
    
    denstiy<-sapply(1:9999/10000,g_DR_g_n_1)
    density_range<-which(0<denstiy)
    max_density<-max(denstiy[density_range])
    
    
    check<-TRUE
    while(check){
      y<-runif(1,density_range[1]/10000,density_range[length(density_range)]/10000)
      if(g_DR_g_n_1(y)/
         (max_density)>runif(1,0,1)){
        DR<-c(DR,y)
        print(DR[i])
        check=FALSE
      }
    }
  }
  plot(DR,type="l")
  out=data.frame(x,DR,beta,mu)
}
setwd(old_d)