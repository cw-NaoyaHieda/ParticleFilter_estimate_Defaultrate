H1=function(x) x
H2=function(x) (x^2-1)
H3=function(x) (x^3-3*x)
H4=function(x) (x^4-6*x^2+3)
H5=function(x) (x^5-10*x^3+15*x)
H6=function(x) (x^6-15*x^4+45*x^2-15)
g_DR.fn=function(rho,PD,DR){
  out=sqrt((1-rho)/rho)*
    exp(0.5*(qnorm(DR)^2-((sqrt(1-rho)*qnorm(DR)-qnorm(PD))/sqrt(rho))^2))
  return(out)
}

GG_DR.fn<-function(rho,PD,DR){
  pnorm((sqrt(1-rho)*qnorm(DR+0.00001)-qnorm(PD))/sqrt(rho)) -
    pnorm((sqrt(1-rho)*qnorm(DR-0.00001)-qnorm(PD))/sqrt(rho))
}

g_DR.fn_grid<-function(data){
  out=sqrt((1-data[2])/data[2])*
    exp(0.5*(qnorm(data_for_Kalman[1,2])^2-((sqrt(1-data[2])*qnorm(data_for_Kalman[1,2])-qnorm(data[1]))/sqrt(data[2]))^2))
  return(out)
}

g_DR.LF_tmp<-function(data){
  out=sqrt((1-data[2])/data[2])*
    exp(0.5*(qnorm(data[3])^2-((sqrt(1-data[2])*qnorm(data[3])-qnorm(data[1]))/sqrt(data[2]))^2))
  return(out)
}

g_DR.LF<-function(parameters,datas){
  parameters<-matrix(parameters,ncol=2,byrow = T)
  data<-data.frame(parameters,datas)
  out=prod(apply(data,1,g_DR.LF_tmp))
  return(out)
}


Kt.fn=function(PD,rho,C3,C4,Ci3,Ci4,n){
  out1 <- qnorm(PD)
  out2 <- ((rho^(3/2)*C3+(1-rho)^(3/2)*Ci3)/(6*sqrt(n))) * ((qnorm(PD))^2-1)
  out3 <- ((rho^2*C4+(1-rho)^2*Ci4)/(24*n)) * ((qnorm(PD))^3-3*qnorm(PD))
  out4 <- ((rho^3*C3^2+(1-rho)^3*(Ci3^2)+2*(rho*(1-rho))^(3/2)*C3*Ci3)/(36*n)) * ((2*qnorm(PD))^3-5*qnorm(PD))
  return(out1+out2+out3-out4)
}
Gi.fn=function(y,Ci3,Ci4,n){
  out=pnorm(y)-dnorm(y)*( (Ci3/(6*sqrt(n)))*H2(y) +(Ci4/(24*n))*H3(y) +(Ci3^2/(72*n))*H5(y) )
  out
}
G0.invfn=function(DR,rho, C3,C4,Ci3, Ci4,n){
  zp=qnorm(DR)
  out=zp+(zp^2-1)*(rho^(3/2)*C3+(1-rho)^(3/2)*Ci3)/(6*sqrt(n))+(zp^3-3*zp)*(rho^2*C4+(1-rho)^2*Ci4)/(24*n)-(2*(zp^3)-5*zp)*(rho^3*C3^2+(1-rho)^3*Ci3^3+2*(rho*(1-rho))^(3/2)*C3*Ci3)/(36*n)
  out
}
DR_density.fn=function(DR,PD,rho,C3,C4,Ci3,Ci4,n){
  #Kt　現在の累積確率から、正規分布ならどこか
  #Goinv 正規分布で表す　現在の累積確率
  Kt = Kt.fn(PD,rho,C3,C4,Ci3,Ci4,n)
  G0inv = G0.invfn(DR,rho,C3,C4,Ci3, Ci4,n)
  cit = (Kt-sqrt(1-rho)*G0inv)/sqrt(rho)
  out = (sqrt(1-rho)/sqrt(rho))*dnorm(cit)*exp(C3/(6*sqrt(n))*H3(cit)+C4/(24*n)*H4(cit))/(dnorm(G0inv)*exp(Ci3/(6*sqrt(n))*H3(G0inv)+Ci4/(24*n)*H4(G0inv)))
  #out=(sqrt(1-rho)/sqrt(rho))*exp((Kt^2 -2*sqrt(1-rho)*G0inv*Kt +(1-2*rho)*G0inv^2)/(-2*rho))*(1+(C3/(6*sqrt(n)))*H3(cit) +(C4/(24*n))*H4(cit) +(C3^2/(72*n))*H6(cit))*exp(-Ci3/(6*sqrt(n))*H3(G0inv)-Ci4/(24*n)*H4(G0inv))
  out
}

Vasicek_MLE.fn=function(data,ini){
  
  #データのデフォルトレートを代入
  DR=data$Default_Rate
  
  f=function(x){#対数尤度関数
    fn2=function(DR){g_DR.fn(rho=x[1],PD=x[2],DR)}
    LL=-sum(log(sapply(DR,fn2)))
    LL
  }
  
  #この関数いらない気がする
  cf=function(x){
    x
  }
  
  #estimates <- solnp(pars=ini, fun=f,ineqfun=cf,ineqLB=c(0,0),ineqUB=c(1,1))$pars
  estimates <- optim(par=ini, fn=f, method="BFGS",
                     control=list(maxit=10000), hessian=FALSE)
  #反復回数maxit
  names(estimates[[1]])=c("rho","PD")
  
  print(estimates[[2]])
  out=estimates[[1]]
}



Expected_loss=function(PD,rho,C3,C4,Ci3,Ci4,n, lower, upper){
  f_m=function(DR,PD,rho,C3,C4,Ci3,Ci4,n) DR*DR_density.fn(DR,PD,rho,C3,C4,Ci3,Ci4,n)
  integrate(f_m, PD=PD,rho=rho,C3=C3,C4=C4,Ci3=Ci3,Ci4=Ci4,n=n,lower=lower, upper=upper )$value
}

Expected_loss2=function(PD,rho,C3,C4,Ci3,Ci4,n, lower, upper){
  f_m=function(DR,PD,rho,C3,C4,Ci3,Ci4,n) DR^2*DR_density.fn(DR,PD,rho,C3,C4,Ci3,Ci4,n)
  integrate(f_m, PD=PD,rho=rho,C3=C3,C4=C4,Ci3=Ci3,Ci4=Ci4,n=n,lower=lower, upper=upper )$value
}

VaR_loss=function(level,PD,rho,C3,C4,Ci3,Ci4,n,lower){
  quant_L <- function(DR){
    f_m=function(DR,PD,rho,C3,C4,Ci3,Ci4,n){ DR_density.fn(DR,PD,rho,C3,C4,Ci3,Ci4,n)}
    out<-integrate(f_m, PD=PD,rho=rho,C3=C3,C4=C4,Ci3=Ci3,Ci4=Ci4,n=n,lower=lower, upper=DR )$value
    out}
  out <-optimize(function(x) (quant_L(x)-level)^2 , interval=c(1e-16,.9))$minimum
  out}

VaR_loss_2=function(level,PD,rho,lower){
  #quant_L <- function(DR){
  #  f_m=function(DR,PD,rho,C3,C4,Ci3,Ci4,n){ g_DR.fn(rho,PD,DR)}
  #  out<-integrate(f_m, PD=PD,rho=rho,lower=lower, upper=DR )$value
  #  out}
  #out <-optimize(function(x) (quant_L(x)-level)^2 , interval=c(1e-16,.9))$minimum
  j<-0
  for(DR in 1:10000/10000){
    if(level<pnorm((sqrt(1-rho)*qnorm(DR)-qnorm(PD))/sqrt(rho)) & j==0)
    {out=DR
    break}
  }  
  out}


Expected_loss_particle=function(parameter, lower, upper){
  f_m=function(DR,PD,rho) DR*g_DR.fn(rho,PD,DR)
  integrate(f_m,PD=parameter$PD,rho=parameter$rho,lower=lower,upper=upper,subdivisions=1000000,stop.on.error = FALSE)$value
}

Expected_loss2_particle=function(parameter, lower, upper){
  f_m=function(DR,PD,rho) DR^2*g_DR.fn(rho,PD,DR)
  integrate(f_m, PD=parameter$PD,rho=parameter$rho,lower=lower, upper=upper, stop.on.error=FALSE)$value
}

Expected_loss2kai=function(PD,rho,C3,C4,Ci3,Ci4,n,lower,upper){
  #rhoとPDは既知。キュミュラントの計算はどうなる？2変量の分散共分散行列どうする？
  #エルミート掛ける標準正規分布の密度関数->標準正規分布のk回微分
  #これを積分すると->k-1回微分の密度関数×(-1)=k-1次エルミート×(-1)×密度関数になる。
  mu1<-Expected_loss(PD,rho,C3,C4,Ci3,Ci4,n,lower,upper)
  mu<-c(mu1,mu1)#平均
  Kt <- Kt.fn(PD,rho,C3,C4,Ci3,Ci4,n)
  Kt<-c(Kt,Kt)
  Ct<-(Kt-sqrt(rho))/sqrt(1-rho)#式中のyにCt代入。
  Ct<-c(Ct,Ct)
  x<-matrix(c(1,0,0,1),ncol=2)#これの値どうする？
  soukan<-matrix(c(1,rho,rho,1),ncol=2)
  ue<-pmvnorm(lower=-Inf,upper=Kt,mean=mu)#,corr=soukan,sigma=xがあると実行できない。
  #インテグラルの部分
  shita<-C3*H2(Ct)*(-1)*dmvnorm(c(Ct,Ct),mean=rep(mu1,2),sigma=x)/6*sqrt(n)+C4*H3(Ct)*(-1)*dmvnorm(c(Ct,Ct),mean=rep(mu1,2),sigma=x)/24*n+C3*Ci3*H5(Ct)*(-1)*dmvnorm(c(Ct,Ct),mean=rep(mu1,2),sigma=x)/72*n
  out<-ue-shita
  out
}

VaR_loss_particle=function(level,parameter,lower){
  quant_L <- function(DR){
    f_m=function(DR,PD,rho){g_DR.fn(rho,PD,DR)}
    out<-integrate(f_m, PD=parameter$PD,rho=parameter$rho,lower=lower, upper=DR , stop.on.error=FALSE)$value
    out}
  out <-optimize(function(x) (quant_L(x)-level)^2 , interval=c(1e-16,.9))$minimum
  out}

#rho,pd,weight
representative_value.fn<-function(data){
  data_rho<-data[order(data[,1]),]
  data_rho_weight_sum<-cumsum(data_rho[,3])
  tf_label<-0
  t_label<-0
  f_label<-0
  s_label<-0
  nf_label<-0
  seventy_fif<-fifty<-twenty_fif<-two_fif<-dim(data_rho)[1]
  nine_fif<-length(data_rho_weight_sum)
  for(j in 1:length(data_rho_weight_sum)){
    if(data_rho_weight_sum[j]>0.05&tf_label==0){two_fif<-j;tf_label<-1}
    if(data_rho_weight_sum[j]>0.25&t_label==0&data_rho[two_fif,1]<data_rho[j,1]){twenty_fif<-j;t_label<-1}
    if(data_rho_weight_sum[j]>0.5&f_label==0&data_rho[twenty_fif,1]<data_rho[j,1]){fifty<-j;f_label<-1}
    if(data_rho_weight_sum[j]>0.75&s_label==0&data_rho[fifty,1]<data_rho[j,1]){seventy_fif<-j;s_label<-1}
    if(data_rho_weight_sum[j]>0.95&nf_label==0&data_rho[seventy_fif,1]<data_rho[j,1]){nine_fif<-j;nf_label<-1}
  }
  min_rho<-min(data_rho[,1])
  two_fif_rho<-data_rho[two_fif,1]
  twenty_fif_rho<-data_rho[twenty_fif,1]
  fifty_rho<-data_rho[fifty,1]
  seventy_fif_rho<-data_rho[seventy_fif,1]
  nine_fif_rho<-data_rho[nine_fif,1]
  max_rho<-max(data_rho[,1])
  
  data_pd<-data[order(data[,2]),]
  data_pd_weight_sum<-cumsum(data_pd[,3])
  tf_label<-0
  t_label<-0
  f_label<-0
  s_label<-0
  nf_label<-0
  seventy_fif<-fifty<-twenty_fif<-two_fif<-dim(data_pd)[1]
  nine_fif<-length(data_rho_weight_sum)
  for(j in 1:length(data_rho_weight_sum)){
    if(data_pd_weight_sum[j]>0.05&tf_label==0){two_fif<-j;tf_label<-1}
    if(data_pd_weight_sum[j]>0.25&t_label==0&data_pd[two_fif,2]<data_pd[j,2]){twenty_fif<-j;t_label<-1}
    if(data_pd_weight_sum[j]>0.5&f_label==0&data_pd[twenty_fif,2]<data_pd[j,2]){fifty<-j;f_label<-1}
    if(data_pd_weight_sum[j]>0.75&s_label==0&data_pd[fifty,2]<data_pd[j,2]){seventy_fif<-j;s_label<-1}
    if(data_pd_weight_sum[j]>0.95&nf_label==0&data_pd[seventy_fif,2]<data_pd[j,2]){nine_fif<-j;nf_label<-1}
  }
  min_pd<-min(data_pd[,2])
  two_fif_pd<-data_pd[two_fif,2]
  twenty_fif_pd<-data_pd[twenty_fif,2]
  fifty_pd<-data_pd[fifty,2]
  seventy_fif_pd<-data_pd[seventy_fif,2]
  nine_fif_pd<-data_pd[nine_fif,2]
  max_pd<-max(data_pd[,2])
  
  out<-c(min_rho=min_rho,two_fif_rho=two_fif_rho,twenty_fif_rho=twenty_fif_rho,fifty_rho=fifty_rho,seventy_fif_rho=seventy_fif_rho,nine_fif_rho=nine_fif_rho,max_rho=max_rho,
         min_pd=min_pd,two_fif_pd=two_fif_pd,twenty_fif_pd=twenty_fif_pd,fifty_pd=fifty_pd,seventy_fif_pd=seventy_fif_pd,nine_fif_pd=nine_fif_pd,max_pd=max_pd)
  out
}
#x,weight
representative_value.fn2<-function(data){
  data_rho<-data[order(data[,1]),]
  data_rho_weight_sum<-cumsum(data_rho[,2])
  tf_label<-0
  t_label<-0
  f_label<-0
  s_label<-0
  nf_label<-0
  seventy_fif<-twenty_fif<-two_fif<-dim(data_rho)[1]
  nine_fif<-length(data_rho_weight_sum)
  for(j in 1:length(data_rho_weight_sum)){
    if(data_rho_weight_sum[j]>0.05&tf_label==0){two_fif<-j;tf_label<-1}
    if(data_rho_weight_sum[j]>0.25&t_label==0&data_rho[two_fif,1]<data_rho[j,1]){twenty_fif<-j;t_label<-1}
    if(data_rho_weight_sum[j]>0.75&s_label==0&data_rho[twenty_fif,1]<data_rho[j,1]){seventy_fif<-j;s_label<-1}
    if(data_rho_weight_sum[j]>0.95&nf_label==0&data_rho[seventy_fif,1]<data_rho[j,1]){nine_fif<-j;nf_label<-1}
  }
  min_rho<-min(data_rho[,1])
  two_fif_rho<-data_rho[two_fif,1]
  twenty_fif_rho<-data_rho[twenty_fif,1]
  seventy_fif_rho<-data_rho[seventy_fif,1]
  nine_fif_rho<-data_rho[nine_fif,1]
  max_rho<-max(data_rho[,1])
  
  
  out<-c(min=min_rho,two_fif=two_fif_rho,twenty_fif=twenty_fif_rho,seventy_fif=seventy_fif_rho,nine_fif=nine_fif_rho,max=max_rho)
  out
}

representative_value_DR.fn<-function(data){
  density<-sapply(1:9999/10000,function(x) g_DR.fn(DR=x,PD=data[2],rho=data[1]))
  
  probability<-cumsum(density)/sum(density)
  
  two_fif_DR<-c(1:9999/10000)[probability>0.05][1]
  twenty_fif_DR<-c(1:9999/10000)[probability>0.25][1]
  seventy_fif_DR<-c(1:9999/10000)[probability>0.75][1]
  nine_fif_DR<-c(1:9999/10000)[probability>0.95][1]
  ninetynine_DR<-c(1:9999/10000)[probability>0.99][1]
  
  
  out<-c(two_fif_DR,twenty_fif_DR,seventy_fif_DR,nine_fif_DR,ninetynine_DR)
  out
}

auxiliary_particle_filtering<-function(sample_N,data_for_Kalman,proposal_theta_m,sigma_1,sigma_2,sigma_3){
  N<- sample_N
  #シグノイド関数にかけてパラメータの初期値をN個取得
  beta_0<-rnorm(N , mean=proposal_theta_m[2] , sd = 0.01)
  mu_0<-rnorm(N , mean=proposal_theta_m[1] , sd = 0.01)
  theta_PD_0<-mu_0+rnorm(N,mean = 0 ,sd=0.01)
  theta_0<-cbind(theta_PD_0,mu_0,beta_0)
  
  weight_0<-rep(1/N,N)
  colnames(theta_0)<-c("PD","mu","beta")
  #t=1
  #シグノイド関数にかけてパラメータをN個取得
  
  I_data<-data.frame(weight=weight_0,state=sig(theta_0[,2]+theta_0[,3]))
  
  cl <- makeCluster(rep('localhost', 4))
  clusterExport(cl, c("GG_DR.fn","rmvnorm","data_for_Kalman","A_r","rho","N"))
  I_probablity<-parApply(cl,I_data,1,function(x){
    x[1]*GG_DR.fn(rho=rho,PD=x[2],DR=data_for_Kalman[1,2])
  })
  
  probablity<<-cumsum(I_probablity)/sum(I_probablity)
  
  clusterExport(cl, c("probablity"))
  I_k<-parSapply(cl,runif(N,0,1),A_r)
  
  
  theta_k_I_k<<-data.frame(theta_0[I_k,2]+theta_0[I_k,3],theta_0[I_k,2]+theta_0[I_k,3],theta_0[I_k,3])
  clusterExport(cl, c("theta_k_I_k","sigma_1","sigma_2","sigma_3"))
  beta_k<-parApply(cl,theta_k_I_k,1,function(x) rnorm(n=1,mean=c(x[3]),sd=sigma_3))
  mu_k<-parApply(cl,theta_k_I_k,1,function(x) rnorm(n=1,mean=c(x[2]),sd=sigma_2))
  theta_PD_k<-parApply(cl,data.frame(mu_k),1,function(x) rnorm(n=1,mean=c(x),sd=sigma_1))
  theta_k<-cbind(theta_PD_k,mu_k,beta_k)
  
  
  state_data<-data.frame(sig(theta_k[,1]),sig(theta_k_I_k[,1]))
  weight<-parApply(cl,state_data,1,function(x) GG_DR.fn(rho=rho,PD=x[1],DR=data_for_Kalman[1,2])/
                     GG_DR.fn(rho=rho,PD=x[2],DR=data_for_Kalman[1,2]))
  weight<-weight/sum(weight)
  stopCluster(cl)
  
  state<-data.frame(sig(theta_k[,1]),weight,ruiseki=cumsum(weight))
  
  head(state,n=10)
  tail(state,n=10)
  
  #リサンプルでもってくる番号の取得
  re_state_num<-c()
  probablity<-state$ruiseki
  if(sum(weight^2)^-1 < N/10){
    for(i in  1:N){
      re_state_num<-c(re_state_num,A_r((i-1+runif(1))/N))
    }
    re_state<-state[re_state_num,c(1,2)]
    re_theta<-theta_k[re_state_num,]
    re_state[,2]<-rep(1/N,N)
  } else{
    re_state_num<-seq(1,N)
    re_state<-state[,c(1,2)]
    re_theta<-theta_k
  }
  
  
  #一応リサンプルする前もとっておく
  state_100<-list(data.frame(PD=state[,1],weight=weight,re_PD=re_state[,1],re_weight=re_state[,2]))
  theta_100<-list(data.frame(re_theta))
  
  cl <- makeCluster(rep('localhost', 4))
  clusterExport(cl, c("GG_DR.fn","data_for_Kalman","N","rho","rmvnorm"))
  
  ###2~
  for(j in 2:dim(data_for_Kalman)[1]){
    post_theta<-re_theta
    theta<-c()
    state<-c()
    
    I_data<-data.frame(weight=weight,state=sig(post_theta[,2]+post_theta[,3]),data_for_Kalman[j,2])
    
    I_probablity<-parApply(cl,I_data,1,function(x){
      x[1]*GG_DR.fn(rho=rho,PD=x[2],DR=x[3])
    })
    
    probablity<-cumsum(I_probablity)/sum(I_probablity)
    
    tmp<-runif(N,0,1)
    I_k<-pforeach(i = 1:N)({
      A_r(tmp[i])
    })
    
    theta_k_I_k<-data.frame(post_theta[I_k,2]+post_theta[I_k,3],post_theta[I_k,2]+post_theta[I_k,3],post_theta[I_k,3])
    clusterExport(cl, c("theta_k_I_k","sigma_1","sigma_2","sigma_3"))
    beta_k<-parApply(cl,theta_k_I_k,1,function(x) rnorm(n=1,mean=c(x[3]),sd=sigma_3))
    mu_k<-parApply(cl,theta_k_I_k,1,function(x) rnorm(n=1,mean=c(x[2]),sd=sigma_2))
    theta_PD_k<-parApply(cl,data.frame(mu_k),1,function(x) rnorm(n=1,mean=x,sd=sigma_1))
    theta_k<-cbind(theta_PD_k,mu_k,beta_k)
    
    
    
    state_data<-data.frame(sig(theta_k[,1]),sig(theta_k_I_k[,1]),data_for_Kalman[j,2])
    weight<-parApply(cl,state_data,1,function(x) GG_DR.fn(rho=rho,PD=x[1],DR=x[3])/
                       GG_DR.fn(rho=rho,PD=x[2],DR=x[3]))
    weight<-weight/sum(weight)
    
    
    
    state<-data.frame(sig(theta_k[,1]),weight=weight,ruiseki=cumsum(weight))
    re_state_num<-c()
    probablity<-state$ruiseki
    if(sum(weight^2)^-1 < N/10){
      re_state_num<-pforeach(i = 1:N)({
        A_r((i-1+runif(1))/N)
      })
      print("resumple")
      re_state<-state[re_state_num,c(1,2)]
      re_theta<-theta_k[re_state_num,]
      re_state[,2]<-rep(1/N,N)
    } else{
      re_state<-state[,c(1,2)]
      re_theta<-theta_k
    }
    
    state_100<- c(state_100,list(data.frame(PD=state[,1],weight=weight,re_PD=re_state[,1],re_weight=re_state[,2])))
    theta_100<-c(theta_100,list(data.frame(re_theta)))
    
  }
  stopCluster(cl)
  out<-list(state_100,theta_100)
}

particle_smoothing<-function(sample_N,data_for_Kalman,state_100,theta_100,sigma_1,sigma_2,sigma_3){
  N<-sample_N
  weight_T<-list(state_100[[length(data_for_Kalman[,1])]][,4])
  sm_state<-list(data.frame(state_100[[length(data_for_Kalman[,1])]][,c(3,4)]))
  #A関数の設定
  A_r<<-function(r){
    for(i in N:1){
      if(r>=weight_n_cumsum[i]){
        return(i+1)
      }
    }
    return(1)
  }
  
  cl <- makeCluster(rep('localhost', 4))
  for(t in c(length(data_for_Kalman[,1])-1):1){
    weight_n<-c()
    X_n_1<<-data.frame(theta_100[[t+1]],weight=state_100[[t+1]][,4])
    X_n<<-data.frame(theta_100[[t]],weight=state_100[[t]][,4])
    X<-data.frame(X_n_1,X_n,weight_T[[length(data_for_Kalman[,1])-t]])
    colnames(X)<-c("n1_theta_PD","n1_mu","n1_beta","n1_weight","n_theta_PD","n_mu","n_beta","n_weight","weight_T")
    bunbo<-c()
    for(j in 1:N){
      j<<-j
      clusterExport(cl, c("dmvnorm","X_n","X_n_1","sigma_1","sigma_2","sigma_3","j"))
      bunbo<-c(bunbo,sum(parApply(cl,X,1,
                                  function(x) x[4]*dmvnorm(as.numeric(X_n_1[j,c(2,3)]),
                                                           mean=c(x[6]+x[7],x[7]),
                                                           sigma=diag(c(sigma_2^2,sigma_3^2)))*
                                    dnorm(X_n_1[j,1],
                                          mean=X_n_1[j,2],
                                          sd=sigma_1))))
    }
    for(i in 1:N){
      i<<-i
      clusterExport(cl,"i")
      bunsi<-parApply(cl,X,1,function(x) 
        x[9]*dmvnorm(as.numeric(x[c(2,3)]),
                     mean=c(as.numeric(X_n[i,2])+as.numeric(X_n[i,3]),as.numeric(X_n[i,3])),
                     sigma=diag(c(sigma_2^2,sigma_3^2)))*
          dnorm(as.numeric(x[1]),
                mean=as.numeric(x[2]),
                sd=sigma_1)
      )
      weight_n<-c(weight_n,X_n[i,2]*sum(bunsi/bunbo))
    }
    
    
    weight_n<-weight_n/sum(weight_n)
    weight_n_cumsum<-cumsum(weight_n)
    tmp<<-runif(N,0,1)
    sm_state_num<-pforeach(i = 1:N)({
      A_r(tmp[i])
    })
    sm_state<-c(sm_state,list(data.frame(state_100[[t]][c(sm_state_num),3],weight_n[sm_state_num]/sum(weight_n[sm_state_num]))))
    weight_T<-c(weight_T,list(weight_n/sum(weight_n)))
  }
  out<-list(sm_state,weight_T)
}

qq<-function(y,theta0,state_0,sigma_1=sigma_1,sigma_2=sigma_2,sigma_3=sigma_3){
  out=g_DR.fn(rho=0.08,PD=theta0,DR=y)*
    dnorm(as.numeric(state_0[3]),as.numeric(c(state_1[3])),sd=sigma_3)
  print(as.numeric(state_0[1]),as.numeric(c(state_1[3])),sd=sigma_3)
}

Q<-function(par,sm){
 n<-length(sm[[1]])
 nParticle<-dim(sm[[1]][[1]])[1]
 
 sigma_1<-sig(par[1])
 sigma_2<-sig(par[2])
 sigma_3<-sig(par[3])
 
 lik.out<-0
 for(k in n-1:1){
   fn2.a<-function(i,y,q1,q2,q3,sigma_1,sigma_2,sigma_3){
     qq(y,theta0=q1[i],state_0=q2[i,],state_1=q3[i,],sigma_1=sigma_1,sigma_2=sigma_2,sigma_3=sigma_3)
   }
   lik2.sm<-log(sapply(1:nParticle,fn2.a,
              y=data_for_Kalman[n+2-k,2],
              q1=sm[[1]][[k]][,1],
              q2=sm[[2]][[k]],
              q3=sm[[2]][[k-1]],
              sigma_1=sigma_1,sigma_2=sigma_2,sigma_3=sigma_3))
   lik.out[k]<-sum(lik2.sm*sm[[1]][[k]][,2])
 }
 print(sum(lik.out))
 
 return(-sum(lik.out))
   
}

EMestimate<-function(par,maxit,data_for_Kalman,proposal_theta_m){
  #出力準備
  out<-matirx(NA,nrow=maxit+1,ncol=(length(par)+3))
  i<-1
  continue<-T
  old<-par
  
  sigma_1<-sig(old[1])
  sigma_2<-sig(old[2])
  sigma_3<-sig(old[3])
  
  pf<-auxiliary_particle_filtering(data_for_Kalman,proposal_theta_m,sigma_1,sigma_2,sigma_3)
  sm<-particle_smoothing(data_for_Kalman,pf[[1]],pf[[2]],sigma_1,sigma_2,sigma_3)
  
  out[1,1:length(par)+1]=c(old,Q(par=old,sm=sm))
  while(continue){
  # Mstep  
    res <- optim(par=old,fn=Q,method="BFGS",sm=sm)
    print(res)
    new<-res$par
    sigma_1<-sig(new[1])
    sigma_2<-sig(new[2])
    sigma_3<-sig(new[3])
    
    pf<-auxiliary_particle_filtering(data_for_Kalman,proposal_theta_m,sigma_1,sigma_2,sigma_3)
    sm<-particle_smoothing(data_for_Kalman,pf[[1]],pf[[2]],sigma_1,sigma_2,sigma_3)
    
    i<-i+1
    old<-new
    continue<-(i<=maxit&abs(out[i,(length(old)+3)])>eps)
  }
  out=out[1:i,]
  colnames(out)<-c("sigma_1","sigma_2","sigma_3","ll","|diff|","error")
  return(list(est=new,trace=out,pf=pf,sm=sm))
}

get_rho<-function(pd,sigma){
  
  estim_rho<-function(rho){
    rho<-sig(rho)
    v.dr = integrate(function(x) { x^2 * g_DR.fn(rho=rho,PD=pd,x) },
                     lower=0, upper=1)$value - pd^2
    out=(v.dr-sigma)^2
  }
  
  rho<-optim(0.1,estim_rho,method = c("BFGS"))$par
  out<-sig(rho)
  out
}




