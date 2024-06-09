library("foreach")
library("doParallel")
require("ggplot2")
library("expm")
library("mvtnorm")
library("mvnfast")

#kernel
Ker=function(s,h){
  temp=ifelse((s-h)>0,0,ifelse((s+h)<0,0,((1-(abs(s/h))**3)**3)*70/(h*81)))
  return(temp);
}

#kernel density for X
obtain_fX=function(n,X,Z,Y,h,x0){
  temp=sum(apply(Ker(X-x0,h),2,prod), na.rm=TRUE)/n;
  return(temp);
}

#local mean for ghat
obtain_ghat=function(n,X,Z,Y,h,x0,z0){
  numerator=sum(as.numeric(apply(Ker(X-x0,h),2,prod))*as.numeric(apply(Ker(Z-z0,h),2,prod))*Y, na.rm=TRUE);
  denominator=sum(as.numeric(apply(Ker(X-x0,h),2,prod))*as.numeric(apply(Ker(Z-z0,h),2,prod)), na.rm=TRUE);
  if(denominator==0){return(0);}
  return(numerator/denominator);
}
obtain_ghat_cverror=function(n,X,Z,Y,h){
  n=dim(X)[2];
  cverror=0;
  for(j in ((floor(n/2)+1):n)){
    if((j%%10)==0){cverror=cverror+(obtain_ghat(n/2,X[,1:floor(n/2),drop=FALSE],Z[,1:floor(n/2),drop=FALSE],Y[1:floor(n/2)],h,X[,j],Z[,j])-Y[j])**2;
    }}
  return(cverror);
}

#local mean for alphahat
obtain_alphahat=function(n,X,Z,Y,h,x0,ghat){
  numerator=sum(as.numeric(apply(Ker(X-x0,h),2,prod))*(Y-ghat)*(Y-ghat), na.rm=TRUE);
  denominator=sum(apply(Ker(X-x0,h),2,prod), na.rm=TRUE);
  if(denominator==0){return(0);}
  return(numerator/denominator);
}
obtain_alphahat_cverror=function(n,X,Z,Y,h,ghat){
  n=dim(X)[2];
  cverror=0;
  for(j in ((floor(n/2)+1):n)){
    if((j%%10)==0){cverror=cverror+( obtain_alphahat(n/2,X[,1:floor(n/2),drop=FALSE],Z[,1:floor(n/2),drop=FALSE],Y[1:floor(n/2)],h,X[,j],ghat[1:floor(n/2)]) - (Y[j]-ghat[j])*(Y[j]-ghat[j]) )**2;
    }}
  return(cverror);
}

#local mean for dhat
obtain_dhat=function(n,X,Z,Y,h,x0){
  numerator=sum(as.numeric(apply(Ker(X-x0,h),2,prod))*Y, na.rm=TRUE);
  denominator=sum(as.numeric(apply(Ker(X-x0,h),2,prod)), na.rm=TRUE);
  if(denominator==0){return(0);}
  return(numerator/denominator);
}
obtain_dhat_cverror=function(n,X,Z,Y,h){
  n=dim(X)[2];
  cverror=0;
  for(j in ((floor(n/2)+1):n)){
    if((j%%10)==0){cverror=cverror+(obtain_dhat(n/2,X[,1:floor(n/2),drop=FALSE],Z[,1:floor(n/2),drop=FALSE],Y[1:floor(n/2)],h,X[,j])-Y[j])**2
    }}
  return(cverror);
}

#local mean for gammahat
obtain_gammahat=function(n,X,Z,Y,h,x0,dhat){
  numerator=sum(as.numeric(apply(Ker(X-x0,h),2,prod))*(Y-dhat)*(Y-dhat), na.rm=TRUE);
  denominator=sum(apply(Ker(X-x0,h),2,prod), na.rm=TRUE);
  if(denominator==0){return(0);}
  return(numerator/denominator);
}
obtain_gammahat_cverror=function(n,X,Z,Y,h,dhat){
  n=dim(X)[2];
  cverror=0;
  for(j in ((floor(n/2)+1):n)){
    if((j%%10)==0){cverror=cverror+( obtain_gammahat(n/2,X[,1:floor(n/2),drop=FALSE],Z[,1:floor(n/2),drop=FALSE],Y[1:floor(n/2)],h,X[,j],dhat[1:floor(n/2)]) - (Y[j]-dhat[j])*(Y[j]-dhat[j]) )**2;
    }}
  return(cverror);
}


### import data
dat <- readRDS("~/Documents/Project 6 IVI with DaiGuorong/code/dat.RDS"); # this data is normalized
N=dim(dat)[1];



### thetahat
cat("thetahat:","\n")
# different xtarget
real_xtarget=c(33,41,49,57,65);
Age_mu=48.94902; Age_sd=15.77189;
xtarget_list=(real_xtarget-Age_mu)/Age_sd;

for(k in c(1:length(xtarget_list))){
  xtarget=xtarget_list[k];
  
  mento=function(it){
    # wash data
    randat=dat[sample(1:N,N,replace=FALSE),]
    Y=randat[,1]; # Y is SBP
    X=t(as.matrix(randat[,3])); # X is AGE
    Z=t(as.matrix(randat[,2])); # Z is BMI
    p1=dim(X)[1];
    p2=dim(Z)[1];

    ### estimation
    M=2;
    n=N/M; # sample size for ghat, dhat, alphahat and gammahat   
    hmin=0.5;hmax=2; # candidates of the bandwidth
    candilen=5; # number of candidates
    
    X1=X[,1:n,drop=FALSE];     Z1=Z[,1:n,drop=FALSE];     Y1=Y[1:n];
    X2=X[,(n+1):N,drop=FALSE]; Z2=Z[,(n+1):N,drop=FALSE]; Y2=Y[(n+1):N];
    X3=X1; Z3=Z1; Y3=Y1; X4=X2; Z4=Z2; Y4=Y2;
    X5=X; Z5=Z; Y5=Y;
    X1ex=X1; Z1ex=Z1; Y1ex=Y1;
    X2ex=X2; Z2ex=Z2; Y2ex=Y2;
    X3ex=X3; Z3ex=Z3; Y3ex=Y3;
    X4ex=X4; Z4ex=Z4; Y4ex=Y4;
    
    ### estimation for theta
    # data1ex for ghat
    candidate_h_ghat=exp(seq(log(hmin),log(hmax),length.out=candilen));
    h_ghat_cverror=array(0,candilen);
    for(hit in 1:candilen){
      h_ghat_cverror[hit]=obtain_ghat_cverror(dim(X1ex)[2],X1ex,Z1ex,Y1ex,candidate_h_ghat[hit]);
      if(is.na(h_ghat_cverror[hit])){h_ghat_cverror[hit]=10**8;}
    }
    h_ghat1=candidate_h_ghat[which.min(h_ghat_cverror)];
    ghat1=array(0,dim(X2)[2]); # on (X2,Z2)
    for(j in 1:dim(X2)[2]){
      ghat1[j]=obtain_ghat(dim(X1ex)[2],X1ex,Z1ex,Y1ex,h_ghat1,X2[,j],Z2[,j]);
    }
    # data2ex for ghat
    candidate_h_ghat=exp(seq(log(hmin),log(hmax),length.out=candilen));
    h_ghat_cverror=array(0,candilen);
    for(hit in 1:candilen){
      h_ghat_cverror[hit]=obtain_ghat_cverror(dim(X2ex)[2],X2ex,Z2ex,Y2ex,candidate_h_ghat[hit]);
      if(is.na(h_ghat_cverror[hit])){h_ghat_cverror[hit]=10**8;}
    }
    h_ghat2=candidate_h_ghat[which.min(h_ghat_cverror)];
    ghat2=array(0,dim(X1)[2]); # on (X1,Z1)
    for(j in 1:dim(X1)[2]){
      ghat2[j]=obtain_ghat(dim(X2ex)[2],X2ex,Z2ex,Y2ex,h_ghat2,X1[,j],Z1[,j]);
    }
    # data3ex for dhat
    candidate_h_dhat=exp(seq(log(hmin),log(hmax),length.out=candilen));
    h_dhat_cverror=array(0,candilen);
    for(hit in 1:candilen){
      h_dhat_cverror[hit]=obtain_dhat_cverror(dim(X3ex)[2],X3ex,Z3ex,Y3ex,candidate_h_dhat[hit]);
      if(is.na(h_dhat_cverror[hit])){h_dhat_cverror[hit]=10**8;}
    }
    h_dhat3=candidate_h_dhat[which.min(h_dhat_cverror)];
    dhat3=array(0,dim(X4)[2]); # on (X4,Z4)
    for(j in 1:dim(X4)[2]){
      dhat3[j]=obtain_dhat(dim(X3ex)[2],X3ex,Z3ex,Y3ex,h_dhat3,X4[,j]);
    }
    # data4ex for dhat
    candidate_h_dhat=exp(seq(log(hmin),log(hmax),length.out=candilen));
    h_dhat_cverror=array(0,candilen);
    for(hit in 1:candilen){
      h_dhat_cverror[hit]=obtain_dhat_cverror(dim(X4ex)[2],X4ex,Z4ex,Y4ex,candidate_h_dhat[hit]);
      if(is.na(h_dhat_cverror[hit])){h_dhat_cverror[hit]=10**8;}
    }
    h_dhat4=candidate_h_dhat[which.min(h_dhat_cverror)];
    dhat4=array(0,dim(X3)[2]); # on (X3,Z3)
    for(j in 1:dim(X3)[2]){
      dhat4[j]=obtain_dhat(dim(X4ex)[2],X4ex,Z4ex,Y4ex,h_dhat4,X3[,j]);
    }
    
    ###
    h2=max(min(h_ghat1,h_ghat2,h_dhat3,h_dhat4)*0.5,hmin);
    
    # ghat1 and data2 for alphahat3
    h_alphahat3=h2;
    alphahat3_on_xtarget=obtain_alphahat(dim(X2)[2],X2,Z2,Y2,h_alphahat3,xtarget,ghat1);
    # ghat2 and data1 for alphahat4
    h_alphahat4=h2;
    alphahat4_on_xtarget=obtain_alphahat(dim(X1)[2],X1,Z1,Y1,h_alphahat4,xtarget,ghat2);
    # dhat3 and data4 for gammahat3
    h_gammahat3=h2;
    gammahat3_on_xtarget=obtain_gammahat(dim(X4)[2],X4,Z4,Y4,h_gammahat3,xtarget,dhat3);
    # dhat4 and data3 for gammahat4
    h_gammahat4=h2;
    gammahat4_on_xtarget=obtain_gammahat(dim(X3)[2],X3,Z3,Y3,h_gammahat4,xtarget,dhat4);
    
    # thetahat
    thetahat_on_xtarget=(alphahat3_on_xtarget+alphahat4_on_xtarget)/(gammahat3_on_xtarget+gammahat4_on_xtarget);
    
    return(as.list(c(thetahat_on_xtarget,thetahat_on_xtarget)));
  }# end of mento
  
  #main function
  closeAllConnections();
  closeAllConnections();
  cl <- makeCluster(100);
  registerDoParallel(cl);
  itermax=100;
  result <- foreach(it=c(1:itermax), .combine='c') %dopar% mento(it);
  stopCluster(cl);
  thetahat_mento=array(0,itermax);
  for(it in 1:itermax){
    thetahat_mento[it]=result[[2*it-1]];
  }
  output=c(real_xtarget[k], xtarget,
           round(mean(thetahat_mento,na.rm=TRUE),4),
           round(sd(thetahat_mento,na.rm=TRUE),4)  );
  cat(output,"\n");
}# end of xtarget



### thetahatstar
cat("thetahatstar:","\n")
# different xtarget
real_xtarget=c(33,41,49,57,65);
Age_mu=48.94902; Age_sd=15.77189;
xtarget_list=(real_xtarget-Age_mu)/Age_sd;

for(k in c(1:length(xtarget_list))){
  xtarget=xtarget_list[k];
  # wash data
  set.seed(1093);
  randat=dat[sample(1:N,N,replace=FALSE),]
  Y=randat[,1]; # Y is SBP
  X=t(as.matrix(randat[,3])); # X is AGE
  Z=t(as.matrix(randat[,2])); # Z is BMI
  p1=dim(X)[1];
  p2=dim(Z)[1];
  A=runif(10**7,-1,1);
  Knorm=(mean(Ker(A,1)*Ker(A,1))*2)**p1; #square of L2 norm
  rm(A);
  tau=0.025; # confidence level 1-2tau
  
  M=2;
  n=N/(2*M); # sample size for ghat, dhat, alphahat and gammahat 
  hmin=0.5;hmax=2; # candidates of the bandwidth
  candilen=5; # number of candidates
  X1=X[,(0*n+1):(1*n),drop=FALSE]; Z1=Z[,(0*n+1):(1*n),drop=FALSE]; Y1=Y[(0*n+1):(1*n)];
  X2=X[,(1*n+1):(2*n),drop=FALSE]; Z2=Z[,(1*n+1):(2*n),drop=FALSE]; Y2=Y[(1*n+1):(2*n)];
  X3=X[,(2*n+1):(3*n),drop=FALSE]; Z3=Z[,(2*n+1):(3*n),drop=FALSE]; Y3=Y[(2*n+1):(3*n)];
  X4=X[,(3*n+1):(4*n),drop=FALSE]; Z4=Z[,(3*n+1):(4*n),drop=FALSE]; Y4=Y[(3*n+1):(4*n)];
  X5=X; Z5=Z; Y5=Y;
  X1ex=cbind(X1,X3,X4);Z1ex=cbind(Z1,Z3,Z4);Y1ex=c(Y1,Y3,Y4);
  X2ex=cbind(X2,X3,X4);Z2ex=cbind(Z2,Z3,Z4);Y2ex=c(Y2,Y3,Y4);
  X3ex=cbind(X1,X2,X3);Z3ex=cbind(Z1,Z2,Z3);Y3ex=c(Y1,Y2,Y3);
  X4ex=cbind(X1,X2,X4);Z4ex=cbind(Z1,Z2,Z4);Y4ex=c(Y1,Y2,Y4);
    
  # data1ex for ghat
  candidate_h_ghat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_ghat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_ghat_cverror[hit]=obtain_ghat_cverror(dim(X1ex)[2],X1ex,Z1ex,Y1ex,candidate_h_ghat[hit]);
    if(is.na(h_ghat_cverror[hit])){h_ghat_cverror[hit]=10**8;}
  }
  h_ghat1=candidate_h_ghat[which.min(h_ghat_cverror)];
  ghat1=array(0,dim(X2)[2]); # on (X2,Z2)
  for(j in 1:dim(X2)[2]){
    ghat1[j]=obtain_ghat(dim(X1ex)[2],X1ex,Z1ex,Y1ex,h_ghat1,X2[,j],Z2[,j]);
  }
  # data2ex for ghat
  candidate_h_ghat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_ghat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_ghat_cverror[hit]=obtain_ghat_cverror(dim(X2ex)[2],X2ex,Z2ex,Y2ex,candidate_h_ghat[hit]);
    if(is.na(h_ghat_cverror[hit])){h_ghat_cverror[hit]=10**8;}
  }
  h_ghat2=candidate_h_ghat[which.min(h_ghat_cverror)];
  ghat2=array(0,dim(X1)[2]); # on (X1,Z1)
  for(j in 1:dim(X1)[2]){
    ghat2[j]=obtain_ghat(dim(X2ex)[2],X2ex,Z2ex,Y2ex,h_ghat2,X1[,j],Z1[,j]);
  }
  # data3ex for dhat
  candidate_h_dhat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_dhat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_dhat_cverror[hit]=obtain_dhat_cverror(dim(X3ex)[2],X3ex,Z3ex,Y3ex,candidate_h_dhat[hit]);
    if(is.na(h_dhat_cverror[hit])){h_dhat_cverror[hit]=10**8;}
  }
  h_dhat3=candidate_h_dhat[which.min(h_dhat_cverror)];
  dhat3=array(0,dim(X4)[2]); # on (X4,Z4)
  for(j in 1:dim(X4)[2]){
    dhat3[j]=obtain_dhat(dim(X3ex)[2],X3ex,Z3ex,Y3ex,h_dhat3,X4[,j]);
  }
  # data4ex for dhat
  candidate_h_dhat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_dhat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_dhat_cverror[hit]=obtain_dhat_cverror(dim(X4ex)[2],X4ex,Z4ex,Y4ex,candidate_h_dhat[hit]);
    if(is.na(h_dhat_cverror[hit])){h_dhat_cverror[hit]=10**8;}
  }
  h_dhat4=candidate_h_dhat[which.min(h_dhat_cverror)];
  dhat4=array(0,dim(X3)[2]); # on (X3,Z3)
  for(j in 1:dim(X3)[2]){
    dhat4[j]=obtain_dhat(dim(X4ex)[2],X4ex,Z4ex,Y4ex,h_dhat4,X3[,j]);
  }
  
  #
  h2=max(min(h_ghat1,h_ghat2,h_dhat3,h_dhat4)*0.5,hmin);
  # ghat1 and data2 for alphahat3
  h_alphahat3=h2;
  alphahat3_on_xtarget=obtain_alphahat(dim(X2)[2],X2,Z2,Y2,h_alphahat3,xtarget,ghat1);
  # ghat2 and data1 for alphahat4
  h_alphahat4=h2;
  alphahat4_on_xtarget=obtain_alphahat(dim(X1)[2],X1,Z1,Y1,h_alphahat4,xtarget,ghat2);
  # dhat3 and data4 for gammahat3
  h_gammahat3=h2;
  gammahat3_on_xtarget=obtain_gammahat(dim(X4)[2],X4,Z4,Y4,h_gammahat3,xtarget,dhat3);
  # dhat4 and data3 for gammahat4
  h_gammahat4=h2;
  gammahat4_on_xtarget=obtain_gammahat(dim(X3)[2],X3,Z3,Y3,h_gammahat4,xtarget,dhat4);
  
  # thetastarhat
  thetastarhat_on_xtarget=(alphahat3_on_xtarget+alphahat4_on_xtarget)/(gammahat3_on_xtarget+gammahat4_on_xtarget);
  
  # data5 for gtilde
  candidate_h_ghat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_ghat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_ghat_cverror[hit]=obtain_ghat_cverror(dim(X5)[2],X5,Z5,Y5,candidate_h_ghat[hit]);
    if(is.na(h_ghat_cverror[hit])){h_ghat_cverror[hit]=10**8;}
  }
  h_ghat=candidate_h_ghat[which.min(h_ghat_cverror)];
  gtilde=array(0,dim(X5)[2]); # on (X5,Z5)
  for(j in 1:dim(X5)[2]){
    gtilde[j]=obtain_ghat(dim(X5)[2],X5,Z5,Y5,h_ghat,X5[,j],Z5[,j]);
  }
  # data5 for alphatilde
  candidate_h_alphahat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_alphahat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_alphahat_cverror[hit]=obtain_alphahat_cverror(dim(X5)[2],X5,Z5,Y5,candidate_h_alphahat[hit],gtilde);
    if(is.na(h_alphahat_cverror[hit])){h_alphahat_cverror[hit]=10**8;}
  }
  h_alphahat=candidate_h_alphahat[which.min(h_alphahat_cverror)];
  alphatilde=array(0,dim(X5)[2]); # on (X5,Z5)
  for(j in 1:dim(X5)[2]){
    alphatilde[j]=obtain_alphahat(dim(X5)[2],X5,Z5,Y5,h_alphahat,X5[,j],gtilde);
  }
  alphatilde_on_xtarget=obtain_alphahat(dim(X5)[2],X5,Z5,Y5,h_alphahat,xtarget,gtilde);
  # data5 for dtilde
  candidate_h_dhat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_dhat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_dhat_cverror[hit]=obtain_dhat_cverror(dim(X5)[2],X5,Z5,Y5,candidate_h_dhat[hit]);
    if(is.na(h_dhat_cverror[hit])){h_dhat_cverror[hit]=10**8;}
  }
  h_dhat=candidate_h_dhat[which.min(h_dhat_cverror)];
  dtilde=array(0,dim(X5)[2]); # on (X5,Z5)
  for(j in 1:dim(X5)[2]){
    dtilde[j]=obtain_dhat(dim(X5)[2],X5,Z5,Y5,h_dhat,X5[,j]);
  }
  # data5 for gammatilde
  candidate_h_gammahat=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_gammahat_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_gammahat_cverror[hit]=obtain_gammahat_cverror(dim(X5)[2],X5,Z5,Y5,candidate_h_gammahat[hit],dtilde);
    if(is.na(h_gammahat_cverror[hit])){h_gammahat_cverror[hit]=10**8;}
  }
  h_gammahat=candidate_h_gammahat[which.min(h_gammahat_cverror)];
  gammatilde=array(0,dim(X5)[2]); # on (X5,Z5)
  for(j in 1:dim(X5)[2]){
    gammatilde[j]=obtain_gammahat(dim(X5)[2],X5,Z5,Y5,h_gammahat,X5[,j],dtilde);
  }
  gammatilde_on_xtarget=obtain_gammahat(dim(X5)[2],X5,Z5,Y5,h_gammahat,xtarget,dtilde);
  
  #kernel density for f_{X} on x_target
  fXtilde_on_xtarget=obtain_fX(dim(X5)[2],X5,Z5,Y5,h_dhat,xtarget);
  
  #V1tilde on x_target, similar formular to dhat
  equivY1=((Y5-gtilde)*(Y5-gtilde)-alphatilde)*((Y5-gtilde)*(Y5-gtilde)-alphatilde);
  candidate_h_V1tilde=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_V1tilde_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_V1tilde_cverror[hit]=obtain_dhat_cverror(dim(X5)[2],X5,Z5,equivY1,candidate_h_V1tilde[hit]);
    if(is.na(h_V1tilde_cverror[hit])){h_V1tilde_cverror[hit]=10**8;}
  }
  h_V1tilde=candidate_h_V1tilde[which.min(h_V1tilde_cverror)];
  V1tildet_on_xtarget=obtain_dhat(dim(X5)[2],X5,Z5,equivY1,h_V1tilde,xtarget);
  #V2tilde on x_target, similar formular to dhat
  equivY2=((Y5-dtilde)*(Y5-dtilde)-gammatilde)*((Y5-dtilde)*(Y5-dtilde)-gammatilde);
  candidate_h_V2tilde=exp(seq(log(hmin),log(hmax),length.out=candilen));
  h_V2tilde_cverror=array(0,candilen);
  for(hit in 1:candilen){
    h_V2tilde_cverror[hit]=obtain_dhat_cverror(dim(X5)[2],X5,Z5,equivY2,candidate_h_V2tilde[hit]);
    if(is.na(h_V2tilde_cverror[hit])){h_V2tilde_cverror[hit]=10**8;}
  }
  h_V2tilde=candidate_h_V2tilde[which.min(h_V2tilde_cverror)];
  V2tildet_on_xtarget=obtain_dhat(dim(X5)[2],X5,Z5,equivY2,h_V2tilde,xtarget);
  #Vtilde on x_target
  Vtildet_on_xtarget=(2*Knorm*V1tildet_on_xtarget)/(fXtilde_on_xtarget*gammatilde_on_xtarget*gammatilde_on_xtarget)+(2*Knorm*V2tildet_on_xtarget*alphatilde_on_xtarget*alphatilde_on_xtarget)/(fXtilde_on_xtarget*gammatilde_on_xtarget*gammatilde_on_xtarget*gammatilde_on_xtarget*gammatilde_on_xtarget)
  #confidence interval
  radius=qnorm(1-tau,0,1)*sqrt(Vtildet_on_xtarget/(N*(h2**p1)));
  #cut off
  ciupper=min(thetastarhat_on_xtarget+radius,1);
  cilower=max(thetastarhat_on_xtarget-radius,0);
   
  output=c(real_xtarget[k], xtarget,
         round(thetastarhat_on_xtarget,4),
         round(cilower,4),
         round(ciupper,4)  );
  cat(output,"\n");
}
# end of xtarget



