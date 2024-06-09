library("foreach")
library("doParallel")
require("ggplot2")
library("expm")
library("mvtnorm")

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


itermax=1000;
for(N in c(1000,2000,3000)){
for(xtarget_number in c(1:5)){
for(scorenumber in c(1:2)){ # 1 for uniform and 2 for Normal
for(mit in c(1:2)){ # marginal IVI, only model 2 for Y=X+Z, only p2=2
modelnumber=2;
p2=2;
if(xtarget_number==1){p1=1;xtarget=0;} # dimension for X, estimate theta(xtarget)
if(xtarget_number==2){p1=1;xtarget=1;}
if(xtarget_number==3){p1=2;xtarget=c(0,0);}
if(xtarget_number==4){p1=2;xtarget=c(1,0);}
if(xtarget_number==5){p1=2;xtarget=c(0,1);}

  M=2;
  n=N/M; # sample size for ghat, dhat, alphahat and gammahat   
  hmin=0.5;hmax=2; # candidates of the bandwidth
  candilen=5; # number of candidates
  sigma=1; # sd of the Y-g(U) 
  Sigma=array(0,c(p1+p2,p1+p2)); # variance matrix for normal distribution
  for(i in c(1:(p1+p2))){for(j in c(1:(p1+p2))){Sigma[i,j]=0.2**(abs(i-j));}}
  
  ################################### mento carlo run
  mento=function(it){
    library("foreach")
    library("doParallel")
    require("ggplot2")
    library("expm")
    library("mvtnorm")

    # generate data
    # data1, the first half
    if(scorenumber==1){
      X=matrix(runif(n*p1,-sqrt(3),sqrt(3)), nrow=p1, ncol=n); 
      Z=matrix(runif(n*p2,-sqrt(3),sqrt(3)), nrow=p2, ncol=n);
    }
    if(scorenumber==2){
      XZ=t(rmvnorm(n, as.numeric(array(0,p1+p2)), Sigma))
      X=XZ[1:p1,,drop=FALSE];Z=XZ[(p1+1):(p1+p2),,drop=FALSE];
    }
    e=rnorm(n,0,sigma);
    Y=array(0,n);
    for(j in 1:n){
      if(modelnumber==1){Y[j]=sum(X[,j])+e[j];}
      if(modelnumber==2){Y[j]=sum(X[,j])+sum(Z[,j])+e[j];}
      if(modelnumber==3){Y[j]=sum(X[,j])+sum(Z[,j])+e[j]*2*pnorm(xtarget[1]);}
    }
    X1=X;rm(X);
    Z1=Z[mit,,drop=FALSE];rm(Z);
    Y1=Y;rm(Y);
    # data2, the second half
    if(scorenumber==1){
      X=matrix(runif(n*p1,-sqrt(3),sqrt(3)), nrow=p1, ncol=n); 
      Z=matrix(runif(n*p2,-sqrt(3),sqrt(3)), nrow=p2, ncol=n);
    }
    if(scorenumber==2){
      XZ=t(rmvnorm(n, as.numeric(array(0,p1+p2)), Sigma))
      X=XZ[1:p1,,drop=FALSE];Z=XZ[(p1+1):(p1+p2),,drop=FALSE];
    }
    e=rnorm(n,0,sigma);
    Y=array(0,n);
    for(j in 1:n){
      if(modelnumber==1){Y[j]=sum(X[,j])+e[j];}
      if(modelnumber==2){Y[j]=sum(X[,j])+sum(Z[,j])+e[j];}
      if(modelnumber==3){Y[j]=sum(X[,j])+sum(Z[,j])+e[j]*2*pnorm(xtarget[1]);}
    }
    X2=X;rm(X);
    Z2=Z[mit,,drop=FALSE];rm(Z);
    Y2=Y;rm(Y);
    # data3, the first half
    X3=X1;
    Z3=Z1;
    Y3=Y1;
    # data4, the second half
    X4=X2;
    Z4=Z2;
    Y4=Y2;
    # data5=data1+data2, the whole data
    X5=cbind(X1,X2);
    Z5=cbind(Z1,Z2);
    Y5=c(Y1,Y2);
    #
    X1ex=X1;Z1ex=Z1;Y1ex=Y1;
    X2ex=X2;Z2ex=Z2;Y2ex=Y2;
    X3ex=X3;Z3ex=Z3;Y3ex=Y3;
    X4ex=X4;Z4ex=Z4;Y4ex=Y4;
    
    ### estimation
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
    
    ### output
    return(as.list(c(thetahat_on_xtarget,thetahat_on_xtarget)));
}#end of mento
  
  
  #################################### main function
  closeAllConnections();
  closeAllConnections();
  cl <- makeCluster(100);
  registerDoParallel(cl);
  result <- foreach(it=c(1:itermax), .combine='c') %dopar% mento(it);
  stopCluster(cl);
  closeAllConnections();
  closeAllConnections();
  thetahat_mento=array(0,itermax);
  for(it in 1:itermax){
    thetahat_mento[it]=result[[2*it-1]];
  }
  
  # true theta
  if(scorenumber==1){true_theta=(1+sigma*sigma)/(p2+sigma*sigma);}
  if(scorenumber==2){
    Sigma11=Sigma[1:p1,1:p1,drop=FALSE];
    Sigma12=Sigma[1:p1,(p1+1):(p1+p2),drop=FALSE];
    Sigma21=Sigma[(p1+1):(p1+p2),1:p1,drop=FALSE];
    Sigma22=Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2),drop=FALSE];
    muZomX=Sigma21%*%solve(Sigma11)%*%xtarget;
    SigmaZonX=Sigma22-Sigma21%*%solve(Sigma11)%*%Sigma12;
    denominator=sum(SigmaZonX)+sigma*sigma;
    if(mit==1){numerator=SigmaZonX[2,2]-SigmaZonX[2,1]*SigmaZonX[1,2]/SigmaZonX[1,1]+sigma*sigma;}
    if(mit==2){numerator=SigmaZonX[1,1]-SigmaZonX[1,2]*SigmaZonX[2,1]/SigmaZonX[2,2]+sigma*sigma;}
    true_theta=numerator/denominator;
    
  } 
  as=abs(thetahat_mento-true_theta);
  output=c(round(mean(as,na.rm=TRUE),4),
           round(sd(as,na.rm=TRUE),4)
           );
  cat(output," ");
  if((scorenumber==2)&&(mit==2)){cat("\n");}
  
}}}} # end of cases







