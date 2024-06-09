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
N=6000;
for(M in c(2,5,10)){
for(scorenumber in c(1:2)){ # 1 for uniform and 2 for Normal
for(xtarget_number in c(3:5)){
for(modelnumber in c(1:3)){ # 1 for Y=X, 2 for Y=X+Z, 3 for Y=X+Z+Phi
for(p2 in c(2)){
          
if(xtarget_number==1){p1=1;xtarget=0;} # dimension for X, estimate theta(xtarget)
if(xtarget_number==2){p1=1;xtarget=1;}
if(xtarget_number==3){p1=2;xtarget=c(0,0);}
if(xtarget_number==4){p1=2;xtarget=c(1,0);}
if(xtarget_number==5){p1=2;xtarget=c(0,1);}

A=runif(10**7,-1,1);
Knorm=(mean(Ker(A,1)*Ker(A,1))*2)**p1; #square of L2 norm
rm(A);

tau=0.025; # confidence level 2tau
n=N/(2*M); # sample size for ghat, dhat, alphahat and gammahat   
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
  totalX=array(0,c(2*M,p1,n));
  totalZ=array(0,c(2*M,p2,n));
  totalY=array(0,c(2*M,n));
  allX=array(0,c(p1,N));
  allZ=array(0,c(p2,N));
  allY=array(0,c(N));
  for(m in 1:(2*M)){
    if(scorenumber==1){
      X=matrix(runif(n*p1,-sqrt(3),sqrt(3)), nrow=p1, ncol=n); 
      Z=matrix(runif(n*p2,-sqrt(3),sqrt(3)), nrow=p2, ncol=n);
    }
    if(scorenumber==2){
      XZ=t(rmvnorm(n, as.numeric(array(0,p1+p2)), Sigma))
      X=XZ[1:p1,,drop=FALSE];
      Z=XZ[(p1+1):(p1+p2),,drop=FALSE];
    }
    e=rnorm(n,0,sigma);
    Y=array(0,n);
    for(j in 1:n){
      if(modelnumber==1){Y[j]=sum(X[,j])+e[j];}
      if(modelnumber==2){Y[j]=sum(X[,j])+sum(Z[,j])+e[j];}
      if(modelnumber==3){Y[j]=sum(X[,j])+sum(Z[,j])+e[j]*2*pnorm(xtarget[1]);}
    }
    for(j in 1:n){
      totalX[m,,j]=X[,j];
      totalZ[m,,j]=Z[,j];
      totalY[m,j]=Y[j];
      allX[,(m-1)*n+j]=X[,j];
      allZ[,(m-1)*n+j]=Z[,j];
      allY[(m-1)*n+j]=Y[j];
    }
    rm(X);rm(Z);rm(Y);
  }
  
  Xex=array(0,c(2*M,p1,N-n));
  Zex=array(0,c(2*M,p2,N-n));
  Yex=array(0,c(2*M,N-n));
  for(m in 1:M){
    temp=c(c(1:((m-1)*n)),c((m*n+1):(M*n)));
    if(m==1){temp=c((m*n+1):(M*n));}
    if(m==M){temp=c(1:((m-1)*n));}
    temp=c(temp,c((M*n+1):N));
    for(k in 1:p1){
      Xex[m,k,]=allX[k,temp];    
    }
    for(k in 1:p2){
      Zex[m,k,]=allZ[k,temp];    
    }
    Yex[m,]=allY[temp];
  }
  for(m in (M+1):(2*M)){
    temp=c(c((M*n+1):((m-1)*n)),c((m*n+1):N));
    if(m==(M+1)){temp=c((m*n+1):N);}
    if(m==(2*M)){temp=c((M*n+1):((m-1)*n));}
    temp=c(c(1:(M*n)),temp);
    for(k in 1:p1){
      Xex[m,k,]=allX[k,temp];    
    }
    for(k in 1:p2){
      Zex[m,k,]=allZ[k,temp];    
    }
    Yex[m,]=allY[temp];
  }

  ### estimation for thetastar
  alphahat3_on_xtarget=array(0,M);
  gammahat3_on_xtarget=array(0,M);
  for(m in 1:M){
    X1ex=array(0,c(p1,N-n));
    Z1ex=array(0,c(p2,N-n));
    Y1ex=array(0,c(N-n));
    X2=array(0,c(p1,n));
    Z2=array(0,c(p2,n));
    Y2=array(0,c(n));
    X3ex=array(0,c(p1,N-n));
    Z3ex=array(0,c(p2,N-n));
    Y3ex=array(0,c(N-n));
    X4=array(0,c(p1,n));
    Z4=array(0,c(p2,n));
    Y4=array(0,c(n));
    for(k in 1:p1){
      X1ex[k,]=Xex[m,k,];
      X2[k,]=totalX[m,k,];
      X3ex[k,]=Xex[M+m,k,];
      X4[k,]=totalX[M+m,k,];
    }
    for(k in 1:p2){
      Z1ex[k,]=Zex[m,k,];
      Z2[k,]=totalZ[m,k,];
      Z3ex[k,]=Zex[M+m,k,];
      Z4[k,]=totalZ[M+m,k,];
    }
    Y1ex=Yex[m,];
    Y2=totalY[m,];
    Y3ex=Yex[M+m,];
    Y4=totalY[M+m,];
    #ghat
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
    #dhat
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
    ### h2
    h2=max(min(h_ghat1,h_dhat3)*0.5,hmin);
    # ghat1 and data2 for alphahat3
    h_alphahat3=h2;
    alphahat3_on_xtarget[m]=obtain_alphahat(dim(X2)[2],X2,Z2,Y2,h_alphahat3,xtarget,ghat1);
    # dhat3 and data4 for gammahat3
    h_gammahat3=h2;
    gammahat3_on_xtarget[m]=obtain_gammahat(dim(X4)[2],X4,Z4,Y4,h_gammahat3,xtarget,dhat3);
  }
  
  # thetastarhat
  thetastarhat_on_xtarget=sum(alphahat3_on_xtarget)/sum(gammahat3_on_xtarget);
  
  
  ### estimation for asymptotic variance from data5
  X5=allX;
  Z5=allZ;
  Y5=allY;
  
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
  
  ### output
  return(as.list(c(thetastarhat_on_xtarget,radius)));
}#end of mento
        
#################################### main function
timestart<-Sys.time();
closeAllConnections();
closeAllConnections();
cl <- makeCluster(100);
registerDoParallel(cl);
result <- foreach(it=c(1:itermax), .combine='c') %dopar% mento(it);
stopCluster(cl);
closeAllConnections();
closeAllConnections();
timeend<-Sys.time();
runningtime<-as.numeric(timeend-timestart);

thetahat_mento=array(0,itermax);
radius_mento=array(0,itermax);
for(it in 1:itermax){
  thetahat_mento[it]=result[[2*it-1]];
  radius_mento[it]=result[[2*it]];
}

# ture value of theta
true_theta_on_xtarget=1;
if(modelnumber==2){
  if(scorenumber==1){true_theta_on_xtarget=(sigma*sigma)/(p2+sigma*sigma);}
  if(scorenumber==2){
    Sigma11=Sigma[1:p1,1:p1,drop=FALSE];
    Sigma12=Sigma[1:p1,(p1+1):(p1+p2),drop=FALSE];
    Sigma21=Sigma[(p1+1):(p1+p2),1:p1,drop=FALSE];
    Sigma22=Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2),drop=FALSE];
    SigmaZonX=Sigma22-Sigma21%*%solve(Sigma11)%*%Sigma12;
    true_theta_on_xtarget=(sigma*sigma)/(sum(SigmaZonX)+sigma*sigma);
  }
}
if(modelnumber==3){
  rs=sigma*sigma*2*pnorm(xtarget[1])*2*pnorm(xtarget[1]);
  if(scorenumber==1){true_theta_on_xtarget=(rs)/(p2+rs);}
  if(scorenumber==2){
    Sigma11=Sigma[1:p1,1:p1,drop=FALSE];
    Sigma12=Sigma[1:p1,(p1+1):(p1+p2),drop=FALSE];
    Sigma21=Sigma[(p1+1):(p1+p2),1:p1,drop=FALSE];
    Sigma22=Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2),drop=FALSE];
    SigmaZonX=Sigma22-Sigma21%*%solve(Sigma11)%*%Sigma12;
    true_theta_on_xtarget=(rs)/(sum(SigmaZonX)+rs);
  }
}

cover=as.numeric(abs(thetahat_mento-true_theta_on_xtarget)<radius_mento);
cilength=array(0,itermax);
for(it in 1:itermax){
  upper=min(thetahat_mento[it]+radius_mento[it],1);
  lower=max(thetahat_mento[it]-radius_mento[it],0);
  cilength[it]=upper-lower;
}
output=c(round(mean(cover,na.rm=TRUE),4),
         round(mean(cilength,na.rm=TRUE),4),
         round(runningtime,4));
cat(output,' ');
if((modelnumber==3)&&(p2==2)){cat("\n");}

}}}}} # end of cases

