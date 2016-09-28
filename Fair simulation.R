library("smoothmest")
library("HiDimDA") #FAIR
library("glmnet") #use Logistic regression with Lasso
set.seed(1000)
#Genarate alpha
p=4500;
ntrain=30;
ntest=200;
error=matrix(0,100,5)   #change
grouping=as.factor(c(rep(1,ntrain),rep(-1,ntrain)))
data1=matrix(0,nrow=ntrain,ncol=p);
data2=matrix(0,nrow=ntrain,ncol=p);
tdata1=matrix(0,nrow=ntest,ncol=p);
tdata2=matrix(0,nrow=ntest,ncol=p);





Fairselect<-function(data,grouping){
  nvkpt=SelectV(data,grouping,Selmethod="Fair")$nvkpt;
  vkptInd=SelectV(data,grouping,Selmethod="Fair")$vkptInd
  return(list(nvkpt=nvkpt,vkptInd=vkptInd))
}





dis.CD.hardthresh = function(x, prior.mass, h=1){ #t is the original value
  n=length(prior.mass)
  A=as.data.frame(table(prior.mass));
  freq=A$Freq;
  uniq=as.numeric(levels(A$prior.mass))[A$prior.mass]
  
  
  if(sum(uniq==0)==0){
    return(x)
  }else {
    location=which(uniq==0)
    tmp = outer(x, uniq, '-'); 
    tmp = exp(-tmp^2/(2*h));
    tmp=t(t(tmp)*freq)
    tmp = tmp/rowSums(tmp); 
    post=ifelse(tmp[,location]>0.5,0,x)
    return(post)
  }
}





dis.CD.sparse = function(x, prior.mass, h=1){
  # prior on mu is a discrete mixture over values provided in prior.mass
  # bandwith = h
  #check whether it is 0 or not
  n=length(prior.mass)
  A=as.data.frame(table(prior.mass));
  freq=A$Freq;
  uniq=as.numeric(levels(A$prior.mass))[A$prior.mass]
  
  if(sum(uniq==0)==0){
    tmp = outer(x, uniq, '-'); 
    tmp = exp(-tmp^2/(2*h));
    tmp=t(t(tmp)*freq)
    tmp = tmp/rowSums(tmp); 
    return(tmp %*% uniq)
  }else {
    location=which(uniq==0)
    tmp = outer(x, uniq, '-'); 
    tmp = exp(-tmp^2/(2*h));
    tmp=t(t(tmp)*freq)
    tmp = tmp/rowSums(tmp); 
    post=tmp %*% uniq;
    post[which(tmp[,location]>0.5)]=0;
    return(post)
  }
  
}



multisparseVBDP=function(x,alpha,sigma, w, T0=10,nfolds=10){
  n=length(x);
  wholeprior=rep(0,n) #record whole prior 
  wholeprob=rep(0,n);
  if(n%%nfolds!=0){stop("vector length should be a multiple of nfolds.")}
  foldid=sample(rep(seq(nfolds), times =n%/%nfolds))
  for(i in seq(nfolds)){
    which= foldid==i
    result=sparseVBDP(x[which],alpha,sigma, w, T0=T0)
    prior=result$prior;
    wholeprior[which]=prior
    prob=result$prob;
    wholeprob[which]=prob;
  }
  csize=length(unique(wholeprior));
  return(list(prior=wholeprior,csize=csize,prob=wholeprob));
}



sparseVBDP=function(x,alpha, sigma, w, T0=10){
  n=length(x);
  lambda=1/sigma^2;
  phi0=matrix(rep(0,n*T0),nrow=n,ncol=T0);
  phi=matrix(rep(1/T0,n*T0),nrow=n,ncol=T0);
  while(max(abs(phi-phi0))>10^-3){
    phi0=phi;
    gamma1=1+colSums(phi)[1:(T0-1)];
    gamma2=alpha+rev(cumsum(rev(colSums(phi)))[1:(T0-1)]);
    tau1=(t(phi)%*%x)[1:T0,];
    tau2=lambda+colSums(phi);
    odds=log(w)-log(1-w)+log(sqrt(1/lambda*colSums(phi)+1))-tau1^2/(2*tau2);
    p=exp(odds)/(1+exp(odds))
    d1=c(digamma(gamma1)-digamma(gamma1+gamma2),0);
    d2=digamma(gamma2)-digamma(gamma1+gamma2);
    d3=d1+c(0,cumsum(d2));
    d=d3-0.5*(1-p)*((tau1/tau2)^2+1/tau2);
    S=outer(x,(1-p)*tau1/tau2,'*')+outer(rep(1,n),d,'*');
    E=exp(S);
    phi=diag((1/rowSums(E)))%*%E;
  }
  mean=c(0,tau1/tau2);
  mean[which(abs(mean)<0.6)]=0;
  newphi=cbind(phi%*%p,phi%*%diag(1-p));
  newphi[,1]=rowSums(newphi[,which(abs(mean)<0.6)]);
  zeroprob=newphi[,1];   
  number=max.col(newphi);
  prior=mean[number];
  csize=length(unique(number));
  return(list(prior=prior,csize=csize,prob=zeroprob));
}





for (time in 1:100){ # Different sparcity
  #data generating process
  B=rbinom(n=p,size=1,prob=rep(0.02,p));
  dx=rdoublex(p,lambda=0.5);
  mu=B*dx;  #alpha is fixed
  tau=0;
  
  a=runif(p,max=0.4);
  a1=c(a[1:(p/3)],rep(0,2*p/3));
  a2=c(rep(0,p/3),a[(p/3+1):(2*p/3)],rep(0,p/3));
  a3=c(rep(0,2*p/3),a[(2*p/3+1):p]);
  rm(a)
  b=runif(p,max=0.2);
  chi1=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  chi2=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  chi3=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  chi4=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  for(i in 1:ntrain){
    nz=rnorm(p);
    data1[i,]=mu+(nz+chi1[i]*a1+chi2[i]*a2+chi3[i]*a3+chi4[i]*b)/sqrt(1+a1^2+a2^2+a3^2+b^2)
    
  }
  
  
  
  chi1=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  chi2=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  chi3=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  chi4=(rchisq(ntrain,df=6)-6)/sqrt(2*6);
  for(i in 1:ntrain){
    nz=rnorm(p);
    data2[i,]=(nz+chi1[i]*a1+chi2[i]*a2+chi3[i]*a3+chi4[i]*b)/sqrt(1+a1^2+a2^2+a3^2+b^2)
  }
  
  
  
  chi1=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  chi2=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  chi3=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  chi4=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  for(i in 1:ntest){
    nz=rnorm(p);
    tdata1[i,]=mu+(nz+chi1[i]*a1+chi2[i]*a2+chi3[i]*a3+chi4[i]*b)/sqrt(1+a1^2+a2^2+a3^2+b^2)
  }
  
  
  
  chi1=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  chi2=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  chi3=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  chi4=(rchisq(ntest,df=6)-6)/sqrt(2*6);
  for(i in 1:ntest){
    nz=rnorm(p);
    tdata2[i,]=(nz+chi1[i]*a1+chi2[i]*a2+chi3[i]*a3+chi4[i]*b)/sqrt(1+a1^2+a2^2+a3^2+b^2)
  }
  
  
#Sparse DP model
  zerolocation=which(B==0);
  meandiff=apply(data1,2,mean)-apply(data2,2,mean);
  s1=apply(data1,2,var);
  s2=apply(data2,2,var);
  S=sqrt((s1+s2)/ntrain)
  u1=t(t(data1)/S);
  u2=t(t(data2)/S);
  
  z=meandiff/S
  tu1=t(t(tdata1)/S);
  tu2=t(t(tdata2)/S);
  
  
  ####Oracle procedure
  
  a=z;
  #Oracle Procedure
  a[zerolocation]=0;
  a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
  pre1=tdata1%*%a+a0;
  pre2=tdata2%*%a+a0;
  error[time,1]=(sum(pre1<0)+sum(pre2>0))/(2*ntest);
  
  
  
  ###Hard Thresh DP
  
  result=multisparseVBDP(z,1,6,0.01,nfolds=50)
  prior=result$prior;
  v=dis.CD.hardthresh(z,prior);
  if(sqrt(sum(v^2))==0){
    error[time,2]=0.5;
  } else{
  a=v/sqrt(sum(v^2))
  a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
  
  
  pre1=tu1%*%a+a0;
  pre2=tu2%*%a+a0;
  error[time,2]=(length(which(pre1<0))+length(which(pre2>0)))/(2*ntest);
  }

  
  ###Sparse DP
  result=multisparseVBDP(z,1,6,0.01,nfolds=50)
  prior=result$prior;
  v=dis.CD.sparse(z,prior);
  if(sqrt(sum(v^2))==0){
    error[time,3]=0.5;
  } else{
    a=v/sqrt(sum(v^2))
    a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
    
    
    pre1=tu1%*%a+a0;
    pre2=tu2%*%a+a0;
    error[time,3]=(sum(pre1<0)+sum(pre2>0))/(2*ntest);
    
  }

  
  
#FAIR Approach
  data=rbind(data1,data2);
  result=Dlda(data,grouping,VSelfunct=Fairselect,ldafun="classification");
  coef=result$coef;
  vkpt=result$vkpt;
  if(length(vkpt)==0) 
  {error[time,4]=0.5;} else{
    v=rep(0,p);   #save the coefficients
    v[vkpt]=coef;
    a=v/sqrt(sum(v^2));
    a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
    pre1=tdata1%*%a+a0;
    pre2=tdata2%*%a+a0;
    error[time,4]=(sum(pre1<0)+sum(pre2>0))/(2*ntest);
  }
  
  
  #Logistic regression with Lasso
  obj=cv.glmnet(data,grouping,family="binomial")
  lassocoef=coef(obj)[2:(p+1)]
  if (sum(lassocoef^2)==0)
  {error[time,5]=0.5;} else {
    a=lassocoef/sqrt(sum(lassocoef^2));
    a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
    pre1=tdata1%*%a+a0;
    pre2=tdata2%*%a+a0;
    error[time,5]=(sum(pre1<0)+sum(pre2>0))/(2*ntest);
    
  }

  print(time)
}

save(error,file="/home/youyang4/fast classification/Fair simulation error.RData")




boxplot.matrix(error,names=c("Oracle","Hard Thresh DP","Sparse DP","FAIR","glmnet"),outline=FALSE, range=0)
stripchart(data.frame(error),method="jitter",vertical=TRUE,pch=20,add=T)


