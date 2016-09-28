library("HiDimDA") #FAIR
library("glmnet") #use Logistic regression with Lasso
library("foreach")

set.seed(100)
#general parameters
p=10^3*10;
delta=c(1,1,1,1.5,2,2.5,3,3.5,4)
l=c(200,100,50,30,20,10,5,5,4)*10
h=0.3 #bandwidth parameter
sparseerror=matrix(0,nrow=9,ncol=7)
sumerror=matrix(0,nrow=9,ncol=7)
aveerror=matrix(0,nrow=9,ncol=7)





s=5/sqrt(2)
grouping=as.factor(c(rep(1,25),rep(-1,25)))



#f-modeling procedure, where h is bandwidth parameter, which returns a value.
fmodel<-function(z,h){
  denom=colSums(dnorm(outer(z,z,"-")/h));
  numer=colSums(outer(z,z,"-")*dnorm(outer(z,z,"-")/h)/h^2);
  v=z+numer/denom;
  a=v/sqrt(sum(v^2));
  return(a);
}

Fairselect<-function(data,grouping){
  nvkpt=SelectV(data,grouping,Selmethod="Fair")$nvkpt;
  vkptInd=SelectV(data,grouping,Selmethod="Fair")$vkptInd
  return(list(nvkpt=nvkpt,vkptInd=vkptInd))
}


dis.CD = function(x, prior.mass, h=1){
  # prior on mu is a discrete mixture over values provided in prior.mass
  # bandwith = h
  A=as.data.frame(table(prior.mass));
  freq=A$Freq;
  uniq=as.numeric(levels(A$prior.mass))[A$prior.mass]
  n=length(prior.mass)
  tmp = outer(x, uniq, '-'); 
  tmp = exp(-tmp^2/(2*h));
  tmp=t(t(tmp)*freq)
  tmp = tmp/rowSums(tmp); 
  return(tmp %*% uniq)
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
    phi=E/rowSums(E);
  }
  mean=c(0,tau1/tau2);
  #mean[which(abs(mean)<0.25)]=0;
  
  newphi=cbind(phi%*%p,phi%*%diag(1-p));
  #newphi[,1]=rowSums(newphi[,which(abs(mean)<0.25)]);
  zeroprob=newphi[,1];   
  number=max.col(newphi);
  prior=mean[number];
  csize=length(unique(number));
  return(list(prior=prior,csize=csize,prob=zeroprob));
}











for(time in 1:100){
  
  
  for (j in 1:9){ # Different sparcity
    #data generating process
    
    #mu=c(rep(delta[j],l[j]),rep(0,p-l[j])); #sparse model
    mu=c(rep(delta[j],l[j]),rnorm(p-l[j],sd=0.1)); # nonsparse model
    tau=rep(0,p);
    data1=matrix(0,nrow=25,ncol=p);  # data group 1, 25*10000 matrix
    data2=matrix(0,nrow=25,ncol=p);  # data group 2
    for (n1 in 1:25){                  #25 data points for each group
      data1[n1,]=mu+rnorm(p,sd=s)
      data2[n1,]=rnorm(p,sd=s)
    }
    
    
    #Empirical Bayes method (EB)
    meandiff=apply(data1,2,mean)-apply(data2,2,mean);
    s1=apply(data1,2,var);
    s2=apply(data2,2,var);
    S=sqrt((s1+s2)/25)
    u1=data1%*%diag(1/S);
    u2=data2%*%diag(1/S);
    
    z=meandiff/S
    v=fmodel(z,h);
    a=v/sqrt(sum(v^2));
    a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
    sparseerror[j,4]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))
    
    
    
    #DP model
    
    results=multisparseVBDP(z,1,4,0.01,nfolds=10)
    prior=results$prior;
    v=dis.CD(z,prior);
        if(v==0){sparseerror[j,3]=0.5;}
        else { a=v/sqrt(sum(v^2))
    a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
    sparseerror[j,3]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))}
    
    
    #Sparse DP model
    
    results=multisparseVBDP(z,1,4,0.01,nfolds=50)
    prior=results$prior;
    v=dis.CD.sparse(z,prior);
    a=v/sqrt(sum(v^2))
    a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
    sparseerror[j,2]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))
    
    
    #Sparse DP model
    
    results=multisparseVBDP(z,1,4,0.01,nfolds=50)
    prior=results$prior;
    v=dis.CD.hardthresh(z,prior);
    a=v/sqrt(sum(v^2))
    a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
    sparseerror[j,1]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))
    
    
    
    
    
    #Independence rule
    data=rbind(data1,data2);
    coef=Dlda(data, grouping, VSelfunct="none",ldafun="classification")$coef
    a=coef/sqrt(sum(coef^2))
    a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
    
    sparseerror[j,5]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))
    
    
    
    #FAIR Approach
    result=Dlda(data,grouping,VSelfunct=Fairselect,ldafun="classification");
    coef=result$coef;
    vkpt=result$vkpt;
    if(length(vkpt)==0) 
    {sparseerror[j,6]=0.5;} else{
      v=rep(0,p);   #save the coefficients
      v[vkpt]=coef
      a=v/sqrt(sum(v^2))
      
      a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
      sparseerror[j,6]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))
      
    }
    
    #Logistic regression with Lasso
    obj=cv.glmnet(data,grouping,family="binomial")
    lassocoef=coef(obj)[2:(p+1)]
    if (sum(lassocoef^2)==0)
    {sparseerror[j,7]=0.5;} else {
      a=lassocoef/sqrt(sum(lassocoef^2));
      a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
      sparseerror[j,7]=0.5*pnorm((-t(a/S)%*%mu-a0)/s)+0.5*(1-pnorm((-t(a/S)%*%tau-a0)/s))
    }
    
  }
  
  print(time)
  sumerror=sumerror+sparseerror
}



aveerror=sumerror/100
aveerror=round(aveerror,4)

save(aveerror,file="/home/youyang4/fast classification/simu12.RData")



