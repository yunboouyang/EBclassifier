library("HiDimDA") #FAIR
library("devtools")
library("roxygen2")
install_github("yunboouyang/VBDP")
library("VBDP")
###############functions##################################
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


##########data preparation #####################
set.seed(100)
data(leukemia)
Train=as.matrix(leukemia[1:38,-1]);
Test=as.matrix(leukemia[-(1:38),-1]);
y=leukemia[1:38,1];
trainlabel=y;
testlabel=leukemia[-(1:38),1];

label=levels(y)


n1=27; n2=11;
data1=Train[which(y == label[1]),];
data2=Train[which(y == label[2]),];
n1 = sum(y == label[1])
n2 = sum(y == label[2])

mean1=apply(data1,2,mean)
mean2=apply(data2,2,mean)
var1=apply(data1,2,var)
var2=apply(data2,2,var)
S=sqrt(var1/n1+var2/n2)

u1=t(t(data1)/S);
u2=t(t(data2)/S);
u=t(t(Train)/S);

t=(mean1-mean2)/S


#####################DP########################################
table(DPclassifier(Train,Test,y,1,4,0.9,nfolds=7,sparse=FALSE),testlabel)
table(DPclassifier(Train,Train,y,1,4,0.9,nfolds=7,sparse=FALSE),trainlabel)
# testlabel
# ALL AML
# ALL  20   2
# AML   0  12


# trainlabel
# ALL AML
# ALL  26   0
# AML   1  11


#####################Sparse DP########################################
table(DPclassifier(Train,Test,y,1,4,0.9,nfolds=7),testlabel)
table(DPclassifier(Train,Train,y,1,4,0.9,nfolds=7),trainlabel)
# testlabel
# ALL AML
# ALL  20   2
# AML   0  12

# trainlabel
# ALL AML
# ALL  26   0
# AML   1  11


##############Hard Thresh DP#############################

results=multisparseVBDP(t,1,4,0.995,nfolds=15)
prior=results$prior;
vsparse=dis.CD.hardthresh(t, prior)
a=vsparse/sqrt(sum(vsparse^2))
a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;

tdata = t(t(Test)/S)

pred=ifelse(tdata %*% a + a0 > 0, label[1], label[2])
table(pred,testlabel)
# testlabel
# pred  ALL AML
# ALL  20   2
# AML   0  12

pred=ifelse(u %*% a + a0 > 0, label[1], label[2])
table(pred,trainlabel)

# trainlabel
# pred  ALL AML
# ALL  26   0
# AML   1  11


############################FAIR Approach############################
p=7129

result=Dlda(Train,y,VSelfunct=Fairselect,ldafun="classification");
coef=result$coef;
vkpt=result$vkpt;
if(length(vkpt)==0) 
{print("No features left");} else{
  v=rep(0,p);   #save the coefficients
  v[vkpt]=coef;
  a=v/sqrt(sum(v^2));
  a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
  pred=ifelse(Test %*% a + a0 < 0, label[1], label[2])
  table(pred,testlabel)  
}

pred=ifelse(Train %*% a + a0 < 0, label[1], label[2])
table(pred,trainlabel)
# testlabel
# pred  ALL AML
# ALL  20   6
# AML   0   8

# trainlabel
# pred  ALL AML
# ALL  27   0
# AML   0  11


#############Independence rule##################################
a=t/((var1*n1+var2*n2)/(n1+n2))
a0=-(mean(as.vector(data1%*%a))+mean(as.vector(data2%*%a)))/2;
tdata = t(t(Test)/S)
pred=ifelse(Test %*% a + a0 < 0, label[1], label[2])
table(pred,testlabel)

pred=ifelse(Train %*% a + a0 > 0, label[1], label[2])
table(pred,trainlabel)

# testlabel
# pred  ALL AML
# ALL  20   6
# AML   0   8

# trainlabel
# pred  ALL AML
# ALL  27   1
# AML   0  10


##############EB#############################
h=0.3 #bandwidth parameter

a=fmodel(t,h);
a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;
a0=-(mean(as.vector(u1%*%a))+mean(as.vector(u2%*%a)))/2;

tdata = t(t(Test)/S)

pred=ifelse(tdata %*% a + a0 > 0, label[1], label[2])
table(pred,testlabel)

pred=ifelse(u %*% a + a0 > 0, label[1], label[2])
table(pred,trainlabel)
# testlabel
# pred  ALL AML
# ALL  20   3
# AML   0  11

# trainlabel
# pred  ALL AML
# ALL  27   0
# AML   0  11







