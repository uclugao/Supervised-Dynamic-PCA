library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
library(lars)
library(pls)
##################
#####factors
FT=function(TT,r){
  ft=mvrnorm(TT,rep(0,r),diag(r))
  return(ft)
}
#######X-predictor
Xt=function(TT,NN,ft,B,sigma){
  U=mvrnorm(TT,rep(0,NN),sigma)
  xt=B%*%t(ft)+t(U)
  return(t(xt))
}
######
yt=function(TT,ft,h,sgm,bbeta){
  yt=numeric(TT+h)
  eps=rnorm(TT+h,0,sgm)
  for (i in 1:(TT-1)){
    yt[i+h+1]=t(bbeta)%*%c(ft[i+1,],ft[i,])+eps[i]
  }
  return(yt[(h+2):(TT+h)])
}
####################
TT=200#sample size
NN=500#predictors
r=2#factors
n=490#zero rows
set.seed(1234)
#a1=runif(NN*r/2,0.5,1.5)
#b1=runif(NN*r/2,-1.5,-0.5)
#c1=c(a1,b1)
#BB=matrix(sample(c1,NN*r),nrow=NN,ncol=r)
BB=matrix(runif(NN*r,-2,2),nrow=NN,ncol=r)
#BB=matrix(numeric(NN*r),NN,r)
#BB[,1]=runif(NN,0.5,2)
#BB[,2]=runif(NN,-2,2)
#BB=svd(BB)$u%*%diag(c(2*sqrt(NN),sqrt(NN)))
aa=sample(1:NN,n)
BB[aa,]=0
bbeta=c(1,-0.8,-1,2)
sigma=diag(NN)
#diag(runif(NN,1,3))#predictor
sgm=1#yt
h=1#step
###########test
#ft=FT(TT,r)
#XX=Xt(TT,NN,ft,BB,sigma)
#yy=yt(TT,ft,h,sgm,bbeta)
######errors
rp=100
SDPCA=numeric(rp)
SPCA=numeric(rp)
PCA=numeric(rp)
SW=numeric(rp)
PLSV=numeric(rp)#newly added
ft=matrix(numeric(TT*2),TT,2)
SDPCA9=numeric(rp)
SDPCA10=numeric(rp)
######################
for (kk in 1:rp){
  ft=FT(TT,r)
  XX1=Xt(TT,NN,ft,BB,sigma)
  XX=scale(XX1, center = F, scale = TRUE)
  yy=yt(TT,ft,h,sgm,bbeta)
 
  XNN=matrix(numeric((TT-1)*NN),TT-1,NN)
  for (ii in 1:NN){
    xds=cbind(XX[2:TT,ii],XX[1:(TT-1),ii])
    xds=scale(xds,center=F,scale = TRUE)
    sdmd=lm(yy~xds-1)
    XNN[,ii]=xds%*%sdmd$coefficients
  }
  
  md.lasso=lars(XNN,yy,type="lasso",trace=FALSE,normalize=TRUE,intercept=FALSE)
  bet49=md.lasso$beta[5,]#number of nonzero
  SDPCA9[kk]=sqrt(sum((yy-XNN%*%bet49)^2)/TT)
  indx=which(md.lasso$beta[5,]!=0)
  XNN1=XNN[,indx]
  beta50=lm(yy~XNN1-1)$coefficients
  SDPCA10[kk]=sqrt(sum((yy-XNN1%*%beta50)^2)/TT)
  
}
##########calculate mean median
mean(SDPCA9)
###
median(SDPCA9)
##########
mean(SDPCA10)
###
median(SDPCA10)
####



##############################################################
###out-of sample
####################
TT=200#sample size
NN=300#predictors
r=2#factors
n=260#zero rows
set.seed(1234)
#a1=runif(NN*r/2,0.5,1.5)
#b1=runif(NN*r/2,-1.5,-0.5)
#c1=c(a1,b1)
#BB=matrix(sample(c1,NN*r),nrow=NN,ncol=r)
BB=matrix(runif(NN*r,-2,2),nrow=NN,ncol=r)
#BB=svd(BB)$u%*%diag(c(2*sqrt(NN),sqrt(NN)))
aa=sample(1:NN,n)
BB[aa,]=0
bbeta=c(1,-0.8,-1,2)
sigma=diag(NN)
#diag(runif(NN,1,3))#predictor
sgm=1#yt
h=1#step
###########test
#ft=FT(TT,r)
#XX=Xt(TT,NN,ft,BB,sigma)
#yy=yt(TT,ft,h,sgm,bbeta)
######errors
######errors
rp=100
T1=round(3*TT/5)
T2=TT-T1
PLSO=matrix(numeric(rp*T2),rp,T2)#newly added
SDPCA9=matrix(numeric(rp*T2),rp,T2)#newly added
SDPCA10=matrix(numeric(rp*T2),rp,T2)#newly added
######################
for (kk in 1:rp){
  ft=FT(TT,r)
  XX1=Xt(TT,NN,ft,BB,sigma)
  XX=scale(XX1, center = F, scale = TRUE)
  yy=yt(TT,ft,h,sgm,bbeta)

  for (si in 1:(T2)){
    syy=yy[1:(T1-1+si-1)]
    sXX=XX[1:(T1+1+si-1),]
    XNN=matrix(numeric((T1-1+si-1)*NN),T1+si-1-1,NN)
    sXNN=matrix(numeric((T1-1+si-1+1)*NN),T1+si-1-1+1,NN)
    for (ii in 1:NN){
      xds=cbind(XX[2:(T1+si-1+1),ii],XX[1:(T1+si-1-1+1),ii])
      xds=scale(xds,center=F,scale = TRUE)
      sdmd=lm(syy~xds[1:(T1-1+si-1),]-1)
      #XNN[,ii]=xds%*%sdmd$coefficients
      # sxds=cbind(XX[2:(T1+si-1+1),ii],XX[1:(T1+si-1-1+1),ii])
      #sxds=scale(sxds,center=TRUE,scale = TRUE)
      sXNN[,ii]=xds%*%sdmd$coefficients
    }
    
    md.lasso=lars(sXNN[1:(T1-1+si-1),],syy,type="lasso",trace=FALSE,normalize=TRUE,intercept=FALSE)
    bet49=md.lasso$beta[5,]#number of nonzero
    SDPCA9[kk,si]=yy[T1-1+si-1+1]-t(sXNN[T1+si-1,])%*%bet49
    indx=which(md.lasso$beta[5,]!=0)
    sXNN1=sXNN[,indx]
    beta50=lm(syy~sXNN1[1:(T1-1+si-1),]-1)$coefficients
    SDPCA10[kk,si]=yy[T1-1+si-1+1]-t(sXNN1[T1+si-1,])%*%beta50
  }
}
##########
mean(sqrt(rowMeans(SDPCA9^2)))
############
median(sqrt(rowMeans(SDPCA9^2)))
##########
mean(sqrt(rowMeans(SDPCA10^2)))
############
median(sqrt(rowMeans(SDPCA10^2)))

