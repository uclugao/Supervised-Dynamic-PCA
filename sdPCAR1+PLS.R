library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
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
######################
for (kk in 1:rp){
  ft=FT(TT,r)
  XX1=Xt(TT,NN,ft,BB,sigma)
  XX=scale(XX1, center = F, scale = TRUE)
  yy=yt(TT,ft,h,sgm,bbeta)
  PLSM=plsr(yy~XX[2:TT,], scale=TRUE, ncomp=2,validation="CV")
  plsft=PLSM$fitted.values[,1,1]
  PLSV[kk]=sqrt(sum((yy-plsft)^2)/TT)
  
}
##########
mean(PLSV)
###
median(PLSV)

####

###out-of sample
####################
TT=200#sample size
NN=500#predictors
r=2#factors
n=470#zero rows
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
######################
for (kk in 1:rp){
  ft=FT(TT,r)
  XX1=Xt(TT,NN,ft,BB,sigma)
  XX=scale(XX1, center = F, scale = TRUE)
  yy=yt(TT,ft,h,sgm,bbeta)
  X1=XX[2:T1,]
  y1=yy[1:(T1-1)]
  
  X2=XX[(T1+1):TT,]
  y2=yy[T1:(TT-1)]
  
  mod1 <- plsr(y1~X1, scale=TRUE, ncomp=2,validation="CV")
  pcr_pred <- predict(mod1, X2, ncomp=1)
  PLSO[kk,]=y2-pcr_pred
  
}
##########
mean(sqrt(rowMeans(PLSO^2)))
############
median(sqrt(rowMeans(PLSO^2)))


