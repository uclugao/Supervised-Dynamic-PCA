library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
##################
r=2
TT=500
ft=mvrnorm(TT,rep(0,r),diag(r))
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
NN=300#predictors
r=2#factors
n=260#zero rows
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
ft=matrix(numeric(TT*2),TT,2)
dd=matrix(numeric(rp*NN),rp,NN)
######################
qmax=5
for (kk in 1:rp){
  ft=FT(TT,r)
  XX1=Xt(TT,NN,ft,BB,sigma)
  XX=scale(XX1, center = F, scale = TRUE)
  yy=yt(TT,ft,h,sgm,bbeta)
  
  for (ii in 1:NN){
    XC=NULL
    for (jj in 1:qmax){
      XC=cbind(XC,XX[(qmax-jj+1):(TT-jj+1),ii])
    }
    
    aicc=numeric(qmax+1)
    aicc[1]=AIC(lm(yy[(qmax-1):(TT-1)]~1))
    for (ll in 1:qmax){
      aicc[ll+1]=AIC(lm(yy[(qmax-1):(TT-1)]~XC[,1:ll]))
    }
    #model=lm(yy[(qmax-1):(TT-1)]~XC)
    
   # AIC(lm(yy[(qmax-1):(TT-1)]~XC[,1:5]))
   # model=lm(yy[4:(TT-1)]~XX[5:TT,ii]+XX[4:(TT-1),ii]+XX[3:(TT-2),ii]+XX[2:(TT-3),ii]+XX[1:(TT-4),ii])
     #tmm=stepAIC(model, direction="both",trace=FALSE)
    #dd[kk,ii]=length(tmm$coefficients)-1
    dd[kk,ii]=which.min(aicc)-1
  }
  
}
##########

tdd=numeric(NN)
sct=setdiff(c(1:NN),aa)
tdd[sct]=2
dv1=matrix(numeric(rp*NN),rp,NN)
dv2=matrix(numeric(rp*NN),rp,NN)
for (rr in 1:rp ){
 dv1[rr,]=I(dd[rr,]==tdd)
 dv2[rr,]=I(dd[rr,]>=tdd)
}
##############
mean(dv1)
mean(dv2)




