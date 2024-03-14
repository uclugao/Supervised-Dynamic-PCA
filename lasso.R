library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
library(lars)
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
set.seed(1234)
TT=300#sample size
NN=500#predictors
r=2#factors
n=460#zero rows
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
bbeta=c(6,3,-5,0)
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
frq=numeric(rp)
######################
for (kk in 1:rp){
  ft=FT(TT,r)
  XX=Xt(TT,NN,ft,BB,sigma)
  yy=yt(TT,ft,h,sgm,bbeta)
  #XX=scale(XX1, center = F, scale = TRUE)
  XNN=matrix(numeric((TT-1)*NN),TT-1,NN)
  BG=matrix(numeric(4*NN),4,NN)
  for (ii in 1:NN){
    xds=cbind(XX[2:TT,ii],XX[1:(TT-1),ii])
    #xds=scale(xds,center=F,scale = TRUE)
    sdmd=lm(yy~xds-1)
    XNN[,ii]=xds%*%sdmd$coefficients
    BG[,ii]=(sdmd$coefficients)%x%BB[ii,]
  }
  sdPCA=XNN%*%t(XNN)
  ft4=eigen(sdPCA)$vectors[,1:4]
  UP=svd(t(BG))$v
  fts=ft4%*%t(UP)
  sdpd=lm(yy~fts-1)
  md.lasso=lars(fts,yy,type="lasso",trace=FALSE,normalize=TRUE,intercept=FALSE)
  id=which(md.lasso$beta[4,]!=0)
  frq[kk]=I(id[1]==1&id[2]==2&id[3]==3)
}
##########
mean(frq)



