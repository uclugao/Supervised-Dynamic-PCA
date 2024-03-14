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
T1=round(3*TT/5)
T2=TT-T1
SDPCA=matrix(numeric(rp*T2),rp,T2)
SPCA=matrix(numeric(rp*T2),rp,T2)
PCA=matrix(numeric(rp*T2),rp,T2)
SW=matrix(numeric(rp*T2),rp,T2)
######################
for (kk in 1:rp){
  ft=FT(TT,r)
  XX1=Xt(TT,NN,ft,BB,sigma)
  XX=scale(XX1, center = F, scale = TRUE)
  yy=yt(TT,ft,h,sgm,bbeta)
  for (si in 1:(T2)){
    syy=yy[1:(T1-1+si-1)]
    sXX=XX[1:(T1+1+si-1),]
    PCV=sXX%*%t(sXX)
    ft1=eigen(PCV)$vectors[,1:2]
    swd=lm(syy~ft1[2:(T1+si-1),]-1)
    bet1=swd$coefficients
    SW[kk,si]=yy[T1-1+si-1+1]-t(bet1)%*%ft1[T1+si-1+1,]
    
    ft2=cbind(ft1[2:(T1+si-1),],ft1[1:(T1-1+si-1),])
    pcad=lm(syy~ft2-1)
    bet2=pcad$coefficients
    #ft2s=cbind(ft1s[2:(T1+si),],ft1s[1:(T1+si-1),])
    PCA[kk,si]=yy[T1-1+si-1+1]-t(bet2)%*%c(ft1[T1+si,],ft1[T1+si-1,])
   
     Hucef=numeric(NN)
    for (jj in 1:NN){
      mmd=lm(syy~XX[2:(T1+si-1),jj]-1)
      Hucef[jj]=mmd$coefficients
    }
    sXN=XX[1:(T1+1+si-1),]%*%diag(Hucef)
    sPCA=sXN%*%t(sXN)
    ft3=eigen(sPCA)$vectors[,1:2]
    spd=lm(syy~ft3[2:(T1+si-1),]-1)
    bet3=spd$coefficients
    SPCA[kk,si]=yy[T1-1+si-1+1]-t(ft3[T1+1+si-1,])%*%bet3
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
    sdPCA=sXNN%*%t(sXNN)
    ft4=eigen(sdPCA)$vectors[,1:4]
    sdpd=lm(syy~ft4[1:(T1-1+si-1),]-1)
    bet4=sdpd$coefficients
    SDPCA[kk,si]=yy[T1-1+si-1+1]-t(ft4[T1+si-1,])%*%bet4
     
    }
}
##########
mean(sqrt(rowMeans(SDPCA^2)))
mean(sqrt(rowMeans(PCA^2)))
mean(sqrt(rowMeans(SPCA^2)))
mean(sqrt(rowMeans(SW^2)))
############
median(sqrt(rowMeans(SDPCA^2)))
median(sqrt(rowMeans(PCA^2)))
median(sqrt(rowMeans(SPCA^2)))
median(sqrt(rowMeans(SW^2)))

