library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
###########
#Mdata=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/macro-data/FREW-MD-md.csv")
data2=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/sp-data/predictor-Monthly-new.csv")
y5=data2[55:744,15]#variance-volatility
y6=data2[55:744,17]#return
#dim(Mdata)
#tdata=data2[55:744,3:15]
fdata=data2[55:744,3:15]-data2[54:743,3:15]
y5=fdata[,13]
#dim(fdata)
#ts.plot(fdata[,123])
y1=fdata[,6]#IP
y2=fdata[,24]#UNRATE
y3=fdata[,103]#CPI
y4=fdata[,123]#VIX
######
TT=690
#NN=123
NN=13
T1=round(4*TT/5)
T2=TT-T1
SDPCA1=numeric(T2)
SDPCA2=numeric(T2)
SDPCA3=numeric(T2)
SDPCA4=numeric(T2)
SDPCA5=numeric(T2)
SDPCA6=numeric(T2)
SDPCA7=numeric(T2)
SPCA1=numeric(T2)
SPCA2=numeric(T2)
SPCA3=numeric(T2)
SPCA4=numeric(T2)
PCA1=numeric(T2)
PCA2=numeric(T2)
PCA3=numeric(T2)
PCA4=numeric(T2)
SW1=numeric(T2)
SW2=numeric(T2)
SW3=numeric(T2)
SW4=numeric(T2)
ARR=numeric(T2)
ARR=numeric(T2)
ARR2=numeric(T2)
######################
h=1
q=2
XX=fdata
XX=scale(XX, center = F, scale = TRUE)
yy=y6
for (si in 1:(T2)){
  syy=yy[(q+h):(T1+si-1)]#regression
  sXX=XX[1:(T1-h+1+si-1),]
  sXX=scale(sXX, center = F, scale = TRUE)
  PCV=sXX%*%t(sXX)
  ftt=eigen(PCV)$vectors
  ft11=ftt[,1]
  swd11=lm(syy~ft11[(q):(T1-h+si-1)]-1)
  bet11=swd11$coefficients
  SW1[si]=yy[T1+si-1+1]-bet11*ft11[T1-h+si-1+1]
  
  ft12=ftt[,1:2]
  swd12=lm(syy~ft12[(q):(T1-h+si-1),]-1)
  bet12=swd12$coefficients
  SW2[si]=yy[T1+si-1+1]-t(bet12)%*%ft12[T1-h+si-1+1,]
  
  ft13=ftt[,1:3]
  swd13=lm(syy~ft13[(q):(T1-h+si-1),]-1)
  bet13=swd13$coefficients
  SW3[si]=yy[T1+si-1+1]-t(bet13)%*%ft13[T1-h+si-1+1,]
  
  
  ft21=cbind(ftt[q:(T1-h+si-1),1],ftt[(q-1):(T1-h-1+si-1),1])
  pcad21=lm(syy~ft21-1)
  bet21=pcad21$coefficients
  #ft2s=cbind(ft1s[2:(T1+si),],ft1s[1:(T1+si-1),])
  PCA1[si]=yy[T1+si-1+1]-t(bet21)%*%c(ftt[T1-h+1+si-1,1],ftt[T1-h+si-1,1])
  
  ft22=cbind(ftt[q:(T1-h+si-1),1:2],ftt[(q-1):(T1-h-1+si-1),1:2])
  pcad22=lm(syy~ft22-1)
  bet22=pcad22$coefficients
  #ft2s=cbind(ft1s[2:(T1+si),],ft1s[1:(T1+si-1),])
  PCA2[si]=yy[T1+si-1+1]-t(bet22)%*%c(ftt[T1-h+1+si-1,1:2],ftt[T1-h+si-1,1:2])
  
  ft23=cbind(ftt[q:(T1-h+si-1),1:3],ftt[(q-1):(T1-h-1+si-1),1:3])
  pcad23=lm(syy~ft23-1)
  bet23=pcad23$coefficients
  #ft2s=cbind(ft1s[2:(T1+si),],ft1s[1:(T1+si-1),])
  PCA3[si]=yy[T1+si-1+1]-t(bet23)%*%c(ftt[T1-h+1+si-1,1:3],ftt[T1-h+si-1,1:3])
  
  
  Hucef=numeric(NN)
  for (jj in 1:NN){
    mmd=lm(syy~XX[(q):(T1-h+si-1),jj]-1)
    Hucef[jj]=mmd$coefficients
  }
  sXN=XX[1:(T1-h+1+si-1),]%*%diag(Hucef)
  sPCA=sXN%*%t(sXN)
  ft3=eigen(sPCA)$vectors
  ft31=ft3[,1]
  spd31=lm(syy~ft31[q:(T1-h+si-1)]-1)
  bet31=spd31$coefficients
  SPCA1[si]=yy[T1+si-1+1]-ft31[T1-h+1+si-1]*bet31
  
  ft32=ft3[,1:2]
  spd32=lm(syy~ft32[q:(T1-h+si-1),]-1)
  bet32=spd32$coefficients
  SPCA2[si]=yy[T1+si-1+1]-t(ft32[T1-h+1+si-1,])%*%bet32
  
  ft33=ft3[,1:3]
  spd33=lm(syy~ft33[q:(T1-h+si-1),]-1)
  bet33=spd33$coefficients
  SPCA3[si]=yy[T1+si-1+1]-t(ft33[T1-h+1+si-1,])%*%bet33
  ######
  
  sXNN=matrix(numeric((T1-h-q+2+si-1)*NN),T1-h-q+2+si-1,NN)
  #sXNN=matrix(numeric((T1-1+si-1+1)*NN),T1+si-1-1+1,NN)
  for (ii in 1:NN){
    xds=NULL
    for (xd in 1:q){
      xds=cbind(xds,XX[(q-xd+1):(T1-h+1-xd+1+si-1),ii])
    }
    xds=scale(xds,center=F,scale = TRUE)
    sdmd=lm(syy~xds[1:(T1-h-q+1+si-1),]-1)
    #XNN[,ii]=xds%*%sdmd$coefficients
    # sxds=cbind(XX[2:(T1+si-1+1),ii],XX[1:(T1+si-1-1+1),ii])
    #sxds=scale(sxds,center=TRUE,scale = TRUE)
    sXNN[,ii]=xds%*%sdmd$coefficients
  }
  sdPCA=sXNN%*%t(sXNN)
  
  ft4=eigen(sdPCA)$vectors
  ft41=ft4[,1]
  sdpd41=lm(syy~ft41[1:(T1-h-q+1+si-1)]-1)
  bet41=sdpd41$coefficients
  SDPCA1[si]=yy[T1+si-1+1]-ft41[T1-h-q+2+si-1]*bet41
  
  ft42=ft4[,1:2]
  sdpd42=lm(syy~ft42[1:(T1-h-q+1+si-1),1:2]-1)
  bet42=sdpd42$coefficients
  SDPCA2[si]=yy[T1+si-1+1]-t(ft42[T1-h-q+2+si-1,1:2])%*%bet42
  
  ft43=ft4[,1:3]
  sdpd43=lm(syy~ft43[1:(T1-h-q+1+si-1),1:3]-1)
  bet43=sdpd43$coefficients
  SDPCA3[si]=yy[T1+si-1+1]-t(ft43[T1-h-q+2+si-1,1:3])%*%bet43
  
  ft44=ft4[,1:4]
  sdpd44=lm(syy~ft44[1:(T1-h-q+1+si-1),1:4]-1)
  bet44=sdpd44$coefficients
  SDPCA4[si]=yy[T1+si-1+1]-t(ft44[T1-h-q+2+si-1,1:4])%*%bet44
  
  ft45=ft4[,1:5]
  sdpd45=lm(syy~ft45[1:(T1-h-q+1+si-1),1:5]-1)
  bet45=sdpd45$coefficients
  SDPCA5[si]=yy[T1+si-1+1]-t(ft45[T1-h-q+2+si-1,1:5])%*%bet45
  ft46=ft4[,1:6]
  sdpd46=lm(syy~ft46[1:(T1-h-q+1+si-1),1:6]-1)
  bet46=sdpd46$coefficients
  SDPCA6[si]=yy[T1+si-1+1]-t(ft46[T1-h-q+2+si-1,1:6])%*%bet46
  
  md.lasso=lars(ft46[1:(T1-h-q+1+si-1),1:6],syy,type="lasso",trace=FALSE,normalize=TRUE,intercept=FALSE)
  bet47=md.lasso$beta[4,]
  SDPCA7[si]=yy[T1+si-1+1]-t(ft46[T1-h-q+2+si-1,1:6])%*%bet47
  AR_fit <-arima(yy[1:(T1-h+1+si-1)], order  = c(1,0,0))
  arf <- predict(AR_fit, n.ahead = h)$pred
  ARR[si]=yy[T1+si-1+1]-arf[h]
  AR_fit2 <-arima(yy[1:(T1-h+1+si-1)], order  = c(2,0,0))
  arf2 <- predict(AR_fit2, n.ahead = h)$pred
  ARR2[si]=yy[T1+si-1+1]-arf2[h]
  
}
##########1ft
sqrt(mean(sum(SDPCA1^2)))
sqrt(mean(sum(PCA1^2))) 
sqrt(mean(sum(SPCA1^2)))
sqrt(mean(sum(SW1^2)))
####2fts
sqrt(mean(sum(SDPCA2^2)))
sqrt(mean(sum(PCA2^2))) 
sqrt(mean(sum(SPCA2^2)))
sqrt(mean(sum(SW2^2)))
############3fts
sqrt(mean(sum(SDPCA3^2)))
sqrt(mean(sum(PCA3^2))) 
sqrt(mean(sum(SPCA3^2)))
sqrt(mean(sum(SW3^2)))
#######4fts5fts for sdpca
sqrt(mean(sum(SDPCA4^2)))
sqrt(mean(sum(SDPCA5^2)))
sqrt(mean(sum(SDPCA6^2)))
sqrt(mean(sum(SDPCA7^2)))
#######AR1
sqrt(mean(sum(ARR^2)))
sqrt(mean(sum(ARR2^2)))
