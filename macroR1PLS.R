library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
library(pls)
###########
Mdata=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/macro-data/FREW-MD-md.csv")
data2=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/sp-data/predictor-Monthly-new.csv")
y5=sqrt(data2[55:744,15])-sqrt(data2[54:743,15])#changevolatility
y6=data2[55:744,17]#return
#dim(Mdata)
fdata=Mdata[2:691,2:124]
#dim(fdata)
#ts.plot(fdata[,123])
y1=fdata[,6]#IP
y2=fdata[,24]#UNRATE
y3=fdata[,103]#CPI
y4=fdata[,123]#VIX
y7=fdata[,2]
y8=fdata[,4]#Real Manu. and Trade Industries Sales
y9=fdata[,3]
y10=fdata[,32]#nonfarm
######
TT=690
NN=123
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
PLSO1=numeric(T2)
PLSO2=numeric(T2)
PLSO3=numeric(T2)
######################
h=1#change here
nc=1#components
  XX=fdata
  XX=scale(XX, center = F, scale = TRUE)
  yy=y5#change here
  
  X1=XX[1:(T1-h),]
  yy1=yy[(1+h):(T1)]
  
  X2=XX[(T1-h+1):(TT-h),]
  yy2=yy[(T1+1):(TT)]
  
  mod1 <- plsr(yy1~X1, scale=TRUE, ncomp=3,validation="CV")
  pcr_pred1 <- predict(mod1, X2, ncomp=1)
  PLSO1=yy2-pcr_pred1
  pcr_pred2 <- predict(mod1, X2, ncomp=2)
  PLSO2=yy2-pcr_pred2
  
  pcr_pred3 <- predict(mod1, X2, ncomp=3)
  PLSO3=yy2-pcr_pred3
 
##########1ft
sqrt(mean(sum(PLSO1^2)))
sqrt(mean(sum(PLSO2^2))) 
sqrt(mean(sum(PLSO3^2))) 