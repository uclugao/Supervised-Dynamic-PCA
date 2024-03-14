library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
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
SDPCA8=numeric(T2)#debiased
SDPCA9=numeric(T2)#debiased
SDPCA10=numeric(T2)#debiased
######################
h=1
q=2
XX=fdata
XX=scale(XX, center = F, scale = TRUE)
yy=y5
#yy=fdata[,32]
for (si in 1:(T2)){
  syy=yy[(q+h):(T1+si-1)]#regression
  sXX=XX[1:(T1-h+1+si-1),]
  sXX=scale(sXX, center = F, scale = TRUE)
  PCV=sXX%*%t(sXX)
  ftt=eigen(PCV)$vectors
  
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
  
  md.lasso=lars(sXNN[1:(T1-h-q+1+si-1),],syy,type="lasso",trace=FALSE,normalize=TRUE,intercept=FALSE)
  bet49=md.lasso$beta[5,]#number of nonzero
  SDPCA9[si]=yy[T1+si-1+1]-t(sXNN[T1-h-q+2+si-1,])%*%bet49
  indx=which(md.lasso$beta[5,]!=0)
  sXNN1=sXNN[,indx]
  beta50=lm(syy~sXNN1[1:(T1-h-q+1+si-1),]-1)$coefficients
  SDPCA10[si]=yy[T1+si-1+1]-t(sXNN1[T1-h-q+2+si-1,])%*%beta50
}
##########1ft
sqrt(mean(sum(SDPCA9^2)))
####2fts
sqrt(mean(sum(SDPCA10^2)))
############3fts

