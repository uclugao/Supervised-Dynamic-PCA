library(MASS)
library(Matrix)
library(stats)
library(base)
library(MTS)
library(forecast)
library(FactoMineR)
###########
Mdata=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/macro-data/FREW-MD-md.csv")
data2=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/sp-data/predictor-Monthly-new.csv")
gp=read.csv("/Users/jasongao/Documents/Papers/MS-topic/Data/macro-data/FREW-MD-md-gp.csv")
grp=gp[1,2:124]
c1=which(grp==1)
c2=which(grp==2)
c3=which(grp==3)
c4=which(grp==4)
c5=which(grp==5)
c6=which(grp==6)
c7=which(grp==7)
c8=which(grp==8)
length(c1)+length(c2)+length(c3)+length(c4)+length(c5)+length(c6)+length(c7)+length(c8)
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
RR=matrix(numeric(NN*6),NN,6)
#####
h=1
q=5
XX=fdata
XX=scale(XX, center = F, scale = TRUE)
pdt=cbind(y1,y2,y3,y8,y5,y6)
syy=pdt[(q+h):(TT),]
for (i in 1:6){
  for (j in 1:NN){
    xds=NULL
    for (xd in 1:q){
      xds=cbind(xds,XX[(q-xd+1):(TT-h-xd+1),j])
    }
    model=lm(syy[,i]~xds[,1]+xds[,2]+xds[,3]+xds[,4]+xds[,5]-1)
    tmm=stepAIC(model, direction="both",trace=0)
    sm=summary(tmm)
    RR[j,i]=sm$r.squared
  }
}
####test
#model=lm(syy[,i]~xds[,1]+xds[,2]+xds[,3]+xds[,4]+xds[,5]-1)
#tmm=stepAIC(model, direction="both")

#res <- LinearModel(sdmd, selection="AIC")
#add1(lm(syy[,i]~ xds[,1]-1), yy[,i]~ ltakers + income + years + public +
#       expend + rank, test="F")
#modAIC <- MASS::stepAIC(sdmd, direction = 'backward')
#summary(modAIC)$r.squared


###############
nc=c(c1,c2,c3,c4,c5,c6,c7,c8)
lgt=c(length(c1),length(c2),length(c3),length(c4),length(c5),length(c6),length(c7),length(c8))
#ND=numeric(NN)
#ND=c(RR[c1,i],RR[c2,i],RR[c3,i],RR[c4,i],RR[c5,i],RR[c6,i],RR[c7,i],RR[c8,i])
#barplot(ND,xlim=c(0,150),ylim=c(0,1))
RRN=rbind(RR[c1,],RR[c2,],RR[c3,],RR[c4,],RR[c5,],RR[c6,],RR[c7,],RR[c8,])
lb=c("OUT","LM","HS","COI","MC","IER","PR","SM")
GP=c(rep(lb[1],lgt[1]),rep(lb[2],lgt[2]),rep(lb[3],lgt[3]),rep(lb[4],lgt[4]),
     rep(lb[5],lgt[5]),rep(lb[6],lgt[6]),rep(lb[7],lgt[7]),rep(lb[8],lgt[8]))
RRNN=cbind(RRN,GP)
#colnames(RRN)=c(rep('R2',123),'group')
mm=c("IP","UNRATE","CPI-All","M&T Sales",
     "Volatility Change",
     "Return")
####
par(mfrow=c(3,2))
par(mar= c(1.2,4.5,1.5,1))
x <- RRN[,5]
# create labels
x.labels <- paste0(RRNN[,7])
# specify colors for groups
group.cols <- c("darkred", "red", "darksalmon", 
                "darkblue", "blue", "lightblue","blueviolet","chocolate1")
cols <- c(rep(group.cols[1], lgt[1]), rep(group.cols[2],lgt[2]), 
          rep(group.cols[3], lgt[3]), rep(group.cols[4], lgt[4]), 
          rep(group.cols[5], lgt[5]), rep(group.cols[6], lgt[6]),
          rep(group.cols[7], lgt[7]),rep(group.cols[8], lgt[8]))
######R2 by AIC
x <- RRN[,1]
barplot(x, space = 0,  col = cols,ylim=c(0,0.26),
        main=expression(paste("Panel A: IP")),
        ylab=expression(paste(R^2)))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.4)
x <- RRN[,2]
barplot(x, space = 0,  col = cols,ylim=c(0,0.25),
        main=expression(paste("Panel B: UNRATE")),
        ylab=expression(paste(R^2)))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.4)
x <- RRN[,3]
barplot(x, space = 0,  col = cols,ylim=c(0,0.26),
        main=expression(paste("Panel C: CPI-All")),
        ylab=expression(paste(R^2)))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.4)
x <- RRN[,4]
barplot(x, space = 0,  col = cols,ylim=c(0,0.2),
        main=expression(paste("Panel D: M&T Sales")),
        ylab=expression(paste(R^2)))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.4)
x <- RRN[,5]
barplot(x, space = 0,  col = cols,ylim=c(0,0.15),
        main=expression(paste("Panel E: Volatility Change")),
        ylab=expression(paste(R^2)))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.4)
x <- RRN[,6]
barplot(x, space = 0,  col = cols,ylim=c(0,0.1),
        main=expression(paste("Panel F: Return")),
        ylab=expression(paste(R^2)))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.4)
##############loadings
M=t(XX)%*%XX
lds=eigen(M)$vectors[,1:8]
####
d1=eigen(M)$values
d1[1]/sum(d1)
#####
MN=rbind(lds[c1,],lds[c2,],lds[c3,],lds[c4,],
         lds[c5,],lds[c6,],lds[c7,],lds[c8,])
MNN=cbind(MN,GP)
##############plot
par(mfrow=c(3,2))
par(mar= c(1.2,4.5,1.5,1))
x <- as.numeric(MNN[,1])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3),
        main=expression(paste("1st")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.5)
x <- as.numeric(MNN[,2])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,4.5),
        main=expression(paste("2nd")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.5)
x <- as.numeric(MNN[,3])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3.5),
        main=expression(paste("3rd")))
legend("topleft", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.5)
x <- as.numeric(MNN[,4])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3),
        main=expression(paste("4th")))
legend("topleft", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.5)
x <- as.numeric(MNN[,5])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3),
        main=expression(paste("5th")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.5)
x <- as.numeric(MNN[,6])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3.5),
        main=expression(paste("6th")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.5)
###################
h=1
q=1
XX=fdata
XX=scale(XX, center = F, scale = TRUE)
yy=y6
syy=yy[(q+h):(TT)]#regression
sXX=XX[1:(TT-h),]
sXX=scale(sXX, center = F, scale = TRUE)
sXNN=matrix(numeric((TT-h-q+1)*NN),TT-h-q+1,NN)
#sXNN=matrix(numeric((T1-1+si-1+1)*NN),T1+si-1-1+1,NN)
for (ii in 1:NN){
  xds=NULL
  for (xd in 1:q){
    xds=cbind(xds,XX[(q-xd+1):(TT-h+1-xd),ii])
  }
  xds=scale(xds,center=F,scale = TRUE)
  sdmd=lm(syy~xds-1)
  #XNN[,ii]=xds%*%sdmd$coefficients
  # sxds=cbind(XX[2:(T1+si-1+1),ii],XX[1:(T1+si-1-1+1),ii])
  #sxds=scale(sxds,center=TRUE,scale = TRUE)
  sXNN[,ii]=xds%*%sdmd$coefficients
}
sdPCA=t(sXNN)%*%sXNN
ft4=eigen(sdPCA)$vectors[,1:6]
###first eigen exp
d2=eigen(sdPCA)$values
d2[1]/sum(d2)
#####
FT=rbind(ft4[c1,],ft4[c2,],ft4[c3,],ft4[c4,],
         ft4[c5,],ft4[c6,],ft4[c7,],ft4[c8,])
FTT=cbind(FT,GP)
##############plot
par(mfrow=c(6,1))
par(mar= c(1.2,4.5,1.5,1))
x <- as.numeric(FTT[,1])*10
barplot(x, space = 0,  col = cols,ylim=c(-5,3),
        main=expression(paste("Return")),
        ylab=expression(paste("1st")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.3)
x <- as.numeric(FTT[,2])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,4.5),
        ylab=expression(paste("2nd")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.3)
x <- as.numeric(FTT[,3])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,4),
        ylab=expression(paste("3rd")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.3)
x <- as.numeric(FTT[,4])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3),
        ylab=expression(paste("4th")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.3)
x <- as.numeric(FTT[,5])*10
barplot(x, space = 0,  col = cols,ylim=c(-3,3),
        ylab=expression(paste("5th")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.3)
x <- as.numeric(FTT[,6])*10
barplot(x, space = 0,  col = cols,ylim=c(-3.5,3.5),
        ylab=expression(paste("6th")))
legend("topright", legend = c(unique(x.labels)), 
       col = group.cols, pch = 1,cex=0.3)
###################