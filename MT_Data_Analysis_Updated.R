library(readxl)
library(vioplot)
library(matrixStats)
library(stinepack)
library(viridis)
library(fitdistrplus)
library(zoib)
library(gamlss)

setwd("C:/users/jmh22/Documents/GitHub/MASH-Main/MASH-dev/JohnHenry/PDG")
### Extracting Asexual Parasitemia Conditioned on Detection - recorded
### zeroes and non-recordings are both listed as NaN
MT = read.delim('MalariaTherapy.txt')
pID = unique(MT$patient)
UT = (MT$Tretmentfree.==1)
isPf = rep(F,length(MT$strain))
for(i in 1:length(MT$strain)){
  isPf[i] = (MT$strain[i]=="Santee_C" | MT$strain[i]=="ElLimon" | MT$strain[i] == "McLendon")
}

## extract the IDs of those who were untreated at time i
UTd = matrix(NaN,nrow=2000,ncol=length(pID))
for(i in 1:length(pID)){
  UTd[1:length(which(MT$patient==pID[i])),i] = UT[which(MT$patient==pID[i])]
}

## define matrix of all asexual parasitemia (MM) and those untreated (UTM)
MM = matrix(NaN,nrow=2000,ncol=length(pID))
UTM = MM
P = MT$Asexual
P[which(P=='.')]=0
P = as.numeric(as.character(P))
## Remove all non-falciparum infections
P[which(isPf==F)]=NaN
for(i in 1:length(pID)){
  MM[1:length(which(MT$patient==pID[i])),i] = P[which(MT$patient==pID[i])]
  UTM[1:length(which(MT$patient==pID[i])),i] = P[which(MT$patient==pID[i])]*UT[which(MT$patient==pID[i])]
}
MM[which(MM==0)]=NaN
UTM[which(UTM==0)]=NaN


### Repeat for Gametocytes; all gametocyte densities (GG) and untreated (UTG)
GG = matrix(NaN,nrow=2000,ncol=length(pID))
UTG = GG
G = MT$Gametocytes
G[which(G=='.')]=0
G = as.numeric(as.character(G))
## Remove all non-falciparum infections
G[which(isPf==F)]=NaN
for(i in 1:length(pID)){
  GG[1:length(which(MT$patient==pID[i])),i] = G[which(MT$patient==pID[i])]
  UTG[1:length(which(MT$patient==pID[i])),i] = G[which(MT$patient==pID[i])]*UT[which(MT$patient==pID[i])]
}
GG[which(GG==0)]=NaN
UTG[which(UTG==0)]=NaN



hist(max(log10(MM[20,]),na.rm=T)-log10(MM[20,]),breaks=15,freq=F)

mu20 = mean(max(log10(MM[20,]),na.rm=T)-log10(MM[20,]),na.rm=T)
var20 = var(max(log10(MM[20,]),na.rm=T)-log10(MM[20,]),na.rm=T)
a0 = mu20^2/var20
b0 = mu20/var20
lines(seq(0,5,.01),dgamma(seq(0,5,.01),a0,b0))


### Log of Means vs Mean of Logs
### Log of Means
plot(log10(rowMeans(UTM,na.rm=T)[1:300]),type="l",ylim=c(0,5))
lines(log10(rowMeans(UTG,na.rm=T)[1:300]),type="l",col="red")

plot(log10(rowMeans(MM,na.rm=T)[1:300]),type="l")
lines(log10(rowMeans(GG,na.rm=T)[1:300]),type="l",col="red")
##### MV power law
plot(log10(rowMeans(UTM,na.rm=T)),log10(rowVars(MM,na.rm=T)))

### Mean of Logs
plot(rowMeans(log10(UTM),na.rm=T)[1:200],type="l",ylim=c(1,4.2))
lines(rowMeans(log10(UTG),na.rm=T)[1:200],type="l",col="red")
##### MV power law
plot(rowMeans(log10(UTM),na.rm=T),rowVars(log10(MM),na.rm=T))


### Cross correlation between untreated mean asexual parasites and gametocytes
ccf(rowMeans(log10(UTM),na.rm=T)[5:100],rowMeans(log10(UTG),na.rm=T)[5:100],ylab="CCF",main="Cross-Correlation between Mean log10 Asexual Parasites and Gametocytes")



### Duration of Infection -
### include only untreated infections,
### determine the final day of patency
dur = c()
for(i in 1:length(pID)){
  UTdur = (max(which(MM[,i]>0))==max(which(UTM[,i]>0)))
  if(UTdur==1){
    dur = c(dur,max(which(MM[,i]>0)))
  }
}
dur[which(is.infinite(dur))]=NaN
dur = as.numeric(na.exclude(dur))
hist(dur,freq=F,ylim=c(0,.007),xlab="Days",ylab="Density",main="Duration of Untreated Infections")
expfit = fitdist(dur,distr="exp")
gammafit = fitdist(dur,distr="gamma")
lines(seq(0,500),dgamma(seq(0,500),shape=2.058,rate=.0158))
abline(v=mean(dur),col="red")
weibullfit = fitdist(dur,distr="weibull")
lines(seq(0,500),dweibull(seq(0,500),shape=1.51,scale=145.24),col="blue")
lognormalfit = fitdist(dur,distr="lnorm")
lines(seq(0,500),dlnorm(seq(0,500),4.61,.78),col="green")
legend(250,.007,legend=c("Gamma","Weibull","Lognormal"),lty=1,col=c("black","blue","green"))
### Note - aic and bic are both ambivalent to gamma and weibull distributions, and both show
### the two are superior to lognormal

### approximate survival curve, show fitted gamma
surv = rep(0,300)
for(i in 1:300){
  surv[i] = 1-sum(dur<=i)/length(dur)
}
days = seq(1,300)
plot(days,surv,xlab="Days",ylab="Fraction with Active Infection",main="Survival Curve of Infections")
lines(1-pgamma(days,shape=2.058,rate=.0158),lwd=3,col="blue")



### Patent Fraction
##### Patent fraction for untreated infections given continued infection
UTID = rep(0,length(pID))
for(i in 1:length(pID)){
  UTID[i] = (max(which(MM[,i]>0))==max(which(UTM[,i]>0)))
}
pat = rep(0,300)
for(i in 1:300){
  pat[i] = sum(UTM[i,which(UTID==1)]>0,na.rm=T)/sum(dur>i)
}
plot(days,pat,type="l")


plot(rowSums(UTM[,which(UTID==1)]>0,na.rm=T)[1:300]/pat[1:300],type="l",xlab="Day",ylab="Fraction Subpatent",main="Fraction of Untreated Infections Subpatent")
##### Patent fraction for all infections given continued infection
durall = rep(0,length(pID))
for(i in 1:length(pID)){
  durall[i] = max(which(MM[,i]>0))
}
patall = rep(0,300)
for(i in 1:300){
  patall[i] = sum(MM[i,]>0,na.rm=T)/sum(durall>i,na.rm=T)
}
plot(days,patall,type="l",xlab="Days",ylab="Fraction of Active Infections",main="Patent Fraction Conditioned on Active Infection")
lines(days,pat,type="l",col="blue")
legend(x=200,y=.9,legend=c("All Infections","Untreated Only"),lty=1,col=c("black","blue"))

plot(rowSums(MM>0,na.rm=T)[1:300]/patall[1:300],type="l",xlab="Day",ylab="Fraction Patent",main="Fraction of All Infections Patent")
lines(rowSums(UTM[,which(UTID==1)]>0,na.rm=T)[1:300]/pat[1:300],col="red")

### Fever - look at fraction with a fever
Temperature = matrix(NaN,nrow=2000,ncol=length(pID))
Temp = MT$Temp
Temp[which(Temp=='.')]=NaN
Temp = as.numeric(as.character(Temp))
for(i in 1:length(pID)){
  Temperature[1:length(which(MT$patient==pID[i])),i] = Temp[which(MT$patient==pID[i])]
}
Fever = (Temperature>38)

plot(rowMeans(Fever,na.rm=T)[1:200])

plot(log10(rowMeans(UTM,na.rm=T)[1:150]),rowMeans(Fever,na.rm=T)[1:150],ylim=c(0,1),xlab="Mean Log10 Asexual Parasite Densities",ylab="Fraction with Fever",main="Fever Fraction as a function of Parasitemia")
points(log10(rowMeans(UTM,na.rm=T)[1:5]),rowMeans(Fever,na.rm=T)[1:5],col="red")

aa = log10(rowMeans(UTM,na.rm=T)[3:150])
bb = rowMeans(Fever,na.rm=T)[3:150]
sigfitfev = nls(bb~p1*exp(p2*aa)/(p3+exp(p2*aa)),start=list(p1=.9,p2=1,p3=100))
p1 = .8968
p2 = 4.621
p3 = 2.447*10^6
sigmoidFev = function(x,p1,p2,p3){
  p1*exp(p2*x)/(p3+exp(p2*x))
}
lines(seq(1,6,.01),sigmoidFev(seq(1,6,.01),p1,p2,p3))
abline(h=p1,lty=2)

plot(rowMeans(Fever,na.rm=T)[1:200],ylim=c(0,1))
lines(1:180,sigmoidFev(log10(rowMeans(UTM,na.rm=T))[1:180],p1,p2,p3),col="red")



### Transmission Efficiency
Mosquito_Transmission = MT$Mosq.
Mosquito_Transmission[which(Mosquito_Transmission==".")]=NaN
Mosq = as.numeric(as.character(Mosquito_Transmission))
TE = matrix(NaN,2000,length(pID))
UTTE = TE

for(i in 1:length(pID)){
  TE[1:length(which(MT$patient==pID[i])),i] = Mosq[which(MT$patient==pID[i])]
  UTTE[1:length(which(MT$patient==pID[i])),i] = Mosq[which(MT$patient==pID[i])]*UT[which(MT$patient==pID[i])]
}

plot(rowMeans(TE,na.rm=T)[1:150],type="l")
plot(rowMeans(log10(UTG),na.rm=T)[1:150],rowMeans(TE,na.rm=T)[1:150])

plot(log10(G),Mosq,xlim=c(.5,4.5),xlab="Log10 Gametocytemia",ylab="Transmission Efficiency")

### restrict those who were treated
GUT = as.numeric(UT)*G
GUT[which(GUT==0)]=NaN
TEfit = function(x){
  temp = rep(0,length(x)-1)
  weight = temp
  for(i in 1:length(x)){
    temp[i] = mean(Mosq[which(log10(GUT)>=(x[i]-.5) & log10(GUT)<(x[i]+.5))],na.rm=T,weight=dnorm(log10(GUT)[which(log10(GUT)>=(x[i]-.5) & log10(GUT)<(x[i]+.5))]-x[i],x[i],.25))
    weight[i] = length(Mosq[which(log10(GUT)>=(x[i]-.5) & log10(GUT)<(x[i]+.5))])
  }
  output = list(TE = temp/100, weight = weight)
  return(output)
}

TEfitAll = function(x){
  temp = rep(0,length(x)-1)
  weight = temp
  for(i in 1:length(x)){
    temp[i] = mean(Mosq[which(log10(GG)>=(x[i]-.5) & log10(GG)<(x[i]+.5))],na.rm=T,weight=dnorm(log10(GG)[which(log10(GG)>=(x[i]-.5) & log10(GG)<(x[i]+.5))]-x[i],x[i],.25))
    weight[i] = length(Mosq[which(log10(GG)>=(x[i]-.5) & log10(GG)<(x[i]+.5))])
  }
  output = list(TE = temp/100, weight = weight)
  return(output)
}

TEFitAll(x)

TEfit = function(x){
  temp = rep(0,length(x)-1)
  weight = temp
  for(i in 1:length(x)){
    temp[i] = mean(Mosq[which(log10(GUT)>=(x[i]-.5) & log10(GUT)<(x[i]+.5))],na.rm=T,weight=dnorm(log10(GUT)[which(log10(GUT)>=(x[i]-.5) & log10(GUT)<(x[i]+.5))]-x[i],x[i],.25))
    weight[i] = length(Mosq[which(log10(GUT)>=(x[i]-.5) & log10(GUT)<(x[i]+.5))])
  }
  output = list(TE = temp/100, weight = weight)
  return(output)
}

x = seq(1,4,.1)

SmoothTE = TEfit(x)
plot(x,SmoothTE$TE,ylim=c(0,1))

sigfitTE = nls(SmoothTE$TE~p1*exp(p2*x)/(p3+exp(p2*x)),start=list(p1=.6,p2=2.5,p3=100),weights=SmoothTE$weight)
sigfitTE
p1TE = .6938
p2TE = 2.1232
p3TE = 83.9391
sigTE = function(x,p1TE,p2TE,p3TE){
  p1TE*exp(p2TE*x)/(p3TE+exp(p2TE*x))
}
plot(x,SmoothTE$TE,xlab="Smoothed Log10 Gametocytemia",ylab="Transmission Efficiency",ylim=c(0,1))
lines(x,sigTE(x,p1TE,p2TE,p3TE),col="green")

plot(log10(GUT),Mosq/100,xlim=c(.5,4.5),xlab="Log10 Gametocytemia",ylab="Transmission Efficiency",cex=.5,main="Transmission Efficiency Predicted by Gametocytemia")
points(x,SmoothTE$TE,xlab="Smoothed Log10 Gametocytemia",ylab="Transmission Efficiency",pch=19,col="blue")
lines(x,sigTE(x,p1TE,p2TE,p3TE),col="dark green",lwd=5)
abline(h=p1TE,lty=2,col="dark green",lwd=3)

plot(rowMeans(TE,na.rm=T)[1:250]/100)
lines(sigTE(log10(rowMeans(UTG,na.rm=T))[1:250],p1TE,p2TE,p3TE),col="dark green")



###### Organizing Figures


### Duration of Untreated Infections

#par(mfrow=c(1,1))
#hist(dur,freq=F,ylim=c(0,.007),xlab="Days",ylab="Density",main="Duration of Untreated Infections")
#gammafit = fitdist(dur,distr="gamma")
#lines(seq(0,500),dgamma(seq(0,500),shape=2.058,rate=.0158))
#abline(v=mean(dur),col="red")
#weibullfit = fitdist(dur,distr="weibull")
#lines(seq(0,500),dweibull(seq(0,500),shape=1.51,scale=145.24),col="blue")
#lognormalfit = fitdist(dur,distr="lnorm")
#lines(seq(0,500),dlnorm(seq(0,500),4.61,.78),col="green")
#legend(250,.007,legend=c("Gamma","Weibull","Lognormal"),lty=1,col=c("black","blue","green"))

Ptmu = function(a){
  temp = rep(0,length(a))
  temp[which(a<1)] = rep(0,length(which(a<1)))
  temp[which(a>=1 & a < 6)] = 3.1004+.2782*a[which(a>=1 & a < 6)]
  temp[which(a>=6 & a < 18)] = 5.11986-.07425*a[which(a>=6 & a < 18)]
  temp[which(a>=18)] = 3.848168 - .008427*a[which(a>=18)]
  return(temp)
}

Gtmu = function(Pt){
  temp = -.5829+.8579*Pt
  temp = c(rep(NaN,9),temp)
  temp = temp[1:length(Pt)]
  return(temp)
}

## refit patency

plot(17:250,log(pat[17:250]/(1-pat[17:250])))
z = log(pat[17:195]/(1-pat[17:195]))
zt = 17:195
patlm = lm(z~zt)
lines(1:250,1.5211-.0151*1:250)
patlm$coefficients

plot(log(pat[1:250]/(1-pat[1:250])))
lines(1:250,1.5211-.0151*1:250)

##### figure 1
par(mar=c(4,4,2,2))
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
#sub a
plot(log10(rowMeans(MM,na.rm=T)[1:300]),ylim=c(1,4.8),xlim=c(0,250),xlab="Days",ylab="Log10 of the Mean Parasitemia",main="Log10 Average Parasitemia by Day")
rect(par("usr")[1]+16, par("usr")[3], par("usr")[2]-242, par("usr")[4], col = 
       "light gray")
lines(days,Ptmu(days),lwd=3,col="blue")
abline(v=c(6,18),lty=2)
points(log10(rowMeans(MM,na.rm=T)[1:300]))
mtext("A",adj=-.05)
#points(log10(rowMeans(UTG,na.rm=T)[1:300]),col="red")
#sub b
plot(days,surv,xlab="Days",ylab="Proportion",main="Proportion with Persisting Infection",xlim=c(0,250))
lines(1-pgamma(days,shape=2.058,rate=.0158),lwd=3,col="blue")
mtext("B",adj=-.05)
#sub c
plot(days,pat,xlab="Days",ylab="Proportion",main="Proportion Patent, Conditioned on Active Infection",xlim=c(0,250))
rect(par("usr")[1]+16, par("usr")[3], par("usr")[2]-242, par("usr")[4], col = 
       "light gray")
mtext("C",adj=-.05)
points(days,pat)
#lines(days,pat,type="l",col="blue")
#legend(x=200,y=.9,legend=c("All Infections","Untreated Only"),lty=1,col=c("black","blue"))

patf = function(a){
  a00 = which(a >= 0 & a < 6)
  a0 = which(a >= 6 & a <= 17)
  a1 = which(a > 17)
  x0 = a[a00]
  x0 = a[a0]
  x1 = a[a1]
  pat = 0*a
  pat[a00] = 1
  pat[a0] = 1-.02*(x0-6)
  pat[a1] = 1/(1+exp(-1.52+0.0151*x1))
  return(pat)
}
lines(seq(0,250),patf(seq(0,250)),col="blue",lwd=3)

##### Figure 2

vars = rep(0,15)
for(i in 1:15){
  logp = log10(rowMeans(MM,na.rm=T)[5:100])
  logg = log10(rowMeans(GG,na.rm=T)[(4+i):(100+i-1)])
  lmpg = lm(logg~logp)
  vars[i] = var(lmpg$residuals)
}
plot(1:15,sqrt(vars),xlab="Lag (in Days)", ylab="Sum of Squared Residuals",main="Standard Deviation of Residuals of Lagged Regression")

logps = log10(rowMeans(MM,na.rm=T)[5:100])
loggs = log10(rowMeans(GG,na.rm=T)[(4+9):(100+9-1)])
lmpg=lm(loggs~logps)

layout(matrix(c(1,2,1,3),2,2,byrow=T))
plot(days,Ptmu(days),ylim=c(.5,4.8),xlab="Days",ylab="Log10 Mean Parasitemia",type="l",lwd=3,col="blue",main="Log10 Average Parasitemia by Day")
mtext("A",adj=-.05)
legend(160,4.6,col=c("black","green"),legend=c("Asexuals","Gametocytes"),pch=c(1,1))
lines(days,Gtmu(Ptmu(days)),col="blue",lwd=3)
points(days,log10(rowMeans(MM,na.rm=T))[1:300],col="black")
points(days,log10(rowMeans(GG,na.rm=T))[1:300],col="dark green")
plot(log10(rowMeans(MM,na.rm=T))[1:200],log10(rowMeans(GG,na.rm=T))[10:209],xlab="Lagged log10 Average Asexual Parasitemia",ylab="log10 Average Gametocytemia",main="Lagged Mean Asexual-Gametocyte Relationship")
mtext("B",adj=-.1)
lines(seq(1,5,.1),lmpg$coefficients[1]+lmpg$coefficients[2]*seq(1,5,.1),lwd=3,col="blue")
plot(sqrt(vars),xlab="Lag (in Days)",ylab="Standard Deviation of Residuals",main="Residuals of Lagged Regression")
mtext("C",adj=-.1)

##### Figure 3

par(mfrow=c(2,2))
vioplot(log10(colMeans(MM[1:30,],na.rm=T)),log10(colMeans(MM[31:60,],na.rm=T)),log10(colMeans(MM[61:90,],na.rm=T)),log10(colMeans(MM[91:120,],na.rm=T)),log10(colMeans(MM[121:150,],na.rm=T)),log10(colMeans(MM[151:180,],na.rm=T)),log10(colMeans(MM[181:210,],na.rm=T)),xlab="Months Since First Patent",ylab="Log10 Parasitemia per mL",main="Monthly Violin Plot of Asexual Parasitemia",ylim=c(0,5))
mtext("A",adj=-.15)
plot(log10(rowMeans(MM,na.rm=T)),log10(rowVars(MM,na.rm=T)),xlab="log10 Daily Mean",ylab="log10 Daily Variance",main="Asexual Parasitemia Mean-Variance Power Law")
mmmu = log10(rowMeans(MM,na.rm=T))[1:100]
mmvar = log10(rowVars(MM,na.rm=T))[1:100]
mmlm = lm(mmvar~mmmu)
lines(seq(1,6,.1),mmlm$coefficients[1]+mmlm$coefficients[2]*seq(1,6,.1),lwd=3,col="blue")
mtext("B",adj=-.2)
vioplot(log10(colMeans(GG[1:30,],na.rm=T)),log10(colMeans(GG[31:60,],na.rm=T)),log10(colMeans(GG[61:90,],na.rm=T)),log10(colMeans(GG[91:120,],na.rm=T)),log10(colMeans(GG[121:150,],na.rm=T)),log10(colMeans(GG[151:180,],na.rm=T)),log10(colMeans(GG[181:210,],na.rm=T)),xlab="Months Since First Patent",ylab="Log10 Gametocytemia per mL",main="Monthly Violin Plot of Gametocytemia",ylim=c(0,5))
mtext("C",adj=-.15)
plot(log10(rowMeans(GG,na.rm=T)),log10(rowVars(GG,na.rm=T)),xlab="log10 Daily Mean",ylab="log10 Daily Variance",main="Gametocyte Mean-Variance Power Law")
ggmu = log10(rowMeans(GG,na.rm=T))[5:100]
ggvar = log10(rowVars(GG,na.rm=T))[5:100]
gglm = lm(ggvar~ggmu)
lines(seq(0,6,.1),gglm$coefficients[1]+gglm$coefficients[2]*seq(0,6,.1),lwd=3,col="blue")
mtext("D",adj=-.2)


##### Figure 4


par(mfrow=c(2,1))
plot(log10(rowMeans(MM,na.rm=T)[6:150]),rowMeans(Fever,na.rm=T)[6:150],ylim=c(0,1),xlab="Mean Log10 Asexual Parasite Densities",ylab="Fraction with Fever",main="Fever Fraction as a Function of Parasitemia",cex=.75)
points(log10(rowMeans(MM,na.rm=T)[1:5]),rowMeans(Fever,na.rm=T)[1:5],col="purple",cex=.75,pch=13)
lines(seq(1,6,.01),sigmoidFev(seq(1,6,.01),p1,p2,p3),col="red",lwd=3)
#abline(h=p1,lty=2)
mtext("A",adj=-.05)
plot(rowMeans(Fever,na.rm=T)[1:200],ylim=c(0,1),xlab="Days Since First Patent",ylab="Probability of Fever",cex=.75,main="Predicted Average Probability of Fever")
lines(1:200,sigmoidFev(Ptmu(1:200),p1,p2,p3),col="red",lwd=3)
mtext("B",adj=-.05)


##### Figure 5

TE1 = Mosq[which(log10(G)>=1 & log10(G)<1.5 & !is.na(Mosq))]/100
TE2 = Mosq[which(log10(G)>=1.5 & log10(G)<2 & !is.na(Mosq))]/100
TE3 = Mosq[which(log10(G)>=2 & log10(G)<2.5 & !is.na(Mosq))]/100
TE4 = Mosq[which(log10(G)>=2.5 & log10(G)<3 & !is.na(Mosq))]/100
TE5 = Mosq[which(log10(G)>=3 & log10(G)<3.5 & !is.na(Mosq))]/100
TE6 = Mosq[which(log10(G)>=3.5 & log10(G)<4 & !is.na(Mosq))]/100
TE1fitZI = fitdist(TE1[which(TE1>0)],distr="beta")
TE2fitZI = fitdist(TE2[which(TE2>0)],distr="beta","mme")
TE3fitZI = fitdist(TE3[which(TE3>0)],distr="beta","mme")
TE4fitZI = fitdist(TE4[which(TE4>0)],distr="beta","mme")
TE5fitZI = fitdist(TE5[which(TE5>0)],distr="beta","mme")
TE6fitZI = fitdist(TE6[which(TE6>0)],distr="beta","mme")


layout(matrix(c(1,1,2,2,3,3,4,5,6,7,8,9),ncol=2,byrow=F))
par(mar = c(5,4,2,2))

plot(log10(G),Mosq/100,xlim=c(.5,4.5),xlab="Log10 Gametocytemia",ylab="Transmission Efficiency",cex=.5,main="Transmission Efficiency from Gametocytemia")
mtext("A",adj=-.07)
points(x,SmoothTE$TE,xlab="Smoothed Log10 Gametocytemia",ylab="Transmission Efficiency",pch=19,col="blue")
lines(x,sigTE(x,p1TE,p2TE,p3TE),col="dark green",lwd=3)

plot(rowMeans(TE,na.rm=T)[1:250]/100,xlab="Days",ylab="Transmission Efficiency",main="Transmission Efficiency over Time")
mtext("B",adj=-.05)
lines(sigTE(Gtmu(Ptmu(1:250)),p1TE,p2TE,p3TE),col="dark green",lwd=2)

alphas = c(TE1fitZI$estimate[1],TE2fitZI$estimate[1],TE3fitZI$estimate[1],TE4fitZI$estimate[1],TE5fitZI$estimate[1],TE6fitZI$estimate[1])
betas = c(TE1fitZI$estimate[2],TE2fitZI$estimate[2],TE3fitZI$estimate[2],TE4fitZI$estimate[2],TE5fitZI$estimate[2],TE6fitZI$estimate[2])
ZI = c(sum(TE1==0)/length(TE1),sum(TE2==0)/length(TE2),sum(TE3==0)/length(TE3),sum(TE4==0)/length(TE4),sum(TE5==0)/length(TE5),sum(TE6==0)/length(TE6))
xx = c(1.25,1.75,2.25,2.75,3.25,3.75)
sigfitZI = nls(ZI~p1*exp(p2*xx)/(p3+exp(p2*xx)),start=list(p1=.7,p2=-2,p3=.1))
plot(xx,ZI,xlab="Log10 Gametocytemia",ylab="Zero-Inflation in TE",ylim=c(0,.8),main="Degree of Zero Inflation")
mtext("C",adj=-.05)
tt = seq(1,4,.1)
lines(tt,sigmoidFev(tt,.8,-1.6,.029),lwd=3,col="blue")

par(mar = c(4,4,1.1,8))
hist(TE1[which(TE1>0)],freq=F,breaks=20,xlab="Transmission Efficiency",ylim=c(0,8),xlim=c(0,1),main="")
legend(.2, 10, "1 < [Gt] <= 1.5", bty = "n",cex=1.3)
mtext("D",adj=-.05)
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE1fitZI$estimate[1],TE1fitZI$estimate[2]),lwd=3,col="dark green")
hist(TE2[which(TE2>0)],freq=F,breaks=20,xlab="Transmission Efficiency",ylim=c(0,8),main="")
legend(.2, 10, "1.5 < [Gt] <= 2", bty = "n",cex=1.3)
mtext("E",adj=-.05)
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE2fitZI$estimate[1],TE2fitZI$estimate[2]),lwd=3,col="dark green")
hist(TE3[which(TE3>0)],freq=F,breaks=20,xlab="Transmission Efficiency",ylim=c(0,8),main="")
legend(.2, 10, "2 < [Gt] <= 2.5", bty = "n",cex=1.3)
mtext("F",adj=-.05)
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE3fitZI$estimate[1],TE3fitZI$estimate[2]),lwd=3,col="dark green")
hist(TE4[which(TE4>0)],freq=F,breaks=20,xlab="Transmission Efficiency",ylim=c(0,8),main="")
legend(.2, 10, "2.5 < [Gt] <= 3", bty = "n",cex=1.3)
mtext("G",adj=-.05)
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE4fitZI$estimate[1],TE4fitZI$estimate[2]),lwd=3,col="dark green")
hist(TE5[which(TE5>0)],freq=F,breaks=20,xlab="Transmission Efficiency",ylim=c(0,8),main="")
legend(.2, 10, "3 < [Gt] <= 3.5", bty = "n",cex=1.3)
mtext("H",adj=-.05)
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE5fitZI$estimate[1],TE5fitZI$estimate[2]),lwd=3,col="dark green")
hist(TE6[which(TE6>0)],freq=F,breaks=20,xlab="Transmission Efficiency",ylim=c(0,8),main="") 
legend(.2, 10, "3.5 < [Gt]", bty = "n",cex=1.3)
mtext("I",adj=-.05)
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE6fitZI$estimate[1],TE6fitZI$estimate[2]),lwd=3,col="dark green")


### Patent Fraction

par(mfrow=c(1,2))
plot(1-rowSums(UTM[,which(UTID==1)]>0,na.rm=T)[1:300]/pat[1:300],type="l",xlab="Day",ylab="Fraction Subpatent",main="Fraction of Untreated Infections Subpatent")
plot(1-rowSums(MM>0,na.rm=T)[1:300]/patall[1:300],type="l",xlab="Day",ylab="Fraction Subpatent",main="Fraction of All Infections Subpatent")


### Average Parasitemia and Gametocytemia

### All Data
### Mean of Logs
layout(rbind(c(1,2),c(1,3)))
plot(rowMeans(log10(MM),na.rm=T)[1:250],type="l",ylim=c(1,4.2),xlab="Days Since First Patent",ylab="Log10 Parasitemia",main="Averaged Log10 Parasitemia Across All Infections")
lines(rowMeans(log10(GG),na.rm=T)[1:250],type="l",col="dark green")
legend(70,y=3.9,legend=c("Asexual Density","Gametocyte Density"),col=c("black","dark green"),lty=c(1,1))
##### MV power law - mean of log
plot(rowMeans(log10(MM),na.rm=T)[1:150],rowVars(log10(MM),na.rm=T)[1:150],xlim=c(1,4.5),ylim=c(0,1.3),xlab="Means of log10 Parasitemia",ylab="Variances of log10 Parasitemia")
mvpl = lm(rowVars(log10(MM),na.rm=T)[11:150]~rowMeans(log10(MM),na.rm=T)[11:150])
lines(seq(1,5,.01),seq(1,5,.01)*mvpl$coefficients[2]+mvpl$coefficients[1])
points(rowMeans(log10(MM),na.rm=T)[5:10],rowVars(log10(MM),na.rm=T)[5:10],pch=8)
mvplg = lm(rowVars(log10(GG),na.rm=T)[11:150]~rowMeans(log10(GG),na.rm=T)[11:100])
lines(seq(1,5,.01),seq(1,5,.01)*mvplg$coefficients[2]+mvplg$coefficients[1],col="dark green")
points(rowMeans(log10(GG),na.rm=T)[1:150],rowVars(log10(GG),na.rm=T)[1:150],col="dark green")
vars = rep(0,15)
for(i in 1:15){
  logp = rowMeans(log10(MM),na.rm=T)[1:150][which(!is.na(rowMeans(log10(GG),na.rm=T)[i:(150+i)]))]
  logg = rowMeans(log10(GG),na.rm=T)[i:(150+i)][which(!is.na(rowMeans(log10(GG),na.rm=T)[i:(150+i)]))]
  lmpg = lm(logg~logp)
  vars[i] = var(lmpg$residuals)
}
plot(1:15,vars,xlab="Lag (in Days)", ylab="Sum of Squared Residuals",main="Variance of Residuals of Lagged Regression")


### Log of Means
layout(rbind(c(1,2),c(1,3)))
plot(log10(rowMeans(MM,na.rm=T))[1:300],type="l",ylim=c(.5,5),xlab="Days Since First Patent",ylab="Log10 Parasitemia",main="Log10 of Average Parasitemia Across All Infections")
lines(log10(rowMeans(GG,na.rm=T))[1:300],type="l",col="dark green")
##### MV power law - mean of log
plot(log10(rowMeans(MM,na.rm=T))[1:150],log10(rowVars(MM,na.rm=T)[1:150]),xlim=c(1.5,4.7),ylim=c(2,10),xlab="log10 of Average Parasitemia",ylab="log10 of Variance of Parasitemia",main="Mean-Variance Power Law")
points(log10(rowMeans(GG,na.rm=T))[1:150],log10(rowVars(GG,na.rm=T)[1:150]),col="dark green")
vars = rep(0,15)
for(i in 1:15){
  logp = log10(rowMeans(MM,na.rm=T)[5:200])
  logg = log10(rowMeans(GG,na.rm=T)[(4+i):(199+i)])
  lmpg = lm(logg~logp)
  vars[i] = var(lmpg$residuals)
}
plot(1:15,vars,xlab="Lag (in Days)", ylab="Sum of Squared Residuals",main="Variance of Residuals of Lagged Regression")


### Untreated only
layout(rbind(c(1,2),c(1,3)))
plot(rowMeans(log10(UTM),na.rm=T)[1:200],type="l",ylim=c(1,4.2),xlab="Days Since First Patent",ylab="Log10 Parasitemia",main="Averaged Log10 Parasitemia Across Untreated Infections")
lines(rowMeans(log10(UTG),na.rm=T)[1:200],type="l",col="red")
##### MV power law - mean of log
#plot(rowMeans(log10(MM),na.rm=T),rowVars(log10(UTM),na.rm=T))
#plot(rowMeans(log10(GG),na.rm=T),rowVars(log10(UTG),na.rm=T))
##### MV power law - log of mean
plot(log10(rowMeans(UTM,na.rm=T)),log10(rowVars(UTM,na.rm=T)),main="log-log Plot of Mean vs Variance of Parasitemia",xlab="log10 Mean of Asexual Parasitemia",ylab="log10 Variance of Asexual Parasitemia")
points(log10(rowMeans(UTG,na.rm=T)),log10(rowVars(UTG,na.rm=T)),col="dark green")

ccf(rowMeans(log10(UTM),na.rm=T)[5:100],rowMeans(log10(UTG),na.rm=T)[5:100],main="CCF Between Average log10 Parasitemia and Gametocytemia",ylab="CCF")



### Fever Probability and TE

par(mfrow=c(2,2))

plot(log10(rowMeans(UTM,na.rm=T))[1:150],rowMeans(Fever,na.rm=T)[1:150],ylim=c(0,1),xlab="Mean Log10 Asexual Parasite Densities",ylab="Fraction with Fever",main="Fever Fraction as a Function of Parasitemia",cex=.5)
points(log10(rowMeans(UTM,na.rm=T))[1:5],rowMeans(Fever,na.rm=T)[1:5],col="red")
lines(seq(1,6,.01),sigmoidFev(seq(1,6,.01),p1,p2,p3),col="red",lwd=3)
abline(h=p1,lty=2)

plot(log10(G),Mosq/100,xlim=c(.5,4.5),xlab="Log10 Gametocytemia",ylab="Transmission Efficiency",cex=.5,main="Transmission Efficiency Predicted by Gametocytemia")
points(x,SmoothTE$TE,xlab="Smoothed Log10 Gametocytemia",ylab="Transmission Efficiency",pch=19,col="blue")
lines(x,sigTE(x,p1TE,p2TE,p3TE),col="dark green",lwd=3)
abline(h=p1TE,lty=2,col="dark green",lwd=3)

plot(rowMeans(Fever,na.rm=T)[1:200],ylim=c(0,1),xlab="Days Since First Patent",ylab="Probability of Fever",cex=.5,main="Predicted Average Probability of Fever")
lines(1:180,sigmoidFev(log10(rowMeans(UTM,na.rm=T))[1:180],p1,p2,p3),col="red")

plot(rowMeans(TE,na.rm=T)[1:250]/100,cex=.5,xlab="Days Since First Patent",ylab="Transmission Efficiency",main="Predicted Average Transmission Efficiency")
lines(sigTE(rowMeans(log10(GG),na.rm=T)[1:250],p1TE,p2TE,p3TE),col="dark green")


### Beta Binomial TE fits

## all data
par(mfrow=c(2,3))
TE1 = Mosq[which(log10(G)>=1 & log10(G)<1.5 & !is.na(Mosq))]/100
TE2 = Mosq[which(log10(G)>=1.5 & log10(G)<2 & !is.na(Mosq))]/100
TE3 = Mosq[which(log10(G)>=2 & log10(G)<2.5 & !is.na(Mosq))]/100
TE4 = Mosq[which(log10(G)>=2.5 & log10(G)<3 & !is.na(Mosq))]/100
TE5 = Mosq[which(log10(G)>=3 & log10(G)<3.5 & !is.na(Mosq))]/100
TE6 = Mosq[which(log10(G)>=3.5 & log10(G)<4 & !is.na(Mosq))]/100
TE1fit = fitdist(TE1,distr="beta","mme")
hist(TE1,freq=F,breaks=20,xlab="Transmission Efficiency",main="1 <= log10[Gametocyte] < 1.5")
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE1fit$estimate[1],TE1fit$estimate[2]))
TE2fit = fitdist(TE2,distr="beta","mme")
hist(TE2,freq=F,breaks=20,xlab="Transmission Efficiency",main="1.5 <= log10[Gametocyte] < 2") 
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE2fit$estimate[1],TE2fit$estimate[2]))
TE3fit = fitdist(TE3,distr="beta","mme")
hist(TE3,freq=F,breaks=20,xlab="Transmission Efficiency",main="2 <= log10[Gametocyte] < 2.5")
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE3fit$estimate[1],TE3fit$estimate[2]))
TE4fit = fitdist(TE4,distr="beta","mme")
hist(TE4,freq=F,breaks=20,xlab="Transmission Efficiency",main="2.5 <= log10[Gametocyte] < 3")
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE4fit$estimate[1],TE4fit$estimate[2]))
TE5fit = fitdist(TE5,distr="beta","mme")
hist(TE5,freq=F,breaks=20,xlab="Transmission Efficiency",main="3 <= log10[Gametocyte] < 3.5")
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE5fit$estimate[1],TE5fit$estimate[2]))
TE6fit = fitdist(TE6,distr="beta","mme")
hist(TE6,freq=F,breaks=20,xlab="Transmission Efficiency",main="3.5 <= log10[Gametocyte] < 4") 
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE6fit$estimate[1],TE6fit$estimate[2]))

## zeros modeled separately

TE1 = Mosq[which(log10(G)>=1 & log10(G)<1.5 & !is.na(Mosq))]/100
TE2 = Mosq[which(log10(G)>=1.5 & log10(G)<2 & !is.na(Mosq))]/100
TE3 = Mosq[which(log10(G)>=2 & log10(G)<2.5 & !is.na(Mosq))]/100
TE4 = Mosq[which(log10(G)>=2.5 & log10(G)<3 & !is.na(Mosq))]/100
TE5 = Mosq[which(log10(G)>=3 & log10(G)<3.5 & !is.na(Mosq))]/100
TE6 = Mosq[which(log10(G)>=3.5 & log10(G)<4 & !is.na(Mosq))]/100
TE1fitZI = fitdist(TE1[which(TE1>0)],distr="beta")
TE2fitZI = fitdist(TE2[which(TE2>0)],distr="beta","mme")
TE3fitZI = fitdist(TE3[which(TE3>0)],distr="beta","mme")
TE4fitZI = fitdist(TE4[which(TE4>0)],distr="beta","mme")
TE5fitZI = fitdist(TE5[which(TE5>0)],distr="beta","mme")
TE6fitZI = fitdist(TE6[which(TE6>0)],distr="beta","mme")

GtTE = t(matrix(c(1,2,3,4,5,6,7,7,7),ncol=3))
layout(GtTE)
hist(TE1[which(TE1>0)],freq=F,breaks=20,xlab="Transmission Efficiency",main="1 <= log10[Gametocyte] < 1.5",ylim=c(0,8),xlim=c(0,1))
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE1fitZI$estimate[1],TE1fitZI$estimate[2]))
hist(TE2[which(TE2>0)],freq=F,breaks=20,xlab="Transmission Efficiency",main="1.5 <= log10[Gametocyte] < 2",ylim=c(0,8))
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE2fitZI$estimate[1],TE2fitZI$estimate[2]))
hist(TE3[which(TE3>0)],freq=F,breaks=20,xlab="Transmission Efficiency",main="2 <= log10[Gametocyte] < 2.5",ylim=c(0,8))
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE3fitZI$estimate[1],TE3fitZI$estimate[2]))
hist(TE4[which(TE4>0)],freq=F,breaks=20,xlab="Transmission Efficiency",main="2.5 <= log10[Gametocyte] < 3",ylim=c(0,8))
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE4fitZI$estimate[1],TE4fitZI$estimate[2]))
hist(TE5[which(TE5>0)],freq=F,breaks=20,xlab="Transmission Efficiency",main="3 <= log10[Gametocyte] < 3.5",ylim=c(0,8))
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE5fitZI$estimate[1],TE5fitZI$estimate[2]))
hist(TE6[which(TE6>0)],freq=F,breaks=20,xlab="Transmission Efficiency",main="3.5 <= log10[Gametocyte] < 4",ylim=c(0,8)) 
lines(seq(0,1,.01),dbeta(seq(0,1,.01),TE6fitZI$estimate[1],TE6fitZI$estimate[2]))


### Zero Inflation
alphas = c(TE1fitZI$estimate[1],TE2fitZI$estimate[1],TE3fitZI$estimate[1],TE4fitZI$estimate[1],TE5fitZI$estimate[1],TE6fitZI$estimate[1])
betas = c(TE1fitZI$estimate[2],TE2fitZI$estimate[2],TE3fitZI$estimate[2],TE4fitZI$estimate[2],TE5fitZI$estimate[2],TE6fitZI$estimate[2])
ZI = c(sum(TE1==0)/length(TE1),sum(TE2==0)/length(TE2),sum(TE3==0)/length(TE3),sum(TE4==0)/length(TE4),sum(TE5==0)/length(TE5),sum(TE6==0)/length(TE6))
xx = c(1.25,1.75,2.25,2.75,3.25,3.75)
sigfitZI = nls(ZI~p1*exp(p2*xx)/(p3+exp(p2*xx)),start=list(p1=.7,p2=-2,p3=.1))
plot(xx,ZI,xlab="Log10 Gametocytemia",ylab="Zero-Inflation in TE",ylim=c(0,.8))
tt = seq(1,4,.1)
lines(tt,sigmoidFev(tt,.8,-1.6,.029))

ZIsmooth = function(x,dx){
  output = rep(0,length(x))
  for(i in 1:length(x)){
    temp = Mosq[which(log10(G)>=x[i]-dx & log10(G)<x[i]+dx & !is.na(Mosq))]
    output[i] = sum(temp==0)/length(temp)
  }
  return(output)
}
par(mfrow=c(1,1))

plot(seq(1,4,.2),ZIsmooth(seq(1,4,.2),.2),xlab="Log10 Gametocytemia",ylab="Fraction with Zero Transmissions",main="Degree of Zero-Inflation in TE Given Gametocytemia",ylim=c(0,1))

xx = seq(1,4,.2)
ZIy = ZIsmooth(xx,.2)
sigfitZI = nls(ZIy~p1*exp(p2*xx)/(p3+exp(p2*xx))+p4,start=list(p1=.62,p2=-1,p3=.5,p4=.1))
sigmoidZI = function(x,p1,p2,p3,p4){
  p1*exp(p2*x)/(p3+exp(p2*x))+p4
}
lines(seq(1,4,.01),sigmoidZI(seq(1,4,.01),.5816,-2.9586,.001228,.0931))


## cumulative prop treated over time
Treated = matrix(NaN,nrow=2000,ncol=length(pID))
for(i in 1:length(pID)){
  Treated[1:length(which(MT$patient==pID[i])),i] = MT$Tretmentfree.[which(MT$patient==pID[i])]
}
cumTreat = abs(rowMeans(Treated,na.rm=T)-4)[1:300]
plot(cumTreat[1:30])
