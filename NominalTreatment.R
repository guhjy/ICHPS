rm(list = ls())
library(lmtest)
library(MASS)
library(nbpMatching)
library(Kendall)
library(mosaic)
library(multcomp)
library(nnet)
library(twang)
library(dplyr)
library(car)
library(RCurl)
library(Matching)
set.seed(0)


# Our interest is the causal effect of risky behavior (None, Alcohol only, Alcohol and Smoking) on blood pressure.

x <- getURL("https://docs.google.com/spreadsheets/d/1aC7nBSDxXX6UZo5pUGuWD2BNiuapRFwqOTsP4xWLl2w/pub?output=csv")
dat<-read.csv(text = x)
dat<-arrange(dat,Nom.Treatment)
dat<-dat %>%
  dplyr::select(-Race)  #Drop race variable 
Vars<-colnames(dat)[c(1:14)]
Vars[1]<-"as.factor(educO)"
Vars[2]<-"as.factor(incO)"
Vars[3]<-"as.factor(foodscrO)"
dat<-filter(dat,bpsystolic!="NA")
var.list<-paste(Vars,collapse="+")
var.list<-paste("Nom.Treatment",var.list,sep="~")
head(dat)

#Adding a variable with more bias (for sake of comparison)
bias<-1
dat$mealsfafh<-rt(nrow(dat),7,ncp=c(rep(0,1056),rep(bias,1345),rep(2*bias,522)))



#Pairwise comparisons

fit1<-multinom(as.formula(var.list),data=dat)
temp<-data.frame(fitted(fit1))
colnames(temp)<-c("p.r","p.s","p.t")
m.2<-cbind(dat,temp)
m.2$pr.rs<-m.2$p.r/(m.2$p.r+m.2$p.s)
m.2$pr.rt<-m.2$p.r/(m.2$p.r+m.2$p.t)

bwplot(p.r~Nom.Treatment,main="Probability of neither",data=m.2)
bwplot(p.s~Nom.Treatment,main="Probability of Alcohol only",data=m.2)
bwplot(p.t~Nom.Treatment,main="Probability of Alcohol & Smoking",data=m.2)

table(m.2$Nom.Treatment)
index.0<-m.2$Nom.Treatment=="0.None"
index.1<-m.2$Nom.Treatment=="Alcohol.only"
index.2<-m.2$Nom.Treatment=="Alcohol.Smoker"

#GPS distributions of full sample

par(mar=c(4,4,4,1))
par(mfrow=c(1,1))
df<-data.frame(x1=m.2$p.s,x2=m.2$p.t,
    Treatment=c(rep("T=r",sum(index.0)),rep("T=s",sum(index.1)),
    rep("T=t",sum(index.2))))

with(df, dataEllipse(logit(x1), logit(x2), Treatment, id.n=0, pch=15:17, 
                     main="95% ellipsoids",
     group.labels=c("T = r", "T = s", "T = t"), xlab="logit(P(T=s|X))",
     ylab="logit(P(T=t|X))",level=.95, fill=TRUE, fill.alpha=0.1))


plot(-30,1,xlim=c(-7,5),ylim=c(-8,4),xlab="logit(P(T=s|X))",
     ylab="logit(P(T=t|X))",main="95% ellipsoids")

with(df, dataEllipse(logit(x1), logit(x2), Treatment,plot.points=FALSE,
                     group.labels=c("Neither", "A. only", "A & S"),level=.95, fill=TRUE, fill.alpha=0.1))

#Pairwise matching

m.2rs<-filter(m.2,Nom.Treatment!="Alcohol.Smoker")
m.2rt<-filter(m.2,Nom.Treatment!="Alcohol.only")

yy1 = Match(Y=NULL, Tr= m.2rs$Nom.Treatment=="Alcohol.only", m.2rs$pr.rs,replace=FALSE,estimand="ATC",caliper=0.1)
yy2 = Match(Y=NULL, Tr= m.2rt$Nom.Treatment=="Alcohol.Smoker",m.2rt$pr.rt,replace=FALSE,estimand="ATC",caliper=0.1)

oo1 = yy1$index.control[which(yy1$index.control %in% yy2$index.control)]
oo2 = yy1$index.treated[which(yy1$index.control %in% yy2$index.control)]
oo3 = yy2$index.treated[which(yy2$index.control %in% yy1$index.control)]
tab.T<-table(m.2$Nom.Treatment)
n.r<-tab.T[1]
n.s<-tab.T[2]
n.t<-tab.T[3]

#GPS distributions among those matched using pairwise procedure

df<-data.frame(x1=c(m.2$p.s[oo1],m.2$p.s[oo2],m.2$p.s[oo3+n.s]),
               x2=c(m.2$p.t[oo1],m.2$p.t[oo2],m.2$p.t[oo3+n.s]),
     Treatment=c(rep("T=r",length(m.2$pr.rt[oo1])),rep("T=s",length(m.2$pr.rs[oo2])),
                 rep("T=t",length(m.2$pr.rs[oo3+n.s]))),ID = c(oo1,oo2,(oo3++n.s)) )

with(df, dataEllipse(logit(x1), logit(x2), Treatment, id.n=0, pch=15:17, center.pch="+",
    main="95% ellipsoids, binary matching",
    group.labels=c("T = r", "T = s", "T = t"), xlab="P(T=s|X)",ylab="P(T=t|X)",
    level=.95, fill=TRUE, fill.alpha=0.1,xlim=c(-5, 3),ylim=c(-6, 3)))

plot(-30,1,xlim=c(-2.5, 3),ylim=c(-6, 3),xlab="logit(P(T=s|X))",
     ylab="logit(P(T=t|X))",main="95% ellipsoids, binary matching")

with(df, dataEllipse(logit(x1), logit(x2), Treatment, id.n=0, plot.points=FALSE, center.pch="+",
    group.labels=c("T = r", "T = s", "T = t"), level=.95, fill=TRUE, fill.alpha=0.1))


# Matching methods with multiple treatments

## Drop subjects with extreme GPSs

min.max.Ps<-m.2 %>% group_by(Nom.Treatment) %>% summarize(minR=min(p.r),maxR=max(p.r),minS=min(p.s),maxS=max(p.s),minT=min(p.t),maxT=max(p.t))
min.max.Ps

m.2$Eligible<-m.2$p.r>=max(min.max.Ps$minR) & m.2$p.r<=min(min.max.Ps$maxR) &
   m.2$p.s>=max(min.max.Ps$minS) & m.2$p.s<=min(min.max.Ps$maxS) &
   m.2$p.t>=max(min.max.Ps$minT) & m.2$p.t<=min(min.max.Ps$maxT)
table(m.2$Eligible)

m.2<-filter(m.2,Eligible)

##Re-fit multinomial model

fit2<-multinom(as.formula(var.list),data=m.2)
temp<-data.frame(fitted(fit2))
colnames(temp)<-c("p.r","p.s","p.t")
m.3<-m.2 %>%
  dplyr::select(-p.r,-p.s,-p.t)
m.2<-cbind(m.3,temp)
m.2$pr.rs<-m.2$p.r/(m.2$p.r+m.2$p.s)
m.2$pr.rt<-m.2$p.r/(m.2$p.r+m.2$p.t)


### Multiple regression model (no weighting, no matching)

var.list<-paste(Vars, collapse="+")
var.list<-paste(var.list,"Nom.Treatment",sep="+")
var.list<-paste("bpsystolic",var.list,sep="~")
fit.MR<-lm(as.formula(var.list),data=dat)
summary(fit.MR)

### Inverse probability weighting, weights from multinomial logit model

m.2$prob.T<-ifelse(m.2$Nom.Treatment=="0.None",m.2$p.r,ifelse(m.2$Nom.Treatment=="Alcohol.only",
    m.2$p.s,m.2$p.t))
m.2$ipw<-1/m.2$prob.T

#How we'll judge bias: maximum pairwise absolute standardized bias

#First, some functions we'll need
Vars.Bias<-colnames(dat)[c(1:14)]
cols<-c(1:14)
funcMSB<-function(SD,x1,x2,x3){max(abs((mean(x1)-mean(x2))/SD),abs((mean(x1)-mean(x3))/SD),abs((mean(x3)-mean(x2))/SD))}
funcMSB2<-function(SD,x1,x2,x3){max(abs(x1-x2)/SD,abs(x1-x3)/SD,abs(x3-x2)/SD)}
sds<-NULL
for (i in 1:length(Vars.Bias)){sds[i]<-sd(m.2[m.2$Nom.Treatment=="0.None",cols[i]])}


### Covariates bias before weighting, calculated using maximum pairwise SB
Start.msbs<-NA
for (i in 1:length(Vars.Bias)){
  Start.msbs[i]<-funcMSB2(sds[i],mean(m.2[m.2$Nom.Treatment=="0.None",cols[i]]),mean(m.2[m.2$Nom.Treatment=="Alcohol.only",cols[i]]),
                            mean(m.2[m.2$Nom.Treatment=="Alcohol.Smoker",cols[i]]))}
names(Start.msbs)<-colnames(m.2)[1:14]
round(Start.msbs,3)

#### Example of bias:
bwplot(age~Nom.Treatment,data=m.2)
tally(incO~Nom.Treatment,data=m.2,format="proportion")

### After weighting using weights from multinomial logistic regression
IPTW.msbs<-NULL
for (i in 1:length(Vars.Bias)){
  IPTW.msbs[i]<-funcMSB(sds[i],weighted.mean(m.2[m.2$Nom.Treatment=="0.None",cols[i]],m.2[m.2$Nom.Treatment=="0.None",]$ipw),
      weighted.mean(m.2[m.2$Nom.Treatment=="Alcohol.only",cols[i]],m.2[m.2$Nom.Treatment=="Alcohol.only",]$ipw),
      weighted.mean(m.2[m.2$Nom.Treatment=="Alcohol.Smoker",cols[i]],m.2[m.2$Nom.Treatment=="Alcohol.Smoker",]$ipw))}
names(IPTW.msbs)<-colnames(m.2)[1:14]
round(IPTW.msbs,3)


### Weights using Generalized Boosted Models

var.list<-paste(Vars.Bias,collapse="+")
var.list<-paste("Nom.Treatment",var.list,sep="~")
m.3<-m.2[,colnames(m.2)%in%c("Nom.Treatment",Vars.Bias,"bpsystolic")]


ps.out1<- mnps(as.formula(var.list), data=m.2,stop.method = c("es.mean", "ks.mean"),treatATT="0.None",verbose=FALSE,
             n.trees=10000, interaction.depth=2, estimand="ATT")
m.2$ipw.GBM = get.weights(ps.out1)

#### Code to use GBM for ATE
#ps.out1<- mnps(as.formula(var.list), data=m.2,stop.method = c("es.mean", "ks.mean"),verbose=FALSE,
#               n.trees=10000, interaction.depth=2, estimand="ATE")


#### A command using GBM to look at the covariates' bias
#bal.table(ps.out)

#Covariates bias after weighting using GBM
GBM.msbs<-NULL
for (i in 1:length(cols)){
  GBM.msbs[i]<-funcMSB(sds[i],weighted.mean(m.2[m.2$Nom.Treatment=="0.None",cols[i]],m.2[m.2$Nom.Treatment=="0.None",]$ipw.GBM),
                        weighted.mean(m.2[m.2$Nom.Treatment=="Alcohol.only",cols[i]],m.2[m.2$Nom.Treatment=="Alcohol.only",]$ipw.GBM),
                        weighted.mean(m.2[m.2$Nom.Treatment=="Alcohol.Smoker",cols[i]],m.2[m.2$Nom.Treatment=="Alcohol.Smoker",]$ipw.GBM))}
names(GBM.msbs)<-colnames(m.2)[1:14]
round(GBM.msbs,3)


# Vector Matching
quintiles.t<-quantile(m.2$pr.rt,c(.2,.4,.6,.8))
clustnum<-5

m.2$Quint.t<-1
temp.t<-kmeans(logit(m.2$p.t),clustnum)
m.2$Quint.t<-temp.t$cluster

m.2$Quint.s<-1
temp.s<-kmeans(logit(m.2$p.s),clustnum)
m.2$Quint.s<-temp.s$cluster

temp.s<-filter(m.2,Nom.Treatment!="Alcohol.Smoker")
temp.t<-filter(m.2,Nom.Treatment!="Alcohol.only")

#see Abadie_Imbens paper (Table 1)
Ms<-Matchby(Y=temp.s$bpsystolic, Tr=temp.s$Nom.Treatment=="0.None",X=logit(temp.s$p.r),by=temp.s$Quint.t,calip=.1,replace=T,estimand="ATT") 
Mt<-Matchby(Y=temp.t$bpsystolic, Tr=temp.t$Nom.Treatment=="0.None",X=logit(temp.t$p.r),by=temp.t$Quint.s,calip=.1,replace=T,estimand="ATT")


rownames(m.2)<-1:nrow(m.2)
m.2$id<-1:nrow(m.2)
m.2$r.both<-m.2$id%in%Ms$index.treated&m.2$id%in%Mt$index.treated
m.temp<-m.2[m.2$r.both=="TRUE",]
m.2$match.r<-NULL
rsmatch<-cbind(Ms$index.treated,Ms$index.control)
rtmatch<-cbind(Mt$index.treated,(sum(m.2$Nom.Treatment=="Alcohol.only")+Mt$index.control))
rsmatch<-rsmatch[rsmatch[,1]%in%rownames(m.temp),]
rtmatch<-rtmatch[rtmatch[,1]%in%rownames(m.temp),]
triplets<-cbind(rsmatch[order(rsmatch[,1]),],rtmatch[order(rtmatch[,1]),])
triplets<-as.matrix(triplets[,c(1,2,4)])
n.trip<-nrow(triplets)
n.tripP<-nrow(triplets)/sum(m.2$Nom.Treatment=="0.None")


VM.msbs<-NA
for (i in 1:length(cols)){
  VM.msbs[i]<-funcMSB(sds[i],m.2[,i][triplets[,1]],m.2[,i][triplets[,2]],m.2[,i][triplets[,3]])}
names(VM.msbs)<-colnames(m.2)[1:14]
round(VM.msbs,3)


##Comparison of ellipsoids using matched cohort

df<-data.frame(
  x1=c(m.2$p.s[triplets[,1]],m.2$p.s[triplets[,2]],m.2$p.s[triplets[,3]]),
  x2=c(m.2$p.t[triplets[,1]],m.2$p.t[triplets[,2]],m.2$p.t[triplets[,3]]),
  Treatment=c(rep("T=r",n.trip),rep("T=s",n.trip),rep("T=t",n.trip)) )

plot(-30,1,xlim=c(-3,2),ylim=c(-7,2),xlab="logit(P(T=s|X))",
     ylab="logit(P(T=t|X))",main="95% ellipsoids, matched on a vector")

with(df, dataEllipse(logit(x1), logit(x2), Treatment, id.n=0, plot.points=FALSE, center.pch="+",
                     group.labels=c("T = r", "T = s", "T = t"), level=.95, fill=TRUE, fill.alpha=0.1))

## This time, including full reference group distribution
df<-data.frame(
  x1=c(m.2$p.s[m.2$Nom.Treatment=="0.None"],m.2$p.s[triplets[,1]],m.2$p.s[triplets[,2]],m.2$p.s[triplets[,3]]),
  x2=c(m.2$p.t[m.2$Nom.Treatment=="0.None"],m.2$p.t[triplets[,1]],m.2$p.t[triplets[,2]],m.2$p.t[triplets[,3]]),
  Treatment=c(rep("Null",738),rep("T=r",n.trip),rep("T=s",n.trip),rep("T=t",n.trip)) )

plot(-30,1,xlim=c(-3,2),ylim=c(-7,2),xlab="logit(P(T=s|X))",
     ylab="logit(P(T=t|X))",main="95% ellipsoids, matching on a vector")

with(df, dataEllipse(logit(x1), logit(x2), Treatment, id.n=0, plot.points=FALSE, center.pch="+",
                     group.labels=c("T = r", "T = s", "T = t"), level=.95, fill=TRUE, fill.alpha=0.1))


## Overall comparison of covariates' bias

covariates.bias<-round(rbind(Start.msbs,GBM.msbs,VM.msbs,IPTW.msbs),3)
covariates.bias


#Analysis phase


#Estimating treatment effects using MLR weights

EY.r = sum(m.2$bpsystolic*as.numeric(m.2$Nom.Treatment=="0.None")/m.2$p.r)*
  sum(as.numeric(m.2$Nom.Treatment=="0.None")/m.2$p.r)^(-1)
EY.s = sum(m.2$bpsystolic*as.numeric(m.2$Nom.Treatment=="Alcohol.only")/m.2$p.s)*
  sum(as.numeric(m.2$Nom.Treatment=="Alcohol.only")/m.2$p.s)^(-1)
EY.t = sum(m.2$bpsystolic*as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")/m.2$p.t)*
  sum(as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")/m.2$p.t)^(-1)

VarY<-var(m.2$bpsystolic)
ATEW<-cbind(EY.r,EY.s,EY.t)

VY.r = sum(VarY*((as.numeric(m.2$Nom.Treatment=="0.None")/m.2$p.r)/
                   sum(as.numeric(m.2$Nom.Treatment=="0.None")/m.2$p.r))^2)
VY.s =  sum(VarY*((as.numeric(m.2$Nom.Treatment=="Alcohol.only")/m.2$p.s)/
                    sum(as.numeric(m.2$Nom.Treatment=="Alcohol.only")/m.2$p.s))^2)
VY.t =  sum(VarY*((as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")/m.2$p.t)/
                    sum(as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")/m.2$p.t))^2)
Var.W<-cbind(VY.r,VY.s,VY.t)

data.frame(Effect=ATEW[2]-ATEW[1],Var=sqrt(Var.W[1]+Var.W[2])) #Alcohol only versus none
data.frame(Effect=ATEW[3]-ATEW[1],Var=sqrt(Var.W[1]+Var.W[3])) #Alcohol and smoking versus none
data.frame(Effect=ATEW[3]-ATEW[2],Var=sqrt(Var.W[2]+Var.W[3])) #Alcohol and smoking versus alcohol only

IPTW.ATE<-round(c(ATEW[2]-ATEW[1],ATEW[3]-ATEW[1],ATEW[3]-ATEW[2]),3)

#Estimating treatment effects using GBM weights

EY.r = sum(m.2$bpsystolic*as.numeric(m.2$Nom.Treatment=="0.None")*m.2$ipw.GBM)*
  sum(as.numeric(m.2$Nom.Treatment=="0.None")*m.2$ipw.GBM)^(-1)
EY.s = sum(m.2$bpsystolic*as.numeric(m.2$Nom.Treatment=="Alcohol.only")*m.2$ipw.GBM)*
  sum(as.numeric(m.2$Nom.Treatment=="Alcohol.only")*m.2$ipw.GBM)^(-1)
EY.t = sum(m.2$bpsystolic*as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")*m.2$ipw.GBM)*
  sum(as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")*m.2$ipw.GBM)^(-1)

VarY<-var(m.2$bpsystolic)
ATEW<-cbind(EY.r,EY.s,EY.t)

VY.r = sum(VarY*((as.numeric(m.2$Nom.Treatment=="0.None")*m.2$ipw.GBM)/
                   sum(as.numeric(m.2$Nom.Treatment=="0.None")*m.2$ipw.GBM))^2)
VY.s =  sum(VarY*((as.numeric(m.2$Nom.Treatment=="Alcohol.only")*m.2$ipw.GBM)/
                    sum(as.numeric(m.2$Nom.Treatment=="Alcohol.only")*m.2$ipw.GBM))^2)
VY.t =  sum(VarY*((as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")*m.2$ipw.GBM)/
                    sum(as.numeric(m.2$Nom.Treatment=="Alcohol.Smoker")*m.2$ipw.GBM))^2)
Var.W<-cbind(VY.r,VY.s,VY.t)

data.frame(Effect=ATEW[2]-ATEW[1],Var=sqrt(Var.W[1]+Var.W[2])) #Alcohol only versus none
data.frame(Effect=ATEW[3]-ATEW[1],Var=sqrt(Var.W[1]+Var.W[3])) #Alcohol and smoking versus none
data.frame(Effect=ATEW[3]-ATEW[2],Var=sqrt(Var.W[2]+Var.W[3])) #Alcohol and smoking versus alcohol only

GBM.ATT<-round(c(ATEW[2]-ATEW[1],ATEW[3]-ATEW[1],ATEW[3]-ATEW[2]),3)


# Treatment effects using VM

VM.ATT<-round(c(Ms$est,Mt$est,Mt$est-Ms$est),3)

TEs<-rbind(IPTW.ATE,GBM.ATT,VM.ATT)
colnames(TEs)<-c("Treatment 2 versus 1","Treatment 3 versus 1","Treatment 3 versus 2")
TEs





