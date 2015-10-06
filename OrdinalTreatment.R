rm(list = ls())
library(lmtest)
library(MASS)
library(nbpMatching)
library(Kendall)
library(multcomp)
library(dplyr)
library(RCurl)

# Our interest is the causal effect of label use on BMI

x <- getURL("https://docs.google.com/spreadsheets/d/1aC7nBSDxXX6UZo5pUGuWD2BNiuapRFwqOTsP4xWLl2w/pub?output=csv")
dat<-read.csv(text = x)

## Generalized propensity score model (proportional odds)

polr.out <- polr(as.ordered(labeluse)~as.factor(educO) + as.factor(incO)+as.factor(foodscrO)+
                   prior.bmi+age+hsize+mealsfafh+female+healthydiet+Race+depression+docdiabet+
                   doctadv1+doctadv2+knowfgp,data=dat)
summary(polr.out)


# It is generally recommended to drop subjects with extreme linear predictors and to refit model

dat$linpred<-polr.out$lp
minLPs<-dat %>% group_by(labeluse) %>% summarize(minLP=min(linpred))
max.lower<-max(minLPs[,2])
maxLPs<-dat %>% group_by(labeluse) %>% summarize(maxLP=max(linpred))
min.upper<-min(maxLPs[,2])
boxplot(dat$linpred~dat$labeluse,xlab="Label Use (Low to High)")
abline(max.lower,0,lty=2,col="red")
abline(min.upper,0,lty=2,col="red")
dat$Eligible<-dat$linpred<=min.upper&dat$linpred>=max.lower
table(dat$Eligible)

# To drop those ineligible
dat1<-filter(dat,linpred<=min.upper,linpred>=max.lower)



## Re-fit the model
polr.out <- polr(as.ordered(labeluse)~as.factor(educO) + as.factor(incO)+as.factor(foodscrO)+
                   prior.bmi+age+hsize+mealsfafh+female+healthydiet+Race+depression+docdiabet+
                   doctadv1+doctadv2+knowfgp,data=dat1)
dat1$linpred<-polr.out$lp



# Subclassification on the linear predictor

dat1$Subclass<-cut(dat1$linpred, breaks=c(quantile(dat1$linpred, probs = seq(0, 1, by = 0.20))), 
    labels=c("K1","K2","K3","K4","K5"))


### Look at the bias: Kendall's Tau

#### This is done for the original data, as well as an average of the within subclass bias

Vars<-colnames(dat)[c(1:9,11:15)]
Bias.Kendall<-data.frame(Var=Vars,Original=NA,KFive=NA)
for (i in 1:nrow(Bias.Kendall)){
  Bias.Kendall[i,]$Original<-Kendall(dat1$labeluse,dat1[,Vars[i]])$tau[1]
  temp.subclass<-NULL
  for (j in 1:5){
    temp<-filter(dat1,Subclass==unique(dat1$Subclass)[j])
    temp.subclass[j]<-Kendall(temp$labeluse,temp[,Vars[i]])$tau[1]
  }
  Bias.Kendall[i,]$KFive<-mean(temp.subclass)
}

Bias.Kendall<-arrange(Bias.Kendall, Original)


#Plotting the Kendall's Tau, both pre and post subclassification

par(mar=c(3,3,3,3))
par(mfrow=c(1,1))
plot(x=0,y=0,ylim=c(0,14),col="white",xlim=c(-0.40,.3),yaxt='n',xlab="Kendall's Tau",xaxt='n')
axis(1,c(-0.2,-0.1,0,0.1,0.2,0.3))
temp<-seq(1:nrow(Bias.Kendall))
for (i in 1:nrow(Bias.Kendall)){
  lines(x=c(-.2,2),y=c(temp[i],temp[i]),lty="dotted",col="grey")}
for (i in 1:nrow(Bias.Kendall)){
  points(x=c(Bias.Kendall[i,2],Bias.Kendall[i,3]), y=c(temp[i],temp[i]),pch=c(1,8),col="black",xaxt='n')
  lines(x=c(Bias.Kendall[i,2],Bias.Kendall[i,3]), y=c(temp[i],temp[i]), lty="solid",xaxt='n', col="black")
  j<-15-i
  text(-.40,temp[i],Bias.Kendall[i,1],cex=.71,adj=0)
}
legend(.14,3.8,c(expression(paste(tau, ", pre ")),
                   expression(paste(bar(tau), ", post (K = 5) "))),pch=c(1,8),cex=.8)
lines(x=c(0,0),y=c(0,33))


# Analysis 


## Anova
anova(lm(bmi~as.factor(Subclass)+as.factor(labeluse) + as.factor(educO)+as.factor(incO)+as.factor(foodscrO)+prior.bmi+
           age+hsize+mealsfafh+female+healthydiet+depression+docdiabet+doctadv1+doctadv2+knowfgp,data=dat1))

## Subclassification

### Categorize ordinal variables

Vars[1]<-"as.factor(educO)"
Vars[2]<-"as.factor(incO)"
Vars[3]<-"as.factor(foodscrO)"
var.list<-paste(Vars,collapse="+")
var.list2<-paste("as.factor(labeluse)","+",var.list,"-1")
var.list3<-paste("bmi",var.list2,sep="~")

Q.1<-summary(lm(as.formula(var.list3),data=filter(dat1,Subclass=="K1")),corr=TRUE)
Q.2<-summary(lm(as.formula(var.list3),data=filter(dat1,Subclass=="K2")),corr=TRUE)
Q.3<-summary(lm(as.formula(var.list3),data=filter(dat1,Subclass=="K3")),corr=TRUE)
Q.4<-summary(lm(as.formula(var.list3),data=filter(dat1,Subclass=="K4")),corr=TRUE)
Q.5<-summary(lm(as.formula(var.list3),data=filter(dat1,Subclass=="K5")),corr=TRUE)


coeff.Q<-matrix(nrow=5,ncol=5)

coeff.Q[1,]<-Q.1$coeff[1:5]
coeff.Q[2,]<-Q.2$coeff[1:5]
coeff.Q[3,]<-Q.3$coeff[1:5]
coeff.Q[4,]<-Q.4$coeff[1:5]
coeff.Q[5,]<-Q.5$coeff[1:5]

colnames(coeff.Q)<-c("T1","T2","T3","T4","T5")
coeffQ.2<-data.frame(coeff.Q)

#Functions for use:
fun<-function(DF, upp,low){
  c(DF$coeff[upp,1]-DF$coeff[low,1],
    sqrt(DF$coeff[upp,2]^2+DF$coeff[low,2]^2-2*DF$corr[low,upp]*DF$coeff[upp,2]*DF$coeff[low,2]))
}


coeffQ.2$two.one<-c(fun(Q.1,2,1)[1],fun(Q.2,2,1)[1],fun(Q.3,2,1)[1],fun(Q.4,2,1)[1],fun(Q.5,2,1)[1])
coeffQ.2$two.one.var<-c(fun(Q.1,2,1)[2],fun(Q.2,2,1)[2],fun(Q.3,2,1)[2],fun(Q.4,2,1)[2],fun(Q.5,2,1)[2])
coeffQ.2$three.one<-c(fun(Q.1,3,1)[1],fun(Q.2,3,1)[1],fun(Q.3,3,1)[1],fun(Q.4,3,1)[1],fun(Q.5,3,1)[1])
coeffQ.2$three.one.var<-c(fun(Q.1,3,1)[2],fun(Q.2,3,1)[2],fun(Q.3,3,1)[2],fun(Q.4,3,1)[2],fun(Q.5,3,1)[2])
coeffQ.2$four.one<-c(fun(Q.1,4,1)[1],fun(Q.2,4,1)[1],fun(Q.3,4,1)[1],fun(Q.4,4,1)[1],fun(Q.5,4,1)[1])
coeffQ.2$four.one.var<-c(fun(Q.1,4,1)[2],fun(Q.2,4,1)[2],fun(Q.3,4,1)[2],fun(Q.4,4,1)[2],fun(Q.5,4,1)[2])
coeffQ.2$five.one<-c(fun(Q.1,5,1)[1],fun(Q.2,5,1)[1],fun(Q.3,5,1)[1],fun(Q.4,5,1)[1],fun(Q.5,5,1)[1])
coeffQ.2$five.one.var<-c(fun(Q.1,5,1)[2],fun(Q.2,5,1)[2],fun(Q.3,5,1)[2],fun(Q.4,5,1)[2],fun(Q.5,5,1)[2])
coeffQ.2$three.two<-c(fun(Q.1,3,2)[1],fun(Q.2,3,2)[1],fun(Q.3,3,2)[1],fun(Q.4,3,2)[1],fun(Q.5,3,2)[1])
coeffQ.2$three.two.var<-c(fun(Q.1,3,2)[2],fun(Q.2,3,2)[2],fun(Q.3,3,2)[2],fun(Q.4,3,2)[2],fun(Q.5,3,2)[2])
coeffQ.2$four.two<-c(fun(Q.1,4,2)[1],fun(Q.2,4,2)[1],fun(Q.3,4,2)[1],fun(Q.4,4,2)[1],fun(Q.5,4,2)[1])
coeffQ.2$four.two.var<-c(fun(Q.1,4,2)[2],fun(Q.2,4,2)[2],fun(Q.3,4,2)[2],fun(Q.4,4,2)[2],fun(Q.5,4,2)[2])
coeffQ.2$five.two<-c(fun(Q.1,5,2)[1],fun(Q.2,5,2)[1],fun(Q.3,5,2)[1],fun(Q.4,5,2)[1],fun(Q.5,5,2)[1])
coeffQ.2$five.two.var<-c(fun(Q.1,5,2)[2],fun(Q.2,5,2)[2],fun(Q.3,5,2)[2],fun(Q.4,5,2)[2],fun(Q.5,5,2)[2])
coeffQ.2$four.three<-c(fun(Q.1,4,3)[1],fun(Q.2,4,3)[1],fun(Q.3,4,3)[1],fun(Q.4,4,3)[1],fun(Q.5,4,3)[1])
coeffQ.2$four.three.var<-c(fun(Q.1,4,3)[2],fun(Q.2,4,3)[2],fun(Q.3,4,3)[2],fun(Q.4,4,3)[2],fun(Q.5,4,3)[2])
coeffQ.2$five.three<-c(fun(Q.1,5,3)[1],fun(Q.2,5,3)[1],fun(Q.3,5,3)[1],fun(Q.4,5,3)[1],fun(Q.5,5,3)[1])
coeffQ.2$five.three.var<-c(fun(Q.1,5,3)[2],fun(Q.2,5,3)[2],fun(Q.3,5,3)[2],fun(Q.4,5,3)[2],fun(Q.5,5,3)[2])
coeffQ.2$five.four<-c(fun(Q.1,5,4)[1],fun(Q.2,5,4)[1],fun(Q.3,5,4)[1],fun(Q.4,5,4)[1],fun(Q.5,5,4)[1])
coeffQ.2$five.four.var<-c(fun(Q.1,5,4)[2],fun(Q.2,5,4)[2],fun(Q.3,5,4)[2],fun(Q.4,5,4)[2],fun(Q.5,5,4)[2])



ComboQ<-table(dat1$Subclass,dat1$labeluse)
coeffQ.2$weights<-(rowSums(ComboQ))/sum(ComboQ)
EffectQ.21<-weighted.mean(coeffQ.2$two.one,coeffQ.2$weights)
VarQ.21<-sqrt(sum(coeffQ.2$two.one.var^2*coeffQ.2$weights^2))
EffectQ.31<-weighted.mean(coeffQ.2$three.one,coeffQ.2$weights)
VarQ.31<-sqrt(sum(coeffQ.2$three.one.var^2*coeffQ.2$weights^2))
EffectQ.41<-weighted.mean(coeffQ.2$four.one,coeffQ.2$weights)
VarQ.41<-sqrt(sum(coeffQ.2$four.one.var^2*coeffQ.2$weights^2))
EffectQ.51<-weighted.mean(coeffQ.2$five.one,coeffQ.2$weights)
VarQ.51<-sqrt(sum(coeffQ.2$five.one.var^2*coeffQ.2$weights^2))
EffectQ.32<-weighted.mean(coeffQ.2$three.two,coeffQ.2$weights)
VarQ.32<-sqrt(sum(coeffQ.2$three.two.var^2*coeffQ.2$weights^2))
EffectQ.42<-weighted.mean(coeffQ.2$four.two,coeffQ.2$weights)
VarQ.42<-sqrt(sum(coeffQ.2$four.two.var^2*coeffQ.2$weights^2))
EffectQ.52<-weighted.mean(coeffQ.2$five.two,coeffQ.2$weights)
VarQ.52<-sqrt(sum(coeffQ.2$five.two.var^2*coeffQ.2$weights^2))
EffectQ.43<-weighted.mean(coeffQ.2$four.three,coeffQ.2$weights)
VarQ.43<-sqrt(sum(coeffQ.2$four.three.var^2*coeffQ.2$weights^2))
EffectQ.53<-weighted.mean(coeffQ.2$five.three,coeffQ.2$weights)
VarQ.53<-sqrt(sum(coeffQ.2$five.three.var^2*coeffQ.2$weights^2))
EffectQ.54<-weighted.mean(coeffQ.2$five.four,coeffQ.2$weights)
VarQ.54<-sqrt(sum(coeffQ.2$five.four.var^2*coeffQ.2$weights^2))
effects.Q<-c(EffectQ.21, EffectQ.31, EffectQ.41, EffectQ.51,
             EffectQ.32, EffectQ.42, EffectQ.52, EffectQ.43, EffectQ.53, EffectQ.54)

VarQ<-c(VarQ.21, VarQ.31, VarQ.41, VarQ.51, VarQ.32, VarQ.42, VarQ.52, VarQ.43, VarQ.53, VarQ.54)


output<-t(data.frame(effects.Q,VarQ))
colnames(output)<-c("2 v 1","3 v 1","4 v 1","5 v 1","3 v 2","4 v 2","5 v 2","4 v 3","5 v 3","5 v 4")
round(output,2)




#Inverse probability weighting (will be shown with nominal treatments more extensively)
pred0 <- predict(polr.out, type = "p")
colnames(pred0)<-c("Pred1","Pred2","Pred3","Pred4","Pred5")
dat1<-cbind(dat1,pred0)
dat1$weight<-1/dat1$Pred1
dat1[dat1$labeluse==2,]$weight<-1/dat1[dat1$labeluse==2,]$Pred2
dat1[dat1$labeluse==3,]$weight<-1/dat1[dat1$labeluse==3,]$Pred3
dat1[dat1$labeluse==4,]$weight<-1/dat1[dat1$labeluse==4,]$Pred4
dat1[dat1$labeluse==5,]$weight<-1/dat1[dat1$labeluse==5,]$Pred5

sum(dat1$weight>10)
max(dat1$weight)


