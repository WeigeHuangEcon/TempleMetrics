rm(list = ls())
# data=read.csv('/Users/weigehuang/Documents/Mac/Dissertation/JobeMarketPaper/factors.csv')
setwd('/Users/weigehuang/Documents/Mac/Dissertation/JobeMarketPaper')
getwd()
data=read.csv('/Users/weigehuang/Documents/Mac/Dissertation/JobeMarketPaper/F-F_Research_Data_5_Factors_2x3.CSV')
attach(data)
library(cfa)
x=cbind(Mkt.RF,SMB,RMW,CMA)
t=HML
y=Lo.30-RF
# y=Med.40-RF
n_ctvals=20
ctvals=seq(quantile(t,.05), quantile(t, .95), length.out=n_ctvals)
s=20
e=0.05

# to obtain conditional quantiles estimators, we estimate S linear quantile regressions of Y on X
fit=getCondQuants.qr(Y=y,X=x,t=t,s=s,e=e)

# plot coefficients on log family income for each model corresponding to different quantiles "tau"
coefficients=fit$coefficients
coefficients=fit$coefficients[2,]
plot(coefficients,type = 'l',ylab = 'Conefficients', xlab = 's-Quantiles')


# to obtain the conditional counterfactual quantiles
# t is the treatment variable, ct is the counterfactual treatments. Both of them can be a set of varibles, such as t=c('t','age')
b=getCountQuant.qr(Y=y,X=x,t=t,ct=ctvals[2],s=s,e=e )

# to obtain the counterfactual distribution functions
fs=makeCountDisFs.qr(Y=y,X=x,t=t,ct=ctvals[2],s=s,e=e)
fs1=makeCountDisFs.qr(Y=y,X=x,t=t,ct=ctvals[4],s=s,e=e)

ysort=sort(y)
plot(ysort,fs,type = 'l')
lines(ysort,fs1,col='red')

quants=c(0.1,0.25,0.5,0.75,0.9)
# convert counterfactual distribution functions to quantiles
DisFs2Quants.qr(Y=y,X=x,t=t,ct=ctvals[5],quants,s=s,e=e)

# obtain unconditional counterfactual quantiles
cQMat=getCountQuants.qr (Y=y,X=x,t=t,cts =ctvals ,quants=quants,s=s,e=e)
matplot(ctvals,t(cQMat),type="l")

plot(ctvals,cQMat[1,],type = 'l')
plot(ctvals,cQMat[2,],type = 'l')
plot(ctvals,cQMat[3,],type = 'l')
plot(ctvals,cQMat[4,],type = 'l')
plot(ctvals,cQMat[5,],type = 'l')

# bootstrap
boot.getCountQuants.qr(B=1,Y=y,X=x,t=t,cts=ctvals,quants=quants,s=s,e=e)


library(parallel)
pt <- proc.time()

sdMat=mclapply(1:2,boot.getCountQuants,Y=y,X=x,t=t,cts=ctvals,quants=quants,s=s,e=e,mc.cores = detectCores()-1)

proc.time() - pt

# obtain standard deviation
sdmat=getSDMat.qr(B=2,Y=y,X=x,t=t,cts=ctvals,quants=quants,s=s,e=e)

# to plot quantiles against treatments with confidence bands, have to obtain quantiles, standerd deviation
CountQuants.plot(CountQuants=cQMat,SDMat=sdmat)

par(mfrow=c(3,2))
plot(ctvals,cQMat[1,],type = 'l',ylab = 'mean',xlab = 't',ylim = c(min(cQMat[1,]-1.96*sdmat[1,]),max(cQMat[1,]+1.96*sdmat[1,])))
lines(ctvals,cQMat[1,]+1.96*sdmat[1,],lty=2)
lines(ctvals,cQMat[1,]-1.96*sdmat[1,],lty=2)

plot(ctvals,cQMat[2,],type = 'l',ylab = '10%',xlab = 't',ylim = c(min(cQMat[2,]-1.96*sdmat[2,]),max(cQMat[2,]+1.96*sdmat[2,])))
lines(ctvals,cQMat[2,]+1.96*sdmat[2,],lty=2)
lines(ctvals,cQMat[2,]-1.96*sdmat[2,],lty=2)

plot(ctvals,cQMat[3,],type = 'l',ylab = '25%',xlab = 't',ylim = c(min(cQMat[3,]-1.96*sdmat[3,]),max(cQMat[3,]+1.96*sdmat[3,])))
lines(ctvals,cQMat[3,]+1.96*sdmat[3,],lty=2)
lines(ctvals,cQMat[3,]-1.96*sdmat[3,],lty=2)

plot(ctvals,cQMat[4,],type = 'l',ylab = '50%',xlab = 't',ylim = c(min(cQMat[4,]-1.96*sdmat[4,]),max(cQMat[4,]+1.96*sdmat[4,])))
lines(ctvals,cQMat[4,]+1.96*sdmat[4,],lty=2)
lines(ctvals,cQMat[4,]-1.96*sdmat[4,],lty=2)

plot(ctvals,cQMat[5,],type = 'l',ylab = '75%',xlab = 't',ylim = c(min(cQMat[5,]-1.96*sdmat[5,]),max(cQMat[5,]+1.96*sdmat[5,])))
lines(ctvals,cQMat[5,]+1.96*sdmat[5,],lty=2)
lines(ctvals,cQMat[5,]-1.96*sdmat[5,],lty=2)

plot(ctvals,cQMat[6,],type = 'l',ylab = '90%',xlab = 't',ylim = c(min(cQMat[6,]-1.96*sdmat[6,]),max(cQMat[6,]+1.96*sdmat[6,])))
lines(ctvals,cQMat[6,]+1.96*sdmat[6,],lty=2)
lines(ctvals,cQMat[6,]-1.96*sdmat[6,],lty=2)

dev.off()



CountQuants=cQMat
SDMat=sdmat
library(ggplot2)
library(reshape2)
cqs=melt(CountQuants)
sd=melt(SDMat)
sqs=merge(cqs,sd,by=c("Var1","Var2"))
colnames(sqs) <- c("quants","ctvals","cqs","sd")

ggcp=ggplot2::ggplot(data=sqs, aes(ctvals, cqs, ymax=cqs+1.96*sd,
                                               ymin=cqs-1.96*sd),group=quants) +
  ggplot2::geom_line(aes(ctvals, cqs, group=quants, color=quants)) +
  facet_wrap(~quants) +
  # ggplot2::geom_errorbar(size=.3, width=.02) +
  ggplot2::geom_line(aes(ctvals, cqs+1.96*sd,group=quants, color=quants), linetype="dashed") +
  ggplot2::geom_line(aes(ctvals, cqs-1.96*sd,group=quants, color=quants), linetype="dashed") +
  ggplot2::geom_point(aes(ctvals, cqs, group=quants,color=quants))  +
  ggplot2::scale_y_continuous("Unconditional Counterfactual Quantiles") + ##, limits=c(-.4, .4)) +
  ggplot2::scale_x_continuous("Treatments") +   ##limits=c(0,1), breaks=c(.1,.3,.5,.7,.9)) +
  ggplot2::theme_classic() +
  ggplot2::theme(panel.border = element_rect(colour = 'black', size=1,
                                             fill=NA,
                                             linetype='solid'))
ggcp



