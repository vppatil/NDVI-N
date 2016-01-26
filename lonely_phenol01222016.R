#NDVI percent N plots and non-linear curve fits
setwd('~/NDVI_N/code/')
#reading in the data
source('julian.r')
source('phenolcurve1172015.r')
source('aic_code.r')

load('../output_plots/ndvi.n.phenol01222016.R')
dir()

#csv file cleaned up by Dan.
#need to compare to original data file- confirm that removing na's in kyle's analysis resulted in consistent dataset being used for all models compared w/ AIC.
vegplots<-read.csv('../Data/VegPlots_for R.csv')

lonely.ndvi<-read.csv('../Data/Lonelybiomassndvi.csv')
lonely.cn<-read.csv('../Data/lonely_cn.csv')

library(AICcmodavg)
library(MASS) #fitdistr function in MASS allows one to fit good parameters.
library(minpack.lm)
library(gam)
library(ggplot2)

#could merge these two dataframes into one based on date and plot. Shouldn't functionally matter.
#make doy column
lonely.ndvi$Date=as.character(lonely.ndvi$Date)
ndvi.date.split=sapply(lonely.ndvi$Date,function(x) as.numeric(strsplit(x,split='/')[[1]]))
lonely.ndvi$Month=ndvi.date.split[1,]
lonely.ndvi$Day=ndvi.date.split[2,]
lonely.ndvi$Year=ndvi.date.split[3,]
lonely.ndvi$jdate=julian.df(lonely.ndvi)

lonely.cn$Date.collected=as.character(lonely.cn$Date.collected)
cn.Date.collected.split=sapply(lonely.cn$Date.collected,function(x) as.numeric(strsplit(x,split='/')[[1]]))
lonely.cn$Month=cn.Date.collected.split[1,]
lonely.cn$Day=cn.Date.collected.split[2,]
lonely.cn$Year=cn.Date.collected.split[3,]
lonely.cn$jdate=julian.df(lonely.cn)
#######################################################################
#scatterplots
plot(NDVI_WV.2~jdate,data=lonely.ndvi)
plot(X.N~jdate,data=lonely.cn)

#vars of interest are NDVI_WV.2,X.N, and jdate.
####################################################
#########################################################
#creating merged ndvi-n dataset for lonely, using worldview2 ndvi for comparison with colville models
l.ndvi.df<-subset(lonely.ndvi,select=c('NDVI_WV.2','Year','jdate','Plot','Treatment'))
l.pn.df<-subset(lonely.cn,select=c('X.N','Year','jdate','Plot'))
l.ndvi.df<-na.omit(l.ndvi.df)
l.pn.df<-na.omit(l.pn.df)
l.ndvi.df$DOY=l.ndvi.df$jdate
l.pn.df$DOY=l.pn.df$jdate
l.pn.df$PrcntN=l.pn.df$X.N

#excluding weird treatments (greenhouse/shadehouse)
#l.ndvi.df<-subset(l.ndvi.df,l.ndvi.df$Treatment %in% c('Control','Exclosure'))

#merging
l.ndvi.df=merge(l.ndvi.df,l.pn.df,by=c('Year','DOY','Plot','jdate'))
l.ndvi.df$Yr=paste(l.ndvi.df$Plot,l.ndvi.df$Year,sep='_') ###two years is not enough

####
#ndvi models
#using the best models from our site
l.sr<-ndvi.phenolvals.func(l.ndvi.df,'Yr','NDVI_WV.2',qua=F,trd=F,wide.doy=T)

#l.sr2<-ndvi.phenolpoly.func(l.ndvi.df,'Yr','NDVI_WV.2',3,wide.doy=F)
#can't use first obs day for this validation exercise.

##generating prediction based on ndvi emodis from colville
#coefs matrix created in ndvi_n_01192016.r

colville.pred<-function(l.df,coefs)
{
  maxder.coefs=coefs[1,] #coefficients copied from colville models
  halfmax.coefs=coefs[2,]
  maxval.coefs=coefs[3,]

  pred.maxder=maxder.coefs[1]+maxder.coefs[2]*l.df$NDVI_WV.2.doy.maxDeriv
  pred.50max=halfmax.coefs[1]+halfmax.coefs[2]*l.df$NDVI_WV.2.doy.50max
  pred.maxval=maxval.coefs[1]+maxval.coefs[2]*l.df$NDVI_WV.2.doy.maxVal

  l.df<-cbind(l.df,pred.maxder,pred.50max,pred.maxval)
  return(l.df)
}

l.sr=colville.pred(l.sr,coefs)

l.sr=l.sr[l.sr$doy.maxN>145,] #some pcnt N curves were too sparse to fit properly, so thought peak n was the first
#day of observations (may1) which is obviously

maxder.lm<-lm(doy.maxN~pred.maxder,data=l.sr)
halfmax.lm<-lm(doy.maxN~pred.50max,data=l.sr)
maxval.lm<-lm(doy.maxN~pred.maxval,data=l.sr)


#plotting results
par(mfrow=c(3,1))
par(oma = c(5,4,0,0) + 0.1,
    mar = c(2,0,2,1.5) + 0.1)
plot(doy.maxN~pred.maxder,data=l.sr,xlim=c(160,190),ylim=c(160,200),main='colville maxDer pred')
abline(a=0,b=1,lty=2)
lm_eqn(maxder.lm)

plot(doy.maxN~pred.50max,data=l.sr,xlim=c(160,190),ylim=c(160,200),main='colville 50max pred')
abline(a=0,b=1,lty=2)
lm_eqn(halfmax.lm)

plot(doy.maxN~pred.maxval,data=l.sr,xlim=c(160,190),ylim=c(160,200),main='colville doy.maxval pred')
abline(a=0,b=1,lty=2)
lm_eqn(maxval.lm)

title(xlab = "Obs DOY %N",
      ylab = "Pred DOY %N",
      outer = TRUE, line = 2)

names(l.sr)
sr.forbox=subset(l.sr,select=c('doy.maxN','pred.maxder','pred.50max','pred.maxval'))
sr.forbox<-melt(sr.forbox)
par(mfrow=c(1,1))
par(mar=c(3,4,4,3)+.1)
boxplot(value~variable,data=sr.forbox,main='predicted (colville SR_WVnir) vs. \nobs (lonely) doy max %N',ylab='DOY')
sr.forbox$variable=factor(sr.forbox$variable,levels=c('doy.maxN','pred.maxder','pred.50max','pred.maxval'))
summary(glm(value~variable,data=sr.forbox,family='quasi'))
abline(v=1.5,lty=2)

#####


#lonely.emodis validation
lonely.emo<-read.csv('../Data/lonely_emo_transformed.csv')
lonely.emo$Yr=lonely.emo$Year
lonely.emo<-subset(lonely.emo,lonely.emo$Single.Pixel==0,select=c('Year','Yr','DOY','eMO_NDVI','Plot'))
phenol_curve_pts(lonely.emo,'eMO_NDVI',qua=F,trd=F)
lonely.emo<-subset(lonely.emo,lonely.emo$Year==2012)

l.emo1<-phenol_curve_pts(lonely.emo[lonely.emo$Plot=='Pixel1',],'eMO_NDVI',qua=F,trd=F,wide.doy=T)
l.emo2<-phenol_curve_pts(lonely.emo[lonely.emo$Plot=='Pixel2',],'eMO_NDVI',qua=F,trd=F,wide.doy=T)
l.emo3<-phenol_curve_pts(lonely.emo[lonely.emo$Plot=='Pixel3',],'eMO_NDVI',qua=F,trd=F,wide.doy=T)

l.emo=rbind(l.emo1,l.emo2,l.emo3)

# l.emo<-phenol_curve_pts(lonely.emo,'eMO_NDVI',qua=F,trd=F)

colville.pred2<-function(l.emo,emo.coefs)
{
  maxder.coefs=emo.coefs[1,] #coefficients copied from colville models
  halfmax.coefs=emo.coefs[2,]
  maxval.coefs=emo.coefs[3,]
  
  l.emo$pred.maxder=maxder.coefs[1]+maxder.coefs[2]*l.emo$eMO_NDVI.doy.maxDeriv
  l.emo$pred.50max=halfmax.coefs[1]+halfmax.coefs[2]*l.emo$eMO_NDVI.doy.50max
  l.emo$pred.maxval=maxval.coefs[1]+maxval.coefs[2]*l.emo$eMO_NDVI.doy.maxVal
  
  return(l.emo)
}

l.emo<-colville.pred2(l.emo,emo.coefs)
boxplot(l.emo$pred.maxder)
boxplot(l.sr$doy.maxN,add=T,col='blue')
par(mfrow=c(1,1))
boxplot(l.sr$doy.maxN,ylim=c(170,200))
abline(h=l.emo$pred.maxder,lty=2,col='red')
abline(h=l.emo$pred.50max,lty=2,col='blue')
abline(h=l.emo$pred.maxval,lty=2,col='green')

for.box<-subset(l.emo,select=c(pred.maxder,pred.50max,pred.maxval))
for.box=melt(for.box)

maxn.obs=l.sr$doy.maxN
maxn.obs=data.frame(variable='doy.maxN',value=maxn.obs)

for.box=rbind(for.box,maxn.obs)
for.box$variable<-factor(for.box$variable,levels=c('doy.maxN','pred.maxder','pred.50max','pred.maxval'))

par(mfrow=c(1,1))
par(mar=c(3,4,4,3)+.1)
boxplot(value~variable,data=for.box,main='2012 Lonely doy max %N\nPredicted (eMODIS) vs. observed (doy.maxN)',ylab='DOY')
abline(v=1.5,lty=2)

summary(glm(value~variable,data=for.box,family='quasi'))
