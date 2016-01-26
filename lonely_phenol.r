#NDVI percent N plots and non-linear curve fits
setwd('~/code/')
#reading in the data
source('julian.r')
source('phenolcurve1172015.r')
source('aic_code.r')

library(minpack.lm)
library(ggplot2)
library(reshape2)

load('../output_plots/ndvi.n.phenol01192016.R')
dir()

#csv file cleaned up by Dan.
#need to compare to original data file- confirm that removing na's in kyle's analysis resulted in consistent dataset being used for all models compared w/ AIC.
vegplots<-read.csv('../Data/VegPlots_for R.csv')

lonely.ndvi<-read.csv('../Data/Lonelybiomassndvi.csv')
lonely.cn<-read.csv('../Data/lonely_cn.csv')

lonely.emo<-read.csv('../Data/lonely_emo_transformed.csv')
lonely.emo$Yr=lonely.emo$Year
lonely.emo<-subset(lonely.emo,lonely.emo$Single.Pixel==1,select=c('Year','Yr','DOY','eMO_NDVI','Plot'))
phenol_curve_pts(lonely.emo,'eMO_NDVI',qua=F,trd=F)

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
l.ndvi.df<-subset(l.ndvi.df,l.ndvi.df$Treatment %in% c('Control','Exclosure'))

#merging
l.ndvi.df=merge(l.ndvi.df,l.pn.df,by=c('Year','DOY','Plot','jdate'))
l.ndvi.df$Yr=paste(l.ndvi.df$Plot,l.ndvi.df$Year,sep='_') ###two years is not enough

#making second dframe for merge with eMO_NDVI satellite dframe
l.ndvi.df2<-subset(l.ndvi.df,select=c('Year','Plot','DOY','PrcntN'))
l.ndvi.df2<-merge(l.ndvi.df2,lonely.emo,all.x=T,all.y=T)

####
#ndvi models
#using the best models from our site - sr wv2
l.wv<-ndvi.phenolvals.func(l.ndvi.df,'Yr','NDVI_WV.2',qua=F,trd=F,wide.doy=T)

#and again with emodis ndvi (satellite-derived, single pixel)
#l.wv<-ndvi.phenolvals.func(l.ndvi.df2,'Yr','eMO_NDV',qua=F,trd=F,wide.doy=T)


##generating prediction based on ndvi emodis from colville
#coefs matrix created in ndvi_n_01192016.r

colville.pred<-function(l.df,coefs)
{
  maxder.coefs=coefs[1,] #coefficients copied from colville models
  halfmax.coefs=coefs[2,]
  maxval.coefs=coefs[3,]

  l.df$pred.maxder=maxder.coefs[1]+maxder.coefs[2]*l.df$NDVI_WV.2.doy.maxDeriv
  l.df$pred.50max=halfmax.coefs[1]+halfmax.coefs[2]*l.df$NDVI_WV.2.doy.50max
  l.df$pred.maxval=maxval.coefs[1]+maxval.coefs[2]*l.df$NDVI_WV.2.doy.maxVal

  return(l.df)
}

l.wv=colville.pred(l.wv,sr.wv.coefs)

l.wv=l.wv[l.wv$doy.maxN>145,] #some pcnt N curves were too sparse to fit properly, so thought peak n was the first
#day of observations (may1) which is obviously

maxder.lm<-lm(doy.maxN~pred.maxder,data=l.wv)
halfmax.lm<-lm(doy.maxN~pred.50max,data=l.wv)
maxval.lm<-lm(doy.maxN~pred.maxval,data=l.wv)


#plotting results
par(mfrow=c(3,1))
par(oma = c(5,4,0,0) + 0.1,
    mar = c(2,0,2,1.5) + 0.1)
plot(doy.maxN~pred.maxder,data=l.wv,xlim=c(160,190),ylim=c(160,200),main='colville maxDer pred')
abline(a=0,b=1,lty=2)
lm_eqn(maxder.lm)

plot(doy.maxN~pred.50max,data=l.wv,xlim=c(160,190),ylim=c(160,200),main='colville 50max pred')
abline(a=0,b=1,lty=2)
lm_eqn(halfmax.lm)

plot(doy.maxN~pred.maxval,data=l.wv,xlim=c(160,190),ylim=c(160,200),main='colville doy.maxval pred')
abline(a=0,b=1,lty=2)
lm_eqn(maxval.lm)

title(xlab = "Obs DOY %N",
      ylab = "Pred DOY %N",
      outer = TRUE, line = 2)


###################################
#repeating the process with new emodis ndvi data from mike buddy.
