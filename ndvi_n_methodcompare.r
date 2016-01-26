#setwd('akscienceenterfiles/aksciencecenter files from drive/Data/')

setwd('code/')
dir()

#csv file cleaned up by Dan.
#need to compare to original data file- confirm that removing na's in kyle's analysis resulted in consistent dataset being used for all models compared w/ AIC.

#will add code to curve functions to include a column with first day of observations

vegplots<-read.csv('../Data/VegPlots_for R.csv')
source('phenolcurve1172015.r') #fromchromebook
source('aic_code.r')
source('julian.r')
library(minpack.lm)

cbind(names(vegplots))

yearname='Yr'
ndvi.names=c('SR_WVnir_NDVI','eMO_NDVI','SR_eMO_NDVI','SR_WVre_NDVI')
median.ndvi1<-aggregate(SR_WVnir_NDVI~DOY+Yr,data=vegplots,FUN=median)
median.ndvi2<-aggregate(eMO_NDVI~DOY+Yr,data=vegplots,FUN=median)
median.ndvi3<-aggregate(SR_eMO_NDVI~DOY+Yr,data=vegplots,FUN=median)
median.ndvi4<-aggregate(SR_WVre_NDVI~DOY+Yr,data=vegplots,FUN=median)

median.pctn<-aggregate(PrcntN~DOY+Yr,data=vegplots,FUN=median)

medians<-list(median.ndvi1,median.ndvi2,median.ndvi3,median.ndvi4,median.pctn)

merge.all <- function(x, y) {
  merge(x, y, all=TRUE, by=c('DOY','Yr'))
}

median.vals <- Reduce(merge.all, medians)
median.vals$Yr<-as.factor(median.vals$Yr)

##############################
ndvi.df<-subset(vegplots,select=c(ndvi.names,'DOY','Yr','PrcntN')) #making df with all ndvi vars

#old method using median vals dataframe
#emo.poly=ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[2],3)
emo.poly.median.obsdoy<-ndvi.phenolpoly.func(median.vals,yearname,ndvi.names[2],3,wide.doy=F)
emo.poly.median<-ndvi.phenolpoly.func(median.vals,yearname,ndvi.names[2],3,wide.doy=T)

emo.poly.median.obsdoy<-pred.obs.compare(emo.poly.median.obsdoy,'doy.maxDeriv','doy.maxN','maxDer')
emo.poly.maxDer=predn.func(emo.poly.median.obsdoy,'doy.maxDeriv','doy.maxN') #not looking at maxder
lm_eqn(emo.poly.maxDer$model)

emo.poly.median<-pred.obs.compare(emo.poly.median,'doy.50max','doy.maxN','50max')
emo.poly.50max=predn.func(emo.poly.median,'doy.50max','doy.maxN')
lm_eqn(emo.poly.50max$model)

summary(emo.poly.maxDer$model) #showing correspondence with first day of obs
summary(emo.poly.50max$model)
###################################
##comparing the third order models from new methods
#widened doy.pred range.
#gaussian and logn only for ndvi
#poly 3rd for prcnt N

emo.poly2.median<-ndvi.phenolvals.func(median.vals,yearname,ndvi.names[2],qua=F,trd=F,wide.doy=T)
emo.poly2.median<-pred.obs.compare(emo.poly2.median,'eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
emo.poly2.maxDer=predn.func(emo.poly2.median,'eMO_NDVI.doy.maxDeriv','doy.maxN')
lm_eqn(emo.poly2.maxDer$model)

emo.poly2.median<-pred.obs.compare(emo.poly2.median,'eMO_NDVI.doy.50max','doy.maxN','50max')
emo.poly2.50max=predn.func(emo.poly2.median,'eMO_NDVI.doy.50max','doy.maxN')
lm_eqn(emo.poly2.50max$model)
summary(emo.poly2.maxDer$model)
summary(emo.poly2.50max$model)

#a bit of debugging, confirm that emodis vals are not duplicates for each plot/doy/yr
#extract metric estimates, use model coefficients from best model (modav,gau lgn,emodis)
#to predict date of maxN, see how well they do.

#if duplicates use best sr emodis model

#now do for sr_eMO_NDVI
sr.df<-subset(ndvi.df,select=c('SR_eMO_NDVI','DOY','Yr','PrcntN'))
sr.phenol.median<-ndvi.phenolvals.func(sr.df,yearname,ndvi.names[3],qua=F,trd=F,wide.doy=T)
sr.phenol.median<-pred.obs.compare(sr.phenol.median,'SR_eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
sr.phenol.maxDer=predn.func(sr.phenol.median,'SR_eMO_NDVI.doy.maxDeriv','doy.maxN')
lm_eqn(sr.phenol.maxDer$model)

sr.phenol.median<-pred.obs.compare(sr.phenol.median,'SR_eMO_NDVI.doy.50max','doy.maxN','50max')
sr.phenol.50max=predn.func(sr.phenol.median,'SR_eMO_NDVI.doy.50max','doy.maxN')
lm_eqn(sr.phenol.50max$model)
summary(sr.phenol.maxDer$model)
summary(sr.phenol.50max$model)

#polynomial fit
sr.phenol2.median<-ndvi.phenolpoly.func(sr.df,yearname,ndvi.names[3],porder=3,wide.doy=F)
sr.phenol2.median<-pred.obs.compare(sr.phenol2.median,'doy.maxDeriv','doy.maxN','maxDer')
sr.phenol2.maxDer=predn.func(sr.phenol2.median,'doy.maxDeriv','doy.maxN')
lm_eqn(sr.phenol2.maxDer$model)

sr.phenol2.median<-pred.obs.compare(sr.phenol2.median,'doy.50max','doy.maxN','50max')
sr.phenol2.50max=predn.func(sr.phenol2.median,'doy.50max','doy.maxN')
lm_eqn(sr.phenol2.50max$model)

summary(sr.phenol2.maxDer$model)
summary(sr.phenol2.50max$model)

