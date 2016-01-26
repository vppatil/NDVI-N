setwd('code/')
dir()

#csv file cleaned up by Dan.
#need to compare to original data file- confirm that removing na's in kyle's analysis resulted in consistent dataset being used for all models compared w/ AIC.

#will add code to curve functions to include a column with first day of observations

#Data
vegplots<-read.csv('../Data/VegPlots_for R.csv')

#functions
source('phenolcurve1172015.r') #fromchromebook
source('aic_code.r')
source('julian.r')

#packages
library(minpack.lm)
library(ggplot2)

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
median.vals$Year<-as.factor(median.vals$Yr)

####subsetting vegplots for analysis
ndvi.df<-subset(vegplots,select=c(ndvi.names,'DOY','Yr','PrcntN')) #making df with all ndvi vars
ndvi.df$Year<-as.factor(ndvi.df$Yr)
##############################
##scatterplots with ggplot to look at trends in raw data
#eMODIS_NDVI plotted based on median vals due to duplication of vals across plots
vegplots$Year<-as.factor(vegplots$Yr) #making year a grouping factor.

#ndvi over season, by year
#- using emodis ndvi.
bmp('../output_plots/emo_ndvi_season_year.bmp',height=480,width=640)
eMO.ndvi.plot<-ggplot(dat=median.vals, aes(x=DOY,y=eMO_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)
eMO.ndvi.plot
dev.off()

#- using emodis_calc ndvi from specrad readings.
bmp('../output_plots/sr_ emo_ndvi_season_year.bmp',height=480,width=640)
sr.eMO.ndvi.plot<-ggplot(dat=ndvi.df, aes(x=DOY,y=SR_eMO_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)
sr.eMO.ndvi.plot
dev.off()

#% N over season, by Year
bmp('../output_plots/pctn_season_year.bmp',height=480,width=640)
pn.plot<-ggplot(dat=vegplots, aes(x=DOY,y=PrcntN,color=Year))+geom_point()+geom_smooth(se=F)
pn.plot
dev.off()
###########################################end scatterplots

########model selection table to compare curve fits, data pooled over years
emo.aic<-phenol_models_aictab(median.vals,'eMO_NDVI',use.aicc=TRUE)
sr.aic<-phenol_models_aictab(ndvi.df,'SR_eMO_NDVI',use.aicc=TRUE)
pn.aic<-phenol_models_aictab(ndvi.df,'PrcntN',use.aicc=TRUE)

emo.aic$var='eMODIS'
sr.aic$var='SR_eMODIS'
pn.aic$var='PrcntN'

phenol.modeltable<-list(emo.aic,sr.aic,pn.aic)
phenol.modeltable
########################################end model selection table

###################################################################3
#curve fit plots, data pooled over all years.
phenol_curve_pts(dframe=median.vals,varname='eMO_NDVI',qua=F,trd=F,wide.doy=F)
legend(locator(1),pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))
phenol_curve_pts(dframe=ndvi.df,varname='SR_eMO_NDVI',qua=F,trd=F,wide.doy=T)
legend(locator(1),pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))
phenol_poly_pts(dframe=median.vals,varname='PrcntN',porder=3)
legend(locator(1),pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))

########################phenology index table generation.

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

sr.phenol.median=pred.obs.compare(sr.phenol.median,'SR_eMO_NDVI.doy.maxVal','doy.maxN','maxVal')
sr.phenol.maxval=predn.func(sr.phenol.median,'SR_eMO_NDVI.doy.maxVal','doy.maxN')
lm_eqn(sr.phenol.maxval$model)

summary(sr.phenol.maxDer$model)
summary(sr.phenol.50max$model)

#polynomial fit
sr.phenol2.median<-ndvi.phenolpoly.func(sr.df,yearname,ndvi.names[3],porder=3,wide.doy=T)
sr.phenol2.median<-pred.obs.compare(sr.phenol2.median,'doy.maxDeriv','doy.maxN','maxDer')
sr.phenol2.maxDer=predn.func(sr.phenol2.median,'doy.maxDeriv','doy.maxN')
lm_eqn(sr.phenol2.maxDer$model)

sr.phenol2.median<-pred.obs.compare(sr.phenol2.median,'doy.50max','doy.maxN','50max')
sr.phenol2.50max=predn.func(sr.phenol2.median,'doy.50max','doy.maxN')
lm_eqn(sr.phenol2.50max$model)

summary(sr.phenol2.maxDer$model)
summary(sr.phenol2.50max$model)

#####making 2x4 panel plot for pred vs. obs scattereplots
#x axis (2) = eMODIS_NDVI vs. SR_eMO_NDVI
#y axis (4) = maxDer, 50max, maxVal, 1st day obs

#xlab should be %N observed
#ylab = %N pred
#main = ndvi var + phenology index ('eMO_NDVI doy.50max')
#emodis plots first
emo.phenol.median<-ndvi.phenolvals.func(median.vals,yearname,ndvi.names[2],qua=F,trd=F,wide.doy=T)
emo.poly.median<-ndvi.phenolpoly.func(median.vals,yearname,ndvi.names[2],porder=3,wide.doy=F)

sr.phenol<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[3],qua=F,trd=F,wide.doy=T)
sr.poly<-ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[3],porder = 3,wide.doy = F)

plot.new()
par(mfrow=c(2,4))

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
emo.phenol.maxDer=predn.func(emo.phenol.median,'eMO_NDVI.doy.maxDeriv','doy.maxN')
lm_eqn(emo.phenol.maxDer$model)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.50max','doy.maxN','50max')
emo.phenol.50max=predn.func(emo.phenol.median,'eMO_NDVI.doy.50max','doy.maxN')
lm_eqn(emo.phenol.50max$model)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.maxVal','doy.maxN','maxVal')
emo.phenol.maxval=predn.func(emo.phenol.median,'eMO_NDVI.doy.maxVal','doy.maxN')
lm_eqn(emo.phenol.maxval$model)

emo.poly.median<-pred.obs.compare(emo.poly.median,'first.obs.day','doy.maxN','first.obs.day')
emo.poly.firstday=predn.func(emo.poly.median,'first.obs.day','doy.maxN')
lm_eqn(emo.poly.firstday$model)

#####SR emomdis
sr.phenol<-pred.obs.compare(sr.phenol,'SR_eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
emo.phenol.maxDer=predn.func(sr.phenol,'SR_eMO_NDVI.doy.maxDeriv','doy.maxN')
lm_eqn(emo.phenol.maxDer$model)

sr.phenol<-pred.obs.compare(sr.phenol,'SR_eMO_NDVI.doy.50max','doy.maxN','50max')
emo.phenol.50max=predn.func(sr.phenol,'SR_eMO_NDVI.doy.50max','doy.maxN')
lm_eqn(emo.phenol.50max$model)

sr.phenol<-pred.obs.compare(sr.phenol,'SR_eMO_NDVI.doy.maxVal','doy.maxN','maxVal')
emo.phenol.maxval=predn.func(sr.phenol,'SR_eMO_NDVI.doy.maxVal','doy.maxN')
lm_eqn(emo.phenol.maxval$model)

sr.poly<-pred.obs.compare(sr.poly,'first.obs.day','doy.maxN','first.obs.day')
emo.poly.firstday=predn.func(sr.poly,'first.obs.day','doy.maxN')
lm_eqn(emo.poly.firstday$model)


