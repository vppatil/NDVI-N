setwd('code/')
dir()

#csv file cleaned up by Dan.
#need to compare to original data file- confirm that removing na's in kyle's analysis resulted in consistent dataset being used for all models compared w/ AIC.

#will add code to curve functions to include a column with first day of observations

#Data
vegplots<-read.csv('../Data/VegPlots_for R.csv')
emo.long<-read.csv('../Data/eMODIS_long.csv')
emo.long<-subset(emo.long,select=c('eMO_NDVI','Yr','DOY'))


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

vegplots$Year<-as.factor(vegplots$Yr) #making year a grouping factor.

####subsetting vegplots for analysis
ndvi.df<-subset(vegplots,select=c(ndvi.names[c(1,3,4)],'DOY','Yr','PrcntN')) #making df with all ndvi vars
ndvi.df<-merge(ndvi.df,emo.long,all.x=T,all.y=T)
ndvi.df$Year<-as.factor(ndvi.df$Yr)
##############################
##scatterplots with ggplot to look at trends in raw data
#eMODIS_NDVI plotted based on median vals due to duplication of vals across plots

#ndvi over season, by year
#- using emodis ndvi.
#bmp('../output_plots/emo_ndvi_season_year.bmp',height=480,width=640)
eMO.ndvi.plot<-ggplot(dat=ndvi.df, aes(x=DOY,y=eMO_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)
eMO.ndvi.plot
#dev.off()

#- using emodis_calc ndvi from specrad readings.
#bmp('../output_plots/sr_ emo_ndvi_season_year.bmp',height=480,width=640)
sr.eMO.ndvi.plot<-ggplot(dat=vegplots, aes(x=DOY,y=SR_eMO_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)
sr.eMO.ndvi.plot
#dev.off()

#% N over season, by Year
#bmp('../output_plots/pctn_season_year.bmp',height=480,width=640)
pn.plot<-ggplot(dat=vegplots, aes(x=DOY,y=PrcntN,color=Year))+geom_point()+geom_smooth(se=F)
pn.plot
#dev.off()
###########################################end scatterplots

########model selection table to compare curve fits, data pooled over years
emo.aic<-phenol_models_aictab(emo.long,'eMO_NDVI',use.aicc=TRUE)
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
phenol_curve_pts(dframe=emo.long,varname='eMO_NDVI',qua=F,trd=F,wide.doy=F)
legend(locator(1),pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))
phenol_curve_pts(dframe=ndvi.df,varname='SR_eMO_NDVI',qua=F,trd=F,wide.doy=T)
legend(locator(1),pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))
phenol_poly_pts(dframe=median.vals,varname='PrcntN',porder=3)
legend(locator(1),pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))

########################phenology index table generation.


#####making 2x4 panel plot for pred vs. obs scattereplots
#x axis (2) = eMODIS_NDVI vs. SR_eMO_NDVI
#y axis (4) = maxDer, 50max, maxVal, 1st day obs

#xlab should be %N observed
#ylab = %N pred
#main = ndvi var + phenology index ('eMO_NDVI doy.50max')
#emodis plots first
emo.phenol.median<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[2],qua=F,trd=F,wide.doy=F)
emo.poly.median<-ndvi.phenolpoly.func(median.vals,yearname,ndvi.names[2],porder=3,wide.doy=F)

sr.phenol<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[3],qua=F,trd=F,wide.doy=T)
sr.poly<-ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[3],porder = 3,wide.doy = F)


####multi-panel scatter showing predictive performance of various ndvi var and phenol indices
#bmp('../output_plots/pred_obs_plots.bmp',width = 480,height = 640)
#panel plot layout
layout.mat=matrix(1:8,4,2)
widths=c(2,2)
heights=c(4,4,4,4)
layout(layout.mat,widths=widths,heights=heights)

par(oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1.5,1) + 0.1)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
emo.phenol.maxDer=predn.func(emo.phenol.median,'eMO_NDVI.doy.maxDeriv','doy.maxN',xaxt='n')
lm_eqn(emo.phenol.maxDer$model)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.50max','doy.maxN','50max')
emo.phenol.50max=predn.func(emo.phenol.median,'eMO_NDVI.doy.50max','doy.maxN',xaxt='n')
lm_eqn(emo.phenol.50max$model)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.maxVal','doy.maxN','maxVal')
emo.phenol.maxval=predn.func(emo.phenol.median,'eMO_NDVI.doy.maxVal','doy.maxN',xaxt='n')
lm_eqn(emo.phenol.maxval$model)

emo.poly.median<-pred.obs.compare(emo.poly.median,'first.obs.day','doy.maxN','first.obs.day')
emo.poly.firstday=predn.func(emo.poly.median,'first.obs.day','doy.maxN',ndvi='eMO_NDVI')
lm_eqn(emo.poly.firstday$model)

#####SR emomdis
sr.phenol<-pred.obs.compare(sr.phenol,'SR_eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
sr.phenol.maxDer=predn.func(sr.phenol,'SR_eMO_NDVI.doy.maxDeriv','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.phenol.maxDer$model)

sr.phenol<-pred.obs.compare(sr.phenol,'SR_eMO_NDVI.doy.50max','doy.maxN','50max')
sr.phenol.50max=predn.func(sr.phenol,'SR_eMO_NDVI.doy.50max','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.phenol.50max$model)

sr.phenol<-pred.obs.compare(sr.phenol,'SR_eMO_NDVI.doy.maxVal','doy.maxN','maxVal')
sr.phenol.maxval=predn.func(sr.phenol,'SR_eMO_NDVI.doy.maxVal','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.phenol.maxval$model)

sr.poly<-pred.obs.compare(sr.poly,'first.obs.day','doy.maxN','first.obs.day')
sr.poly.firstday=predn.func(sr.poly,'first.obs.day','doy.maxN',ndvi='SR_eMO_NDVI',yaxt='n')
lm_eqn(sr.poly.firstday$model)

title(xlab = "Obs DOY %N",
      ylab = "Pred DOY %N",
      outer = TRUE, line = 3)
#dev.off()

#tables
emo.phenol.median
sr.phenol

#coefficients for lonelycamp
maxDer.coefs=coefficients(sr.phenol.maxDer$model)
halfmax.coefs=coefficients(sr.phenol.50max$model)
maxVal.coefs=coefficients(sr.phenol.maxval$model)
firstday.coefs=coefficients(sr.poly.firstday$model)

coefs<-rbind(maxDer.coefs,halfmax.coefs,maxVal.coefs,firstday.coefs)
colnames(coefs)=c('Int','Slope')
row.names(coefs)=c('maxDer','halfmax','maxval','firstday')
save.image('../output_plots/ndvi.n.phenol01192016.R')
