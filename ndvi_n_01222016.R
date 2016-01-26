setwd('~/code/')
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
library(grid)
library(gridExtra)

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
#bmp('../output_plots/emo_ndvi_season_year.bmp',height=480,width=640)
eMO.ndvi.plot<-ggplot(dat=median.vals, aes(x=DOY,y=eMO_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)
#dev.off()

#- using emodis_calc ndvi from specrad readings.
#bmp('../output_plots/sr_ emo_ndvi_season_year.bmp',height=480,width=640)
sr.eMO.ndvi.plot<-ggplot(dat=ndvi.df, aes(x=DOY,y=SR_eMO_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)
#dev.off()

#- using wv2_calc ndvi from specrad readings.
#bmp('../output_plots/sr_ wv_ndvi_season_year.bmp',height=480,width=640)
sr.wv.ndvi.plot<-ggplot(dat=ndvi.df, aes(x=DOY,y=SR_WVnir_NDVI,color=Year))+geom_point()+geom_smooth(se=FALSE)

#% N over season, by Year
#bmp('../output_plots/pctn_season_year.bmp',height=480,width=640)
pn.plot<-ggplot(dat=vegplots, aes(x=DOY,y=PrcntN,color=Year))+geom_point()+geom_smooth(se=F)
#dev.off()

grid.arrange(eMO.ndvi.plot,sr.eMO.ndvi.plot,sr.wv.ndvi.plot,pn.plot,nrow=2)
###########################################end scatterplots

########model selection table to compare curve fits, data pooled over years
emo.aic<-phenol_models_aictab(median.vals,'eMO_NDVI',use.aicc=TRUE)
sr.aic<-phenol_models_aictab(ndvi.df,'SR_eMO_NDVI',use.aicc=TRUE)
sr.wv.aic<-phenol_models_aictab(ndvi.df,'SR_WVnir_NDVI',use.aicc=TRUE)

pn.aic<-phenol_models_aictab(ndvi.df,'PrcntN',use.aicc=TRUE)

emo.aic$var='eMODIS'
sr.aic$var='SR_eMODIS'
sr.wv.aic$var='SR__WVnir'

pn.aic$var='PrcntN'

phenol.modeltable<-list(emo.aic,sr.aic,sr.wv.aic,pn.aic)
phenol.modeltable
########################################end model selection table
###################################################################3
#curve fit plots, data pooled over all years.
par(mfrow=c(2,2))
par(mar=c(5,4,0.5,.5)+.1)
phenol_curve_pts(dframe=median.vals,varname='eMO_NDVI',qua=F,trd=F,wide.doy=F)
legend(195,.23,pch=c(15,16,17),bty='n',col=c('red','green','blue'),legend=c('maxVal','50% max','maxDeriv'))
phenol_curve_pts(dframe=ndvi.df,varname='SR_eMO_NDVI',qua=F,trd=F,wide.doy=T)
phenol_curve_pts(dframe=ndvi.df,varname='SR_WVnir_NDVI',qua=F,trd=F,wide.doy=T)
phenol_poly_pts(dframe=median.vals,varname='PrcntN',porder=3)

########################phenology index table generation.
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

sr.wv.phenol<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[1],qua=F,trd=F,wide.doy=T)
sr.wv.poly<-ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[1],porder = 3,wide.doy = F)


####multi-panel scatter showing predictive performance of various ndvi var and phenol indices
#bmp('../output_plots/pred_obs_plots.bmp',width = 480,height = 640)
#panel plot layout
layout.mat=matrix(1:9,3,3)
widths=c(2,2,2)
heights=c(4,4,4)
layout(layout.mat,widths=widths,heights=heights)

par(oma = c(5,4,0,0) + 0.1,
    mar = c(0,.5,1.5,1) + 0.1)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.maxDeriv','doy.maxN','maxDer')
emo.phenol.maxDer=predn.func(emo.phenol.median,'eMO_NDVI.doy.maxDeriv','doy.maxN',xaxt='n')
lm_eqn(emo.phenol.maxDer$model)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.50max','doy.maxN','50max')
emo.phenol.50max=predn.func(emo.phenol.median,'eMO_NDVI.doy.50max','doy.maxN',xaxt='n')
lm_eqn(emo.phenol.50max$model)

emo.phenol.median<-pred.obs.compare(emo.phenol.median,'eMO_NDVI.doy.maxVal','doy.maxN','maxVal')
emo.phenol.maxval=predn.func(emo.phenol.median,'eMO_NDVI.doy.maxVal','doy.maxN',xaxt='n')
lm_eqn(emo.phenol.maxval$model)

# emo.poly.median<-pred.obs.compare(emo.poly.median,'first.obs.day','doy.maxN','first.obs.day')
# emo.poly.firstday=predn.func(emo.poly.median,'first.obs.day','doy.maxN',ndvi='eMO_NDVI')
# lm_eqn(emo.poly.firstday$model)

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

# sr.poly<-pred.obs.compare(sr.poly,'first.obs.day','doy.maxN','first.obs.day')
# sr.poly.firstday=predn.func(sr.poly,'first.obs.day','doy.maxN',ndvi='SR_eMO_NDVI',yaxt='n')
# lm_eqn(sr.poly.firstday$model)

####SR worldview2
sr.wv.phenol<-pred.obs.compare(sr.wv.phenol,'SR_WVnir_NDVI.doy.maxDeriv','doy.maxN','maxDer')
sr.wv.phenol.maxDer=predn.func(sr.wv.phenol,'SR_WVnir_NDVI.doy.maxDeriv','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.wv.phenol.maxDer$model)

sr.wv.phenol<-pred.obs.compare(sr.wv.phenol,'SR_WVnir_NDVI.doy.50max','doy.maxN','50max')
sr.wv.phenol.50max=predn.func(sr.wv.phenol,'SR_WVnir_NDVI.doy.50max','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.wv.phenol.50max$model)

sr.wv.phenol<-pred.obs.compare(sr.wv.phenol,'SR_WVnir_NDVI.doy.maxVal','doy.maxN','maxVal')
sr.wv.phenol.maxval=predn.func(sr.wv.phenol,'SR_WVnir_NDVI.doy.maxVal','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.wv.phenol.maxval$model)

# sr.wv.poly<-pred.obs.compare(sr.wv.poly,'first.obs.day','doy.maxN','first.obs.day')
# sr.wv.poly.firstday=predn.func(sr.wv.poly,'first.obs.day','doy.maxN',ndvi='SR_WVnir_NDVI',yaxt='n')
# lm_eqn(sr.wv.poly.firstday$model)

title(xlab = "Obs DOY %N",
      ylab = "Pred DOY %N",
      outer = TRUE, line = 3,cex=3)
#dev.off()


###SR red edge wv
sr.re.phenol<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[4],qua=F,trd=F,wide.doy=T)
sr.re.poly<-ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[4],porder = 3,wide.doy = F)

layout.mat=matrix(1:4,2,2)
widths=c(2,2)
heights=c(4,4)
layout(layout.mat,widths=widths,heights=heights)

par(oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1.5,1) + 0.1)
sr.re.phenol<-pred.obs.compare(sr.re.phenol,'SR_WVre_NDVI.doy.maxDeriv','doy.maxN','maxDer')
sr.re.phenol.maxDer=predn.func(sr.re.phenol,'SR_WVre_NDVI.doy.maxDeriv','doy.maxN',xaxt='n')
lm_eqn(sr.re.phenol.maxDer$model)

sr.re.phenol<-pred.obs.compare(sr.re.phenol,'SR_WVre_NDVI.doy.50max','doy.maxN','50max')
sr.re.phenol.50max=predn.func(sr.re.phenol,'SR_WVre_NDVI.doy.50max','doy.maxN')
lm_eqn(sr.re.phenol.50max$model)

sr.re.phenol<-pred.obs.compare(sr.re.phenol,'SR_WVre_NDVI.doy.maxVal','doy.maxN','maxVal')
sr.re.phenol.maxval=predn.func(sr.re.phenol,'SR_WVre_NDVI.doy.maxVal','doy.maxN',xaxt='n',yaxt='n')
lm_eqn(sr.re.phenol.maxval$model)

sr.re.poly<-pred.obs.compare(sr.re.poly,'first.obs.day','doy.maxN','first.obs.day')
sr.re.poly.firstday=predn.func(sr.re.poly,'first.obs.day','doy.maxN',ndvi='SR_WVre_NDVI',yaxt='n')
lm_eqn(sr.re.poly.firstday$model)

title(xlab = "Obs DOY %N",
      ylab = "Pred DOY %N",
      outer = TRUE, line = 3)

######################################################
#tables
emo.phenol.median
sr.phenol
sr.wv.phenol

write.csv(list(emo.phenol.median,sr.phenol,sr.wv.phenol),'../output_plots/phenol_metric_tables.csv')

sr.re.phenol

# #coefficients for lonelycamp

#now using sr wv2 for comparison with lonely camp because it did the best overall of the two SR ndvi variables we have
#meaning that two of its phenol metrics predicted colville date of peak %N with reasonable R2 and p <.1
maxDer.coefs=coefficients(sr.wv.phenol.maxDer$model)
halfmax.coefs=coefficients(sr.wv.phenol.50max$model)
maxVal.coefs=coefficients(sr.wv.phenol.maxval$model)

sr.wv.coefs<-rbind(maxDer.coefs,halfmax.coefs,maxVal.coefs)
colnames(sr.wv.coefs)=c('Int','Slope')
row.names(sr.wv.coefs)=c('maxDer','halfmax','maxval')

#now emodis ndvi coefficients
maxDer.coefs=coefficients(emo.phenol.maxDer$model)
halfmax.coefs=coefficients(emo.phenol.50max$model)
maxVal.coefs=coefficients(emo.phenol.maxval$model)

emo.coefs<-rbind(maxDer.coefs,halfmax.coefs,maxVal.coefs)
colnames(emo.coefs)=c('Int','Slope')
row.names(emo.coefs)=c('maxDer','halfmax','maxval')
save.image('../output_plots/ndvi.n.phenol01222016.R')
