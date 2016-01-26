#NDVI percent N plots and non-linear curve fits

#reading in the data
setwd('~/Data/')
dir()

#csv file cleaned up by Dan.
#need to compare to original data file- confirm that removing na's in kyle's analysis resulted in consistent dataset being used for all models compared w/ AIC.
vegplots<-read.csv('VegPlots_for R.csv')

#load phenology curve functions
source('../code/phenolcurve.r')
source('../code/aic_code.r')
library(AICcmodavg)
library(MASS) #fitdistr function in MASS allows one to fit good parameters.
library(minpack.lm)
library(ggplot2)

####################################################################
#ggplot scatter plots to look at variation in %N and ndvi over season and by year.
#install ggplot2 if necessary
vegplots$Yr.fac<-as.factor(vegplots$Yr) #making year a grouping factor.

#ndvi over season, by year
#- using worldview2 ndvi from specrad readings.
#bmp('../output_plots/emo_ndvi_season_year.bmp',height=480,width=640)
eMO.ndvi.plot<-ggplot(dat=vegplots, aes(x=DOY,y=eMO_NDVI,color=Yr.fac))+geom_point()+geom_smooth()
eMO.ndvi.plot
#dev.off()
#asymptote with gradual decline on the right side.

#% N over season, by Year
#bmp('../output_plots/pctN_season_year.bmp',height=480,width=640)
pN.plot<-ggplot(dat=vegplots, aes(x=DOY,y=PrcntN,color=Yr.fac))+geom_point()+geom_smooth()
pN.plot
#dev.off()
#looks like a right-tailed curve.