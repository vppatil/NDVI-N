#NDVI percent N plots and non-linear curve fits with phenology metric calculation
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
bmp('../output_plots/emo_ndvi_season_year.bmp',height=480,width=640)
eMO.ndvi.plot<-ggplot(dat=vegplots, aes(x=DOY,y=eMO_NDVI,color=Yr.fac))+geom_point()+geom_smooth()
eMO.ndvi.plot
dev.off()
#asymptote with gradual decline on the right side.

#% N over season, by Year
bmp('../output_plots/pctn_season_year.bmp',height=480,width=640)
wv.sr.ndvi.plot<-ggplot(dat=vegplots, aes(x=DOY,y=PrcntN,color=Yr.fac))+geom_point()+geom_smooth()
wv.sr.ndvi.plot
dev.off()
#looks like a right-tailed curve.

######################################################################
#curve fitting and model selection for all years and by year.
ndvi.names=c('eMO_NDVI','SR_eMO_NDVI')
ndvi.df<-subset(vegplots,select=c(ndvi.names,'DOY','Yr','PrcntN'))
ndvi.df<-na.omit(ndvi.df)

###curve fitting and model averaging all years
#emo_ndvi
emo.quad=quadratic(dframe=ndvi.df,yvar='eMO_NDVI')
emo.third=thirdorder(dframe=ndvi.df,yvar='eMO_NDVI')
emo.gaus=nlsLM(gaus.f('eMO_NDVI'),data=ndvi.df,start=c(mu=200,sigma=20,k=1))
emo.logn=nlsLM(eMO_NDVI~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df) #lognormal

#specrad emo ndvi
sr.emo.quad=quadratic(dframe=ndvi.df,yvar='SR_eMO_NDVI')
sr.emo.third=thirdorder(dframe=ndvi.df,yvar='SR_eMO_NDVI')
sr.emo.gaus=nlsLM(gaus.f('SR_eMO_NDVI'),data=ndvi.df,start=c(mu=200,sigma=20,k=1))
sr.emo.logn=nlsLM(SR_eMO_NDVI~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df) #lognormal

#percent N
prcntN.quad=quadratic(dframe=ndvi.df,yvar='PrcntN')
prcntN.third=thirdorder(dframe=ndvi.df,yvar='PrcntN')
prcntN.gaus=nlsLM(gaus.f('PrcntN'),data=ndvi.df,start=c(mu=200,sigma=20,k=1))
prcntN.logn=nlsLM(PrcntN~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df) #lognormal

#model averaging
emo.list=list(emo.quad,emo.third,emo.gaus,emo.logn)
sr.emo.list=list(sr.emo.quad,sr.emo.third,sr.emo.gaus,sr.emo.logn)
prcntN.list=list(prcntN.quad,prcntN.third,prcntN.gaus,prcntN.logn)

emo.aic<-AIC(emo.quad,emo.third,emo.gaus,emo.logn)
sr.emo.aic=AIC(sr.emo.quad,sr.emo.third,sr.emo.gaus,sr.emo.logn)
prcntN.aic=AIC(prcntN.quad,prcntN.third,prcntN.gaus,prcntN.logn)

emo.aic$n=sapply(emo.list,nobs)
sr.emo.aic$n=sapply(sr.emo.list,nobs)
prcntN.aic$n=sapply(prcntN.list,nobs)

emo.aic=aic.table(emo.aic,'AIC',k='df',n='n',aicc=T)
sr.emo.aic=aic.table(sr.emo.aic,'AIC',k='df',n='n',aicc=T)
prcntN.aic=aic.table(prcntN.aic,'AIC',k='df',n='n',aicc=T)

print(list(emo.aic,sr.emo.aic,prcntN.aic))

#averaged predictions
pred=data.frame(DOY=seq(min(ndvi.df$DOY),max(ndvi.df$DOY),by=.1)) #dataframe of 1/10 day increments

#emodis predictions
pred$eMO_NDVI=model.avg.prediction(emo.list,emo.aic$w,pred)
pred$emo.quad=predict(emo.list[[1]],pred)
pred$emo.third=predict(emo.list[[2]],pred)
pred$emo.gaus=predict(emo.list[[3]],pred)
pred$emo.logn=predict(emo.list[[4]],pred)

#sr_emodis
pred$SR_eMO_NDVI=model.avg.prediction(sr.emo.list,sr.emo.aic$w,pred)
pred$sr.emo.quad=predict(sr.emo.list[[1]],pred)
pred$sr.emo.third=predict(sr.emo.list[[2]],pred)
pred$sr.emo.gaus=predict(sr.emo.list[[3]],pred)
pred$sr.emo.logn=predict(sr.emo.list[[4]],pred)

#pN
pred$PrcntN=model.avg.prediction(prcntN.list,prcntN.aic$w,pred)
pred$prcntN.quad=predict(prcntN.list[[1]],pred)
pred$prcntN.third=predict(prcntN.list[[2]],pred)
pred$prcntN.gaus=predict(prcntN.list[[3]],pred)
pred$prcntN.logn=predict(prcntN.list[[4]],pred)

#pred dataframe now contains predicted values for both ndvi vars and percent n based on model avg and each 
#individual model

#############################
#plotting model-averaged curve along with other curves, with predicted phenology metrics shown
bmp('../output_plots/emo_mavg.bmp',height=480,width=640)
plot(eMO_NDVI~DOY,data=ndvi.df)
lines(eMO_NDVI~DOY,data=pred,lwd=2)
	emo.phenol=phenol_curve_pts(ndvi.df,'eMO_NDVI')
	emo.phenol.df=data.frame(DOY=t(emo.phenol),index=names(emo.phenol))
	emo.phenol.df$yval=pred$eMO_NDVI[match(emo.phenol.df$DOY,pred$DOY)]
	points(emo.phenol.df$DOY[1],emo.phenol.df$yval[1],col='red',pch=15,lwd=3,cex=2)
	points(emo.phenol.df$DOY[2],emo.phenol.df$yval[2],col='green',pch=16,lwd=3,cex=2)
	points(emo.phenol.df$DOY[3],emo.phenol.df$yval[3],col='blue',pch=17,lwd=3,cex=2)

	legend(x=200,y=.2,legend=c('peak','50%max','max Delta'),col=c('red','green','blue'),
       pch=c(15,16,17))
dev.off()

bmp('../output_plots/sr.emo_mavg.bmp',height=480,width=640)
plot(SR_eMO_NDVI~DOY,data=ndvi.df)
lines(SR_eMO_NDVI~DOY,data=pred,lwd=2)
sr.emo.phenol=phenol_curve_pts(ndvi.df,'SR_eMO_NDVI')
	sr.emo.phenol.df=data.frame(DOY=t(sr.emo.phenol),index=names(sr.emo.phenol))
	sr.emo.phenol.df$yval=pred$SR_eMO_NDVI[match(sr.emo.phenol.df$DOY,pred$DOY)]
	points(sr.emo.phenol.df$DOY[1],sr.emo.phenol.df$yval[1],col='red',pch=15,lwd=3,cex=2)
	points(sr.emo.phenol.df$DOY[2],sr.emo.phenol.df$yval[2],col='green',pch=16,lwd=3,cex=2)
	points(sr.emo.phenol.df$DOY[3],sr.emo.phenol.df$yval[3],col='blue',pch=17,lwd=3,cex=2)

	legend(x=205,y=.33,legend=c('peak','50%max','max Delta'),col=c('red','green','blue'),
       pch=c(15,16,17),bty='n')
dev.off()

bmp('../output_plots/pN_mavg.bmp',height=480,width=640)
plot(PrcntN~DOY,data=ndvi.df,ylab= '% N')
lines(prcntN.third~DOY,data=pred,lwd=2)
prcntN.phenol=phenol_curve_pts(ndvi.df,'PrcntN')
	prcntN.phenol.df=data.frame(DOY=t(prcntN.phenol),index=names(prcntN.phenol))
	prcntN.phenol.df$yval=pred$PrcntN[match(prcntN.phenol.df$DOY,pred$DOY)]
	points(prcntN.phenol.df$DOY[1],prcntN.phenol.df$yval[1],col='red',pch=15,lwd=3,cex=2)
	points(prcntN.phenol.df$DOY[3],prcntN.phenol.df$yval[3],col='blue',pch=17,lwd=3,cex=2)
	points(prcntN.phenol.df$DOY[2],prcntN.phenol.df$yval[2],col='green',pch=16,lwd=3,cex=2)

	legend(x=200,y=4.5,legend=c('peak','50%max','max Delta'),col=c('red','green','blue'),
       pch=c(15,16,17))
dev.off()

#################################################
#curve-fitting by years and phenology metric extraction
yrs=unique(ndvi.df$Yr)
phenol.vals=data.frame()
for(i in 1:length(yrs))
{
	temp=ndvi.df[ndvi.df$Yr==yrs[i],]
	
	emo.phenol=phenol_curve_pts(temp,'eMO_NDVI',modavg=FALSE)
	sr.emo.phenol=phenol_curve_pts(temp,'SR_eMO_NDVI',modavg=FALSE)
	pn.phenol=phenol_curve_pts(temp,'PrcntN',modavg=FALSE)
	
	phenol.yr<- cbind(emo.phenol,sr.emo.phenol,pn.phenol,Year=yrs[i])
	phenol.vals=rbind(phenol.vals,phenol.yr)
}
print(phenol.vals)
write.csv(phenol.vals,'../output_plots/phenolvals_table.csv')

############################pred_vs_obs plots (bootstrapped SE?)
pred.vars=names(phenol.vals)[1:6]
obs.n=phenol.vals$PrcntN.doy.maxVal


par(mfrow=c(3,2))
pred.obs.mlist=list()
for(i in 1:length(pred.vars))
{
	pvar=pred.vars[i]
	phenol.lm=lm(obs.n~phenol.vals[[pvar]])
	pred.obs.mlist[[i]]=phenol.lm
	names(pred.obs.mlist)[i]=pvar
	plot(obs.n~phenol.vals[[pvar]],xlab=paste('pred',pvar))
	abline(a=0,b=1,lty=2)
}

r2=sapply(pred.obs.mlist,function(x) summary(x)$adj.r.squared)
p=sapply(pred.obs.mlist,function(x) {fstat=summary(x)$fstatistic;pf(fstat[1],fstat[2],fstat[3],lower.tail=F)})
f=sapply(pred.obs.mlist,function(x) summary(x)$fstatistic)
coefs=sapply(pred.obs.mlist,coefficients)
row.names(coefs)[2]='slope'

pred.obs.summary<-rbind(r2,p,f,coefs)











##################################Lonely camp validation set.









###################################Environmental comparison table

