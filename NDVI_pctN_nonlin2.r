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
######################## Basic Scatterplots##############################################################
# #ggplot scatter plots to look at variation in %N and ndvi over season and by year.
# #install ggplot2 if necessary
# library(ggplot2)
# vegplots$Yr.fac<-as.factor(vegplots$Yr) #making year a grouping factor.
# 
# #ndvi over season, by year
# #- using worldview2 ndvi from specrad readings.
# pdf('../output_plots/ndvi_season_year.pdf',height=4,width=8)
# wv.sr.ndvi.plot<-ggplot(dat=vegplots, aes(x=DOY,y=SR_WVnir_NDVI,color=Yr.fac))+geom_point()+geom_smooth()+facet_grid(~Yr.fac)
# wv.sr.ndvi.plot
# dev.off()
# #asymptote with gradual decline on the right side.
# 
# #% N over season, by Year
# pdf('../output_plots/pctN_season_year.pdf',height=4,width=8)
# wv.sr.ndvi.plot<-ggplot(dat=vegplots, aes(x=DOY,y=PrcntN,color=Yr.fac))+geom_point()+geom_smooth()+facet_grid(~Yr.fac)
# wv.sr.ndvi.plot
# dev.off()
# #looks like a right-tailed curve.


############################################################################################
#calculate median values by year
cbind(names(vegplots))

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

########################################### Median Scatterplots with loess curves ###########################3
# #median.ndvi over season, by year
# pdf('../output_plots/median_ndvi_season_year.pdf',height=4,width=8)
# median.sr.ndvi.plot<-ggplot(dat=median.vals, aes(x=DOY,y=SR_WVnir_NDVI,color=Yr))+geom_point()+geom_smooth(se=FALSE)+facet_grid(~Yr)
# median.sr.ndvi.plot
# dev.off()
# 
# pdf('../output_plots/median_pctn_season_year.pdf',height=4,width=8)
# median.sr.ndvi.plot<-ggplot(dat=median.vals, aes(x=DOY,y=PrcntN,color=Yr))+geom_point()+geom_smooth(se=FALSE)+facet_grid(~Yr)
# median.sr.ndvi.plot
# dev.off()

################################################################################################################3
#phenology metric extraction!!
#first step is creating a subset database so you don't lose data to NA's in other columns
ndvi.df<-subset(vegplots,select=c(ndvi.names,'DOY','Yr','PrcntN'))
ndvi.df<-na.omit(ndvi.df)

poly.m<-thirdorder(dframe=ndvi.df,yvar='PrcntN')

pred.df=ndvi.df[unique(ndvi.df$DOY),]
pred.df=pred.df[order(pred.df$DOY),]

poly.pred=predict(poly.m,pred.df)

############################################################################
###nonlinear curves
plot(SR_WVnir_NDVI~DOY,data=ndvi.df)
#############################################################################
#gaussian process regression
gaus.f=function(yvar){return(formula(paste(yvar,"~ k*exp(-1/2*(DOY-mu)^2/sigma^2)")))}
emo.gaus=gaus.f('SR_WVnir_NDVI') #use mean x, sd x, and y height as starting pars
prn.gaus=gaus.f('PrcntN')

gaus.emo.nls<- nlsLM(emo.gaus, data = ndvi.df,start=c(mu=200,sigma=20,k=1))
gaus.pred=gaus.emo.nls$m$predict(pred.df)
lines(pred.df$DOY,gaus.pred,type='l')

plot(PrcntN~DOY,data=median.vals)
gaus.pn.nls=nlsLM(prn.gaus,data=median.vals,start=c(mu=200,sigma=20,k=1))
gaus.pn.pred=gaus.pn.nls$m$predict(pred.df)
lines(pred.df$DOY,gaus.pn.pred,col='red')
##################################################################
#lognormal process regression
lognormal=function(x, mu, sd, C) {
  ans <- vector(length = length(x), mode = "numeric")
  for (i in 1:length(x)) {
    if (x[i]>0) {
      value <- exp(-(log10(x[i])-mu)^2/(2*sd^2))/(sqrt(2*pi)*sd*x[i])*C
      ans[i] <- value
    }  else {  ans[i] <- 0  }  };  return(ans)  }

wvnir.nls=nlsLM(SR_WVnir_NDVI~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df)
logn.pred=wvnir.nls$m$predict(pred.df)
lines(pred.df$DOY,logn.pred,col='blue')

pn.logn.nls<-nlsLM(PrcntN~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df)
pn.pred=pn.logn.nls$m$predict(pred.df)

plot(PrcntN~DOY,data=ndvi.df)
lines(pred.df$DOY,pn.pred,col='blue')


#########################################################
library(nls2)
weibull.f=function(yvar){formula(paste(yvar,"~ a-b*exp(-c*DOY^d) + c",collapse=''))}
ndvi.w=weibull.f('SR_WVnir_NDVI')
pars <- expand.grid(a=seq(0.1, 300, len=10),
                    b=seq(0.1, 300, len=10),
                    c=seq(0.1,300,len=10),
                    d=seq(10, 300, len=10))

test=nlsLM(ndvi.w,data=ndvi.df,start=pars,algorithm='plinear-random')

#response ~ c * ((temp/a)^(b - 1)) * exp(-(temp/a)^b)
library(minpack.lm)
emo.w=weibull.f('eMO_NDVI')


plot(eMO_NDVI~DOY,data=ndvi.df)

lines(pred.df$DOY,poly.pred,col='red')
lines(pred.df$DOY,logn.pred,col='blue')
lines(pred.df$DOY,gaus.pred,col='green')

legend(x=200,y=.2,legend=c('peak','50%max','max Delta'),col=c('blue','red','red'),
       pch=c(15,16,17))

#fitting and comparing gaussian, weibull, and lognormal curves to percent N data as with doiron
#weibull:  f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a) for x > 0
#cumweibull:F(x) = 1 - exp(- (x/b)^a) on x > 0, the mean is E(X) = b ??(1 + 1/a), and the Var(X) = b^2 * (??(1 + 2/a) - (??(1 + 1/a))^2)
#lognormal: f(x) = 1/(???(2 ??) ?? x) e^-((log x - ??)^2 / (2 ??^2))
#gaussian: f(x) = a*exp(-((x-b)^2)/2*c^2)


 m1<-nlsLM(PrcntN~(c*((DOY/a)^(b-1)))*exp(-(DOY/a)^b),start=list(a=160,b=2,c=1.5),data=ndvi.df,control=nls.lm.control(maxiter=150))


xvar=pred.df$DOY

DOY=1:400
a=193.910;b=1.917;c=6.634
wei.pred=(c*((DOY/a)^(b-1)))*exp(-(DOY/a)^b)
plot(DOY,wei.pred)

############comparing polynomial power models
#can run a simple quadratic or third-order polynomial with quadratic(yvar,df) or thirdorder(yvar,df)

model.compare('PrcntN',ndvi.df)
wv1.curveaic=model.compare('SR_WVnir_NDVI',ndvi.df)
emo.curveaic=model.compare('eMO_NDVI',ndvi.df)
sr_emo.curveaic=model.compare('SR_eMO_NDVI',ndvi.df)
redege.curveaic=model.compare('SR_WVre_NDVI',ndvi.df)

#aictab(emo.curveaic) #model averaging poorly defined for gam- may need to write own function and/or reconsider use.

#aic tables from this output could be used for model averaging- 
#including gams in this process seems dicier, but might still be warranted- need to read more.
#####################################################################3

#gam & polynomial model of ndvi- trying to extract first derivative. using phenolcurve.r

#with polynomial curves in dotted lines.
N.phenol1=phenol_curve_pts(ndvi.df,'PrcntN')[,2]
N.phenol.p3=phenol_poly_pts(ndvi.df,'PrcntN',3,add=T)[,2]
N.phenol.p2=phenol_poly_pts(ndvi.df,'PrcntN',2,add=T)[,2]

ndvi.phenol2=phenol_curve_pts(ndvi.df,ndvi.names[2])
ndvi.phenol2=phenol_poly_pts(ndvi.df,ndvi.names[2],3,add=T)
ndvi.phenol2=phenol_poly_pts(ndvi.df,ndvi.names[2],2,add=T)

#now calculating these vals for all years separately
yearname='Yr'
sr1<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[1])
emo<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[2])
sr2<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[3])
sr3<-ndvi.phenolvals.func(ndvi.df,yearname,ndvi.names[4])

#will need to fix/modify this function so you don't have to fit the same polynomials to the ndvi and n datasets.
#polynomial.curves
sr.poly=ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[1],3)
emo.poly=ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[2],3)
sr.poly2=ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[3],3)
sr.poly3=ndvi.phenolpoly.func(ndvi.df,yearname,ndvi.names[4],3)

# #compare w/ second order
# sr.poly=phenol_poly_pts(ndvi.df,ndvi.names[1],2)
# emo.poly=phenol_poly_pts(ndvi.df,ndvi.names[2],2)
# sr.poly2=phenol_poly_pts(ndvi.df,ndvi.names[3],2)
# sr.poly3=phenol_poly_pts(ndvi.df,ndvi.names[4],2)

#based on third order polynomials
emo.poly<-pred.obs.compare(emo.poly,'doy.maxDeriv','doy.maxN','maxDer')
emo.poly.maxDer=predn.func(emo.poly,'doy.maxDeriv','doy.maxN')
emo.poly<-pred.obs.compare(emo.poly,'doy.50max','doy.maxN','50max')
emo.poly.50max=predn.func(emo.poly,'doy.50max','doy.maxN')

f#based on smooth spline gam
emo.maxDer.df<-pred.obs.compare(emo,'doy.maxDeriv','doy.maxN','maxDer')
emo.maxDer=predn.func(emo,'doy.maxDeriv','doy.maxN')
emo.50max.df<-pred.obs.compare(emo,'doy.50max','doy.maxN','50max')
emo.50max=predn.func(emo,'doy.50max','doy.maxN')

summary(emo.poly.maxDer$model) #not a true measure of predictive power- would need to predict against another dataset and/or do some kind of bootstrapping.
emo.poly.maxDer$coefs

summary(emo.poly.50max$model)
emo.poly.50max$coefs

summary(emo.50max$model)
summary(emo.maxDer$model)

#may be worth using the second order polynomial for ndvi if that works better.
#do a lit review on seasonal trends in ndvi and biomass

#repeating plots with median vals
sr.m1<-ndvi.phenolvals.func(median.vals,yearname,ndvi.names[1]) #sig relationship with maxval
emo.m1<-ndvi.phenolvals.func(median.vals,yearname,ndvi.names[2])
sr.m2<-ndvi.phenolvals.func(median.vals,yearname,ndvi.names[3]) #sig relationship with maxval
sr.m3<-ndvi.phenolvals.func(median.vals,yearname,ndvi.names[4])

emo.m1.poly=ndvi.phenolvals.func(median.vals,yearname,ndvi.names[2])









#todo
#improve resolution of gam function-done

#fit parameteric curves and extract derivative info from them
#compare your methods to kyles w/ third order polynomials- done.

#basic plots for other NDVI measures- overlay on scatterplot.-done

#compare performance of metrics across ndvi measures, and across
#randomforest model using emodis+weather/environment- try predicting lonely camp biomass with this
#rf pct n vs. ndvi model


