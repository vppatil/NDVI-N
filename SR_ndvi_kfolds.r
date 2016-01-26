setwd('~/code/')

#Data
vegplots<-read.csv('../Data/VegPlots_for R.csv')

#functions
source('phenolcurve1172015.r') #fromchromebook
source('aic_code.r')
source('julian.r')

#packages
library(minpack.lm)

yearname='Yr'
ndvi.df<-subset(vegplots,select=c('SR_WVre_NDVI','DOY','Yr','PrcntN')) #making df with ndvi var
ndvi.df<-na.omit(ndvi.df)

n=dim(ndvi.df)[1]

inds=1:n #indices for selecting training and test 

K=5
iter=1000

slice.size=round(n/K)

res=matrix(0,iter,3)

for(i in 1:iter)
{
	test.ind=sample(inds,slice.size,replace=F)
	
	ndvi.test=ndvi.df[test.ind,]
	ndvi.train=ndvi.df[-test.ind,]

	train.phenol=try(ndvi.phenolvals.func(ndvi.train,yearname,'SR_WVre_NDVI',qua=F,trd=F,wide.doy=T,plt=F),silent=T)
	if(class(train.phenol)=='try-error'){next()}
	test.phenol=try(ndvi.phenolvals.func(ndvi.test,yearname,'SR_WVre_NDVI',qua=F,trd=F,wide.doy=T,plt=F),silent=T)
	if(class(test.phenol)=='try-error'){next()}

	train.pred<-predn.func(dframe=train.phenol,xvar='SR_WVre_NDVI.doy.50max',yvar='doy.maxN',ndvi='SR_WVre_NDVI')

	int=c(train.pred$coefs[[1]])
	slope=c(train.pred$coefs[[2]])
	
	validate.pred=int+(slope*test.phenol$SR_WVre_NDVI.doy.50max)
	test.phenol$pred.doy.maxN=validate.pred
	
	maxder.lm<-lm(pred.doy.maxN~doy.maxN,data=test.phenol)
	r2=summary(maxder.lm)$adj.r.squared
	slope=coefficients(maxder.lm)[2]
	p=anova(maxder.lm)[[5]][1]

	res[i,]=c(r2,p,slope)
}

res<-res[apply(res,1,sum)!=0,]
apply(res,2,median)