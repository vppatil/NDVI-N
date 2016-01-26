#should probably actually be bootstrapping the prediction curves to get an idea of the range, median value/diff.


vegplots<-read.csv('VegPlots_for R.csv')
ndvi.df<-subset(vegplots,select=c(ndvi.names,'DOY','Yr','PrcntN'))
source('../Code/phenolcurve.r')

ndvi.names=c('SR_WVnir_NDVI','eMO_NDVI','SR_eMO_NDVI','SR_WVre_NDVI')
yearname='Yr'

b=10000
ind=1:dim(ndvi.df)[1]
boot.out=matrix(0,b,5)

for(i in 1:b)
{
	boot.ind=sample(ind,length(ind),replace=T)
	df.boot=ndvi.df[boot.ind,]

	emo.poly.allyr=phenol_poly_pts(df.boot,ndvi.names[2],3,plt=FALSE)
	emo.poly=ndvi.phenolpoly.func(df.boot,yearname,ndvi.names[2],3,plt=FALSE)

	emo.poly<-pred.obs.compare(emo.poly,'doy.maxDeriv','doy.maxN','maxDer',plt=FALSE)
	emo.poly.maxDer=predn.func(emo.poly,'doy.maxDeriv','doy.maxN',plt=FALSE)

	r2=emo.poly.maxDer$rsq
	p=emo.poly.maxDer$anova.t[5][1,1]
	bias=mean(emo.poly.maxDer$model$resid)
	coefs=sapply(emo.poly.maxDer$coefs,c)
	boot.out.line=c(r2,p,bias,coefs)
	boot.out[i,]=boot.out.line
}

colnames(boot.out)=c('r2','p','bias','int','slope')
boot.out<-data.frame(boot.out)
