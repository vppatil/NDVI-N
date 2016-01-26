plot(SR_WVnir_NDVI~DOY,data=ndvi.df)
#############################################################################
#gaussian process regression
gaus.f=function(yvar){return(formula(paste(yvar,"~ k*exp(-1/2*(DOY-mu)^2/sigma^2)")))}
emo.gaus=gaus.f('SR_WVnir_NDVI') #use mean x, sd x, and y height as starting pars

gaus.emo.nls<- nlsLM(emo.gaus, data = ndvi.df,start=c(mu=200,sigma=20,k=1))
gaus.pred=gaus.emo.nls$m$predict(pred.df)
lines(pred.df$DOY,gaus.pred,type='l')

##############################################################################
#3 parameter weibull from stoneman et al. 1996 with intercept added
#bottom is 3 par formulation from maltamo et al.

weibull.func=function(xvar,a,b,c)
{
	ans<-vector(length=length(xvar),mode='numeric')
	for (i in 1:length(xvar)) 
	{

			#value<- a*(c/b)*((xvar[i]/b)^(c^-1))*(exp(-(xvar[i]/b)^c))
			value=(c/b)*(((xvar[i]-a)/b)^(c-1))*exp(-((xvar[i]-a)/b)^c)
			ans[i] <- value
	}  
	return(ans)  
}

xvar=pred.df$DOY
a=10;b=140;c=20
wei.pred=(c/b)*(((xvar-a)/b)^(c-1))*exp(-((xvar-a)/b)^c)

test<-nls(SR_WVnir_NDVI~weibull.func(DOY,a,b,c),data=ndvi.df,start=c(a=160,b=54,c=16),control=nls.lm.control(maxiter=100))

wei.pred=test$m$predict(pred.df)
plot(xvar,wei.pred)

lines(pred.df$DOY,wei.pred,col='blue')
589.702 228.690   2.078
xvar=160:220

for(i in seq(10,200,by=10))
{
	Sys.sleep(.3)
	for(j in seq(10,200,by=10))
	{
		Sys.sleep(.3)
		plot(weibull.func(xvar,120,i,j),main=paste(c("b=",i,",b=",j),collapse=""))
	}
}

weibull.form<-function(yvar)
{
	rhs="~(c/b)*(((x-a)/b)^(c-1))*exp(-((x-a)/b)^c))"
	f=paste(yvar,rhs,collapse='')
	return(f)
}

lognormal=function(x, mu, sd, C) {
  ans <- vector(length = length(x), mode = "numeric")
  for (i in 1:length(x)) {
    if (x[i]>0) {
      value <- exp(-(log10(x[i])-mu)^2/(2*sd^2))/(sqrt(2*pi)*sd*x[i])*C
      ans[i] <- value
    }  else {  ans[i] <- 0  }  };  return(ans)  }

wvnir.nls=nlsLM(PrcntN~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df)
logn.pred=wvnir.nls$m$predict(pred.df)
lines(pred.df$DOY,logn.pred,col='blue')

wvnir.nls=nlsLM(SR_WVnir_NDVI~lognormal(DOY,mu,sd,C),start=c(mu=log(200),sd=log(37),C=1),data=ndvi.df)
logn.pred=wvnir.nls$m$predict(pred.df)
plot(SR_WVnir_NDVI~DOY,data=ndvi.df,xlim=c(0,400),ylim=c(0,1))
lines(pred.df$DOY,logn.pred,col='blue')
