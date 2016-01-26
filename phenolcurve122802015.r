#gam model of ndvi- trying to extract first derivative
#should use a pared down data.frame with only the vars of interest and no NA's. will remove NA's if necessary.

#parameteric curves for nitrogen- gaussian, weibull, lognormal.

phenol_curve_pts<-function(dframe,varname,modavg=TRUE)
{
    dframe=subset(dframe,select=c(varname,'DOY','Yr'))
    
    if(sum(is.na(dframe))>0){dframe<-na.omit(dframe)}
    
    if('DOY' %in% names(dframe)==FALSE)
    {
        cat("dataframe must have a numeric DOY column")
        break()
    }
    dframe=dframe[order(dframe$DOY),]
    
    quad=quadratic(dframe=dframe,yvar=varname)
    third=thirdorder(dframe=dframe,yvar=varname)
    gaus=nlsLM(gaus.f(varname),data=dframe,start=c(mu=200,sigma=20,k=1),
               control=nls.lm.control(maxiter=150))
    logn=nlsLM(as.formula(paste(varname,"~lognormal(DOY,mu,sd,C)",collapse='')),
               data=dframe,start=c(mu=log(200),sd=log(37),C=1),
               control=nls.lm.control(maxiter=150))
    
    #doy range for prediction
    doy.range=range(dframe$DOY)
    doy.seq=seq(doy.range[1],doy.range[2],by=0.1)
    doy.df=data.frame(DOY=doy.seq)
    
    model.list=list(quad,third,gaus,logn)
    
    aic.tab=AIC(quad,third,gaus,logn)#model selection table
    aic.tab$n=sapply(model.list,nobs)
    aic.tab=aic.table(aic.tab,'AIC',k='df',n='n',aicc=T)
	aic.tab$w<-ifelse(aic.tab$df==5,1,0) #cheater line to force selection of third order poly every time.
    
	if(modavg) #either use model-averaged values or pick best model.
	{
		predvals=model.avg.prediction(model.list,aic.tab$w,doy.df)
	}else
	{
		rows=1:length(model.list)
		best=rows[aic.tab$w==max(aic.tab$w)]
		predvals=predict(model.list[[best]],doy.df)
	}
	
	doy.df[[varname]]=predvals

    
    ####calculating day of year for maximum, max deriv1, and 50%max
    #finite diff approximation for first derivative of curve
    eps=1e-7 #small difference for comparison
    diff.df=data.frame(DOY=doy.df$DOY+eps) 
	
	if(modavg)
	{
		diff.pred=model.avg.prediction(model.list,aic.tab$w,diff.df)
	}else
	{
		rows=1:length(model.list)
		best=rows[aic.tab$w==max(aic.tab$w)]
		diff.pred=predict(model.list[[best]],diff.df)		
	}
   
	der1=(diff.pred-predvals)/eps #approximate instantaneous dy/dx
	
	
	
	doy.maxVal=min(round(doy.seq[abs(der1-0)==min(abs(der1-0))]))
	maxval=predvals[doy.seq==doy.maxVal]
	
	halfmax=maxval/2
    doy.posSlope=doy.seq[doy.seq<=doy.maxVal] #only initial increase
    
    der1.posSlope=der1[doy.seq<=doy.maxVal]      #truncating the first derivative and curve predictions 
    pred.posSlope=predvals[doy.seq<=doy.maxVal]  #to initial uphill slope of curve
    
    #finding the point of maximum rate of change, rounded to nearest whole day
    doy.maxDeriv=min(round(doy.posSlope[der1.posSlope==max(der1.posSlope)]))
    
    #finding whole day with predicted %N as close as possible to half max value
    doy.50max=min(round(doy.posSlope[abs(pred.posSlope-halfmax)==min(abs(pred.posSlope-halfmax))]))
    
    pheno.df=data.frame(doy.maxVal=doy.maxVal,doy.50max=doy.50max,doy.maxDeriv=doy.maxDeriv)

	plot(dframe$DOY,dframe[[varname]],ylab=varname,xlab='DOY')
	lines(doy.df$DOY,doy.df[[varname]])
	points(doy.maxVal,predvals[doy.seq==doy.maxVal],col='red',pch=15,lwd=3,cex=2)
	points(doy.50max,predvals[doy.seq==doy.50max],col='green',pch=16,lwd=3,cex=2)
	points(doy.maxDeriv,predvals[doy.seq==doy.maxDeriv],col='blue',pch=17,lwd=3,cex=2)
	

    names(pheno.df)=paste(varname,names(pheno.df),sep='.')
    return(pheno.df)
}

#polynomial curve fits and comparison functions for regressing variables against DOY
quadratic<-function(yvar,dframe)
{
  f=as.formula(paste(yvar,"~I(DOY^2)+DOY",sep=''))
  mod=glm(f,data=dframe)
  res=list(AIC=AIC(mod),model=summary(mod))
  return(mod)
}

thirdorder<-function(yvar,dframe)
{
  f=as.formula(paste(yvar,"~I(DOY^3)+I(DOY^2)+DOY",sep=''))
  mod=glm(f,data=dframe)
  return(mod)
}

doy.lm<-function(yvar,dframe)
{
  f=as.formula(paste(yvar,"~DOY",sep=''))
  print(f)
  mod=glm(f,data=dframe)
  return(mod)
}

doy.gam<-function(yvar,dframe)
{
  f=as.formula(paste(yvar,"~s(DOY)",sep=''))
  print(f)
  mod=gam(f,data=dframe)
  return(mod)
}

gaus.f=function(yvar){return(formula(paste(yvar,"~ k*exp(-1/2*(DOY-mu)^2/sigma^2)")))} #gaussian FORMULA (not model fit)
#inconsistent with other functions because I got lazy

lognormal=function(x, mu, sd, C) {							#lognormal
  ans <- vector(length = length(x), mode = "numeric")
  for (i in 1:length(x)) {
    if (x[i]>0) {
      value <- exp(-(log10(x[i])-mu)^2/(2*sd^2))/(sqrt(2*pi)*sd*x[i])*C
      ans[i] <- value
    }  else {  ans[i] <- 0  }  };  return(ans)  }

predn.func<-function(dframe,xvar,yvar,plt=TRUE)
{
  pred.lm<-lm(as.formula(paste(yvar,"~",xvar)),data=dframe)
  yvar.pred=predict(pred.lm,type='response')
  print(yvar.pred)
  int=coef(pred.lm)[1]
  slope=coef(pred.lm)[2]
  
  if(plt)
  {
    plot(dframe[[yvar]],yvar.pred)
    abline(b=1,a=0,lty=2)
  }
  
  res=list(model=pred.lm,coefs=data.frame(int=int,slope=slope),predn=yvar.pred,rsq=summary(pred.lm)$adj.r.squared,anova.t=anova(pred.lm))
  return(res)
}

pred.obs.compare=function(pheno.df,xvar,yvar,pheno.met,plt=TRUE)
{
  temp=predn.func(pheno.df,xvar,yvar,plt=plt)
  pheno.df$predn=temp$predn
  names(pheno.df)[names(pheno.df)=='predn']=paste('predn',pheno.met,sep='_')
  pheno.df$resid=temp$model$residuals
  names(pheno.df)[names(pheno.df)=='resid']=paste('resid',pheno.met,sep='_')
  return(pheno.df)
}
