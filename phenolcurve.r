#gam model of ndvi- trying to extract first derivative
#should use a pared down data.frame with only the vars of interest and no NA's. will remove NA's if necessary.

#parameteric curves for nitrogen- gaussian, weibull, lognormal.

phenol_curve_pts<-function(dframe,varname,add=FALSE,plt=TRUE)
{
  dframe=subset(dframe,select=c(varname,'PrcntN','DOY','Yr'),plt=plt)
	library(gam)
	if(sum(is.na(dframe))>0){dframe<-na.omit(dframe)}
	
	if('DOY' %in% names(dframe)==FALSE)
	{
		cat("dataframe must have a numeric DOY column")
		break()
	}
	dframe=dframe[order(dframe$DOY),]
	gam.formula=as.formula(paste(varname,"~s(DOY)"))
	gam1<-gam(gam.formula,data=dframe)
	
	#doy range for prediction
	doy.range=range(dframe$DOY)
	doy.seq=seq(doy.range[1],doy.range[2],by=0.1)
	doy.df=data.frame(DOY=doy.seq)
	
	gam.pred=predict(gam1,doy.df,type='response')
	eps=1e-7 #small difference for comparison
	diff.df=doy.df+eps
	gam.diff=predict(gam1,diff.df,type='response')
	der1=(gam.diff-gam.pred)/eps
	
	####calculating day of year for maximum, max deriv1, and 50%max
	halfmax=max(gam.pred)/2
	doy.maxVal=round(doy.seq[gam.pred==max(gam.pred)])
	doy.posSlope=doy.seq[doy.seq<=doy.maxVal] #only initial increase
	der1.posSlope=der1[doy.seq<=doy.maxVal]
  pred.posSlope=gam.pred[doy.seq<=doy.maxVal]
  
	doy.maxDeriv=min(round(doy.posSlope[der1.posSlope==max(der1.posSlope)]))
  doy.50max=min(round(doy.posSlope[abs(pred.posSlope-halfmax)==min(abs(pred.posSlope-halfmax))]))
  
	
	#plotoutput
  if(plt)
  {
    plot.formula=as.formula(paste(varname,"~DOY"))
    if(!add)
    {
      plot(plot.formula,data=dframe)
    }
    lines(gam.pred~doy.seq,col='blue',lwd=2)
    points(doy.maxVal,gam.pred[doy.seq==doy.maxVal],lwd=5,col='red')
    points(doy.maxDeriv,gam.pred[doy.seq==doy.maxDeriv],lwd=5,col='green')
    points(doy.50max,gam.pred[doy.seq==doy.50max],lwd=5,col='orange')
  }
	
	
	pheno.df=data.frame(doy.maxVal=doy.maxVal,doy.50max=doy.50max,doy.maxDeriv=doy.maxDeriv)
	return(pheno.df)
}
	
ndvi.phenolvals.func<-function(dframe,yearname,ndvi,add=F,plt=TRUE)
{
  yrs=sort(unique(dframe[[yearname]]))
  ndvi.phenolvals=data.frame(Year=NULL,doy.MaxVal=NULL,doy.50Max=NULL,doy.maxDer1=NULL,doy.maxN=NULL)
  for(i in 1:length(yrs))
  {
    yr=data.frame(Year=yrs[i])
    temp=dframe[dframe[[yearname]]==yrs[i],]
    ndvi.phenol<-phenol_curve_pts(temp,ndvi,add=add,plt=plt)
    doy.maxN<-phenol_curve_pts(temp,'PrcntN')$doy.maxVal
    ndvi.phenol.yr=cbind(yr,ndvi.phenol,doy.maxN)
    ndvi.phenolvals=rbind(ndvi.phenolvals,ndvi.phenol.yr)
  }
  ndvi.forcor=ndvi.phenolvals[,2:4]
  print(apply(ndvi.forcor,2,function(x) cor.test(x,ndvi.phenolvals$doy.maxN)))
  
  return(ndvi.phenolvals)
}


ndvi.phenolpoly.func<-function(dframe,yearname,ndvi,porder,add=FALSE,plt=TRUE)
{
  yrs=sort(unique(dframe[[yearname]]))
  ndvi.phenolvals=data.frame(Year=NULL,doy.MaxVal=NULL,doy.50Max=NULL,doy.maxDer1=NULL,doy.maxN=NULL)
  for(i in 1:length(yrs))
  {
    yr=data.frame(Year=yrs[i])
    temp=dframe[dframe[[yearname]]==yrs[i],]
    ndvi.phenol<-phenol_poly_pts(temp,ndvi,porder,add=add,plt=plt)
    doy.maxN<-phenol_poly_pts(temp,'PrcntN',porder,plt=FALSE)$doy.maxVal
    ndvi.phenol.yr=cbind(yr,ndvi.phenol,doy.maxN)
    ndvi.phenolvals=rbind(ndvi.phenolvals,ndvi.phenol.yr)
  }
  ndvi.forcor=ndvi.phenolvals[,2:4]
  print(apply(ndvi.forcor,2,function(x) cor.test(x,ndvi.phenolvals$doy.maxN)))
  
  return(ndvi.phenolvals)
}

phenol_poly_pts<-function(dframe,varname,porder,add=FALSE,plt=TRUE)
{
  #need to modify to allow flexible specification of different polynomials for ndvi vs. N.
  dframe=subset(dframe,select=c(varname,'PrcntN','DOY','Yr'))

  if(sum(is.na(dframe))>0){dframe<-na.omit(dframe)}
  
  if('DOY' %in% names(dframe)==FALSE)
  {
    cat("dataframe must have a numeric DOY column")
    break()
  }
  
  dframe=dframe[order(dframe$DOY),]
  order.seq=1:porder
  
  fo.poly=paste("I(DOY^",order.seq,")",sep='')
  fo.poly=paste(fo.poly,collapse='+')
  poly.formula=as.formula(paste(paste(varname,"~",fo.poly)))
  poly.lm<-lm(poly.formula,data=dframe)
  
  #doy range for prediction
  doy.range=range(dframe$DOY)
  doy.seq=seq(doy.range[1],doy.range[2],by=0.1)
  doy.df=data.frame(DOY=doy.seq)
  
  poly.pred=predict(poly.lm,doy.df,type='response')
  eps=1e-7 #small difference for comparison
  diff.df=doy.df+eps
  poly.diff=predict(poly.lm,diff.df,type='response')
  der1=(poly.diff-poly.pred)/eps
  
  ####calculating day of year for maximum, max deriv1, and 50%max
  halfmax=max(poly.pred)/2
  doy.maxVal=round(doy.seq[poly.pred==max(poly.pred)])
  doy.posSlope=doy.seq[doy.seq<=doy.maxVal] #only initial increase
  der1.posSlope=der1[doy.seq<=doy.maxVal]
  pred.posSlope=poly.pred[doy.seq<=doy.maxVal]
  
  doy.maxDeriv=min(round(doy.posSlope[der1.posSlope==max(der1.posSlope)]))
  doy.50max=min(round(doy.posSlope[abs(pred.posSlope-halfmax)==min(abs(pred.posSlope-halfmax))]))
  
  
  #plotoutput
  plot.formula=poly.formula
  if(plt)
  {
    if(!add)
    {
      plot(as.formula(paste(varname,"~DOY")),data=dframe)
    }
    lines(poly.pred~doy.seq,col='blue',lty=2,lwd=3)
    points(doy.maxVal,poly.pred[doy.seq==doy.maxVal],lwd=5,col='red',pch=2)
    points(doy.maxDeriv,poly.pred[doy.seq==doy.maxDeriv],lwd=5,col='green',pch=2)
    points(doy.50max,poly.pred[doy.seq==doy.50max],lwd=5,col='orange',pch=2)
    
  }
      
  pheno.df=data.frame(doy.maxVal=doy.maxVal,doy.50max=doy.50max,doy.maxDeriv=doy.maxDeriv)
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


lognormal=function(x, mu, sd, C) {							#lognormal
  ans <- vector(length = length(x), mode = "numeric")
  for (i in 1:length(x)) {
    if (x[i]>0) {
      value <- exp(-(log10(x[i])-mu)^2/(2*sd^2))/(sqrt(2*pi)*sd*x[i])*C
      ans[i] <- value
    }  else {  ans[i] <- 0  }  };  return(ans)  }


model.compare=function(yvar,dframe)
{
  library(gam)
  models=list()
  models[[1]]=doy.lm(yvar,dframe)
  models[[2]]=quadratic(yvar,dframe)
  models[[3]]=thirdorder(yvar,dframe)
  models[[4]]=doy.gam(yvar,dframe)
  
  AIC.table=sapply(models,AIC)
  AIC.table=data.frame(model=c('linear','quadratic','thirdorder','gam'),AIC=AIC.table)
  AIC.table=AIC.table[order(AIC.table$AIC),]
  print(models)
  res.list=list(models,AIC.table)
  return(res.list)
}
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
