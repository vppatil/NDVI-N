#gam model of ndvi- trying to extract first derivative
#should use a pared down data.frame with only the vars of interest and no NA's. will remove NA's if necessary.

#parameteric curves for nitrogen- gaussian, weibull, lognormal.

#need to fix phenolcurve so maxval does not end up at the minimum instead.
#ndvi.phenolvals.func constrained to fit 3rd order polynomial for prcnt N.
phenol_models_aictab<-function(dframe,varname,use.aicc=FALSE)
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
    
    model.list=list(quad,third,gaus,logn)
    
	##to make sure you are getting the same model(s)
  
	  aic.tab=AIC(quad,third,gaus,logn)#model selection table
	  aic.tab$n=sapply(model.list,nobs)
	  aic.tab=aic.table(aic.tab,'AIC',k='df',n='n',aicc=use.aicc)
 
	  return(aic.tab)
}



ndvi.phenolvals.func<-function(dframe,yearname,ndvi,modav=TRUE,qua=TRUE,trd=TRUE,gau=TRUE,lgn=TRUE,wide.doy=FALSE,plt=TRUE)
{
  yrs=sort(unique(dframe[[yearname]]))
  ndvi.phenolvals=data.frame(Year=NULL,doy.MaxVal=NULL,doy.50Max=NULL,doy.maxDer1=NULL,doy.maxN=NULL)
  for(i in 1:length(yrs))
  {
    yr=data.frame(Year=yrs[i])
    temp=dframe[dframe[[yearname]]==yrs[i],]
    ndvi.phenol<-phenol_curve_pts(temp,ndvi,modav=modav,qua=qua,trd=trd,gau=gau,lgn=lgn,wide.doy=wide.doy)
    doy.maxN<-phenol_poly_pts(temp,'PrcntN',porder=3,wide.doy=wide.doy)[,1]
    ndvi.phenol.yr=cbind(yr,ndvi.phenol,doy.maxN)
    ndvi.phenolvals=rbind(ndvi.phenolvals,ndvi.phenol.yr)
  }
  ndvi.forcor=ndvi.phenolvals[,2:4]
  print(apply(ndvi.forcor,2,function(x) cor.test(x,ndvi.phenolvals$doy.maxN)))
  
  return(ndvi.phenolvals)
}

phenol_curve_pts<-function(dframe,varname,modav=TRUE,qua=TRUE,trd=TRUE,gau=TRUE,lgn=TRUE,wide.doy=FALSE,plt=TRUE)
{
    modav=modav
    use.vec=c(qua,trd,gau,lgn)

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
    if(wide.doy){doy.range[1]=140} #extending predictions to mid may

    doy.seq=seq(doy.range[1],doy.range[2],by=0.1)
    doy.df=data.frame(DOY=doy.seq)
    
    model.list=list(quad,third,gaus,logn)[use.vec]
    
    for(i in 1:length(model.list))
	{print(summary(model.list[[i]]))}    ##to make sure you are getting the same model(s)
    
     

	if(modav) #either use model-averaged values or pick best model.
	{
	  aic.tab=AIC(quad,third,gaus,logn)[use.vec,]#model selection table
	  aic.tab$n=sapply(model.list,nobs)
	  aic.tab=aic.table(aic.tab,'AIC',k='df',n='n',aicc=F)
	  predvals=model.avg.prediction(model.list,aic.tab$w,doy.df)
   
	}else
	{
		aic.tab<-data.frame(AIC=sapply(model.list,AIC),w=1)
		rows=1:length(model.list)
		best=rows[aic.tab$w==max(aic.tab$w)]
		predvals=predict(model.list[[best]],doy.df)
	}

	
	doy.df[[varname]]=predvals

    
    ####calculating day of year for maximum, max deriv1, and 50%max
    #finite diff approximation for first derivative of curve
    eps=1e-7 #small difference for comparison
    diff.df=data.frame(DOY=doy.df$DOY+eps) 
	
	if(modav)
	{
		diff.pred=model.avg.prediction(model.list,aic.tab$w,diff.df)
	}else
	{
		rows=1:length(model.list)
		best=rows[aic.tab$w==max(aic.tab$w)]
		diff.pred=predict(model.list[[best]],diff.df)		
	}
   
	der1=(diff.pred-predvals)/eps #approximate instantaneous dy/dx


 ##########################################################################################
 ####calculating day of year for maximum, max deriv1, and 50%max old style
  halfmax=max(predvals)/2 ##3
  doy.maxVal=round(doy.seq[predvals==max(predvals)])
  doy.posSlope=doy.seq[doy.seq<=doy.maxVal] #only initial increase
  der1.posSlope=der1[doy.seq<=doy.maxVal]
  pred.posSlope=predvals[doy.seq<=doy.maxVal]

  doy.maxDeriv=min(round(doy.posSlope[der1.posSlope==max(der1.posSlope)]))
  doy.50max=min(round(doy.posSlope[abs(pred.posSlope-halfmax)==min(abs(pred.posSlope-halfmax))]))
  
  first.obs.day=min(doy.seq)

  ########################################################################################	
  #cumbersome approach to identifying peak based on general location of curve apex
  #should avoid the problem of simply picking the right hand side of a thirdorder curve
  #pos.slope.yesterday=vector('logical',length=length(doy.seq))
  #neg.slope.tomorrow=vector('logical',length=length(doy.seq))
  #for(i in 10:(length(doy.seq)-10))
  #{
  #  pos.slope.yesterday[i]=der1[i-9]>0
  #  neg.slope.tomorrow[i]=der1[i+9]<0
  #}
  #
  #pos.slope.yesterday[1:10]=TRUE
  #neg.slope.tomorrow[(length(doy.seq)-10):length(doy.seq)]=neg.slope.tomorrow[length(doy.seq)-11]
  # 
  #peak.seq=doy.seq[pos.slope.yesterday & neg.slope.tomorrow]
  #peak.pred=predvals[doy.seq %in% peak.seq]
  #maxval=max(peak.pred)
  #doy.maxVal=round(peak.seq[peak.pred==maxval]) #date of maximum value #now rounded
  ##doy.maxVal=peak.seq[peak.pred==maxval] #date of maximum value #now rounded
  	
  #halfmax=maxval/2
  #doy.posSlope=doy.seq[doy.seq<=doy.maxVal] #only initial increase
    
  #der1.posSlope=der1[doy.seq<=doy.maxVal]      #truncating the first derivative and curve predictions 
  #pred.posSlope=predvals[doy.seq<=doy.maxVal]  #to initial uphill slope of curve
    
  ##finding the point of maximum rate of change
  #doy.maxDeriv=round(doy.posSlope[der1.posSlope==max(der1.posSlope)])
  ##doy.maxDeriv=doy.posSlope[der1.posSlope==max(der1.posSlope)]

  ##finding day with predicted %N as close as possible to half max value
  #doy.50max=round(doy.posSlope[abs(pred.posSlope-halfmax)==min(abs(pred.posSlope-halfmax))])
  ##doy.50max=doy.posSlope[abs(pred.posSlope-halfmax)==min(abs(pred.posSlope-halfmax))]
  #######################################################################################
  pheno.df=data.frame(doy.maxVal=doy.maxVal,doy.50max=doy.50max,doy.maxDeriv=doy.maxDeriv,first.obs.day=first.obs.day)
  
  if(plt)
  {
	plot(dframe$DOY,dframe[[varname]],ylab=varname,xlab='DOY',xlim=doy.range)
	lines(doy.df$DOY,doy.df[[varname]])
	points(doy.maxVal,predvals[doy.seq==doy.maxVal],col='red',pch=15,lwd=3,cex=2)
	points(doy.50max,predvals[doy.seq==doy.50max],col='green',pch=16,lwd=3,cex=2)
	points(doy.maxDeriv,predvals[doy.seq==doy.maxDeriv],col='blue',pch=17,lwd=3,cex=2)
  }
  
    names(pheno.df)=paste(varname,names(pheno.df),sep='.')
    return(pheno.df)
}

#polynomial curve fits and comparison functions for regressing variables against DOY
quadratic<-function(yvar,dframe)
{
  dframe=subset(dframe,select=c('Yr','DOY',yvar))
  if(sum(is.na(dframe))>0)
  {
    dframe<-na.omit(dframe)
  }
  f=as.formula(paste(yvar,"~I(DOY^2)+DOY",sep=''))
  mod=glm(f,data=dframe)
  res=list(AIC=AIC(mod),model=summary(mod))
  return(mod)
}

thirdorder<-function(yvar,dframe)
{
  dframe=subset(dframe,select=c('Yr','DOY',yvar))
  if(sum(is.na(dframe))>0)
  {
    dframe<-na.omit(dframe)
  }
  f=as.formula(paste(yvar,"~I(DOY^3)+I(DOY^2)+DOY",sep=''))
  mod=glm(f,data=dframe)
  res=list(AIC=AIC(mod),model=summary(mod))
  return(mod)
}

doy.lm<-function(yvar,dframe)
{
  dframe=subset(dframe,select=c('Yr','DOY',yvar))
  if(sum(is.na(dframe))>0)
  {
    dframe<-na.omit(dframe)
  }
  f=as.formula(paste(yvar,"~DOY",sep=''))
  print(f)
  mod=glm(f,data=dframe)
  return(mod)
}

doy.gam<-function(yvar,dframe)
{
  dframe=subset(dframe,select=c('Yr','DOY',yvar))
  if(sum(is.na(dframe))>0)
  {
    dframe<-na.omit(dframe)
  }
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

predn.func<-function(dframe,xvar,yvar,plt=TRUE,xaxt=NULL,yaxt=NULL,ndvi='')
{
  pred.lm<-lm(as.formula(paste(yvar,"~",xvar)),data=dframe)
  yvar.pred=predict(pred.lm,type='response')
  print(yvar.pred)
  int=coef(pred.lm)[1]
  slope=coef(pred.lm)[2]
  
  if(plt)
  {
    plot(dframe[[yvar]],yvar.pred,xlim=c(160,190),ylim=c(160,190),xlab='%N obs',ylab='%N pred',xaxt=xaxt,yaxt=yaxt,main=xvar)
    abline(b=1,a=0,lty=2)
  }

  res=list(model=pred.lm,coefs=data.frame(int=int,slope=slope),predn=yvar.pred,rsq=summary(pred.lm)$adj.r.squared,anova.t=anova(pred.lm))
  return(res)
}

lm_eqn <- function(m){
    a=round(coef(m)[1],2)
    b=round(coef(m)[2],2)
    r2=round(summary(m)$adj.r.squared,2)
    p=round(anova(m)[["Pr(>F)"]][1],3)
    #eq <- expression(paste0("y=",a,"+",b,"*x/n", R^{2},'=',r2,",p=",p))
    
    mylabel1 = bquote(italic(R)^2 == .(format(r2, digits = 3)))
    mylabel2 = bquote(italic(Y) == .(format(a,digits=3))+.(format(b,digits=3))*italic(X))
    mylabel3 = bquote(italic(p)==.(format(p,digits=3)))
    text(182,172,mylabel2)
    text(182,168,mylabel1)
    text(182,164,mylabel3)
}

pred.obs.compare=function(pheno.df,xvar,yvar,pheno.met,plt=FALSE)
{
  temp=predn.func(pheno.df,xvar,yvar,plt=plt)
  pheno.df$predn=temp$predn
  names(pheno.df)[names(pheno.df)=='predn']=paste('predn',pheno.met,sep='_')
  pheno.df$resid=temp$model$residuals
  names(pheno.df)[names(pheno.df)=='resid']=paste('resid',pheno.met,sep='_')
  return(pheno.df)
}

###old polynomial fit code to reproduce those great initial fits.
ndvi.phenolpoly.func<-function(dframe,yearname,ndvi,porder,add=FALSE,plt=TRUE,wide.doy=FALSE)
{
  yrs=sort(unique(dframe[[yearname]]))
  ndvi.phenolvals=data.frame(Year=NULL,doy.MaxVal=NULL,doy.50Max=NULL,doy.maxDer1=NULL,doy.maxN=NULL)
  for(i in 1:length(yrs))
  {
    yr=data.frame(Year=yrs[i])
    temp=dframe[dframe[[yearname]]==yrs[i],]
    ndvi.phenol<-phenol_poly_pts(temp,ndvi,porder,add=add,plt=plt,wide.doy=wide.doy)
    doy.maxN<-phenol_poly_pts(temp,'PrcntN',porder,plt=FALSE,wide.doy)$doy.maxVal
    ndvi.phenol.yr=cbind(yr,ndvi.phenol,doy.maxN)
    ndvi.phenolvals=rbind(ndvi.phenolvals,ndvi.phenol.yr)
  }
  ndvi.forcor=ndvi.phenolvals[,2:4]
  print(apply(ndvi.forcor,2,function(x) cor.test(x,ndvi.phenolvals$doy.maxN)))
  
  return(ndvi.phenolvals)
}

phenol_poly_pts<-function(dframe,varname,porder,add=FALSE,plt=TRUE,wide.doy=FALSE)
{
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
  if(wide.doy){doy.range[1]=140} #extending predictions to mid may
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
  
    first.obs.day=min(doy.seq)

  #plotoutput
  plot.formula=poly.formula
  if(plt)
  {
    if(!add)
    {
      plot(as.formula(paste(varname,"~DOY")),data=dframe)
    }
    lines(poly.pred~doy.seq)
    points(doy.maxVal,poly.pred[doy.seq==doy.maxVal],col='red',lwd=3,cex=2,pch=15)
    points(doy.50max,poly.pred[doy.seq==doy.50max],col='blue',lwd=3,cex=2,pch=17)
    points(doy.maxDeriv,poly.pred[doy.seq==doy.maxDeriv],col='green',lwd=3,cex=2,pch=16)
    
  }
      
  pheno.df=data.frame(doy.maxVal=doy.maxVal,doy.50max=doy.50max,doy.maxDeriv=doy.maxDeriv,first.obs.day=first.obs.day)
  return(pheno.df)
}

