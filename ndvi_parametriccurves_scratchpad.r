####################3Experimental parametric curves section.
#fitting parametric curves to prcntN and to ndvi vals.
#using mmSAR for ndvi curves
install.packages("C:/Users/Vijay/Downloads/mmSAR_1.0.zip", repos = NULL,type='win.binary')
install.packages(c('numDeriv','nortest'))
library(mmSAR)


data(power)
data(expo)
data(negexpo)
data(monod)
data(ratio)
data(logist)
data(lomolino)
data(weibull)
model.list=list(power,expo,negexpo,ratio,logist,lomolino,weibull)
res.list=list()

ndvi.mmsar<-list(name="sr_ndvi",data=data.frame(a=ndvi.df$DOY,s=ndvi.df$SR_WVnir_NDVI))
ndvi.mmsar$data=na.omit(ndvi.mmsar$data)

for(i in 1:length(model.list))
  res.list[[i]]=rssoptim(model.list[[i]],ndvi.mmsar)


mods <- c("power","expo","negexpo","ratio","logist","lomolino","weibull")
rescompare<-data.frame(model=NULL,AICc=NULL,R2=NULL)
for(i in 1:length(mods))
{
  temp=res.list[[i]]
  temp.df=data.frame(model=mods[i],AICc=temp$AICc,R2=temp$R2)
  rescompare=rbind(rescompare,temp.df)
}

Weibull.form=formula(ndvi~c*(1-exp(-z*DOY^f)))
ndvi.df$ndvi=ndvi.df$SR_WVnir_NDVI

ndvi.nls1<-nls(Weibull.form,dat=ndvi.df,start=list(c=.6,z=0.00000000001,f=1))

#fitting all the models to the Galapagos dataset and perform multimodel averaging
resAverage <- multiSAR(modelList=mods,ndvi.mmsar)

logis.starts<-getInitial(ndvi~SSlogis(DOY),data=ndvi.df)
Asym=logis.starts[1]
xmid=logis.starts[2]
scale=logis.starts[3]
ndvi.logis<-nls(ndvi~SSlogis(DOY,Asym,xmid,scale),data=ndvi.df)

coefs<-coef(ndvi.logis)
Asym=coefs[1]
xmid=coefs[2]
scale=coefs[3]


power.func<-function(DOY,coefs)
{
  c=coefs[1]
  z=coefs[2]
  ndvi=c*DOY^z
  return(ndvi)
}
library(minipack.lm)
power.func.exp=as.formula(ndvi~c*DOY^z)
power.nls1=nlsLM(power.func.exp,ndvi.df,start=list(c=.5,z=.5),control=nls.lm.control(maxiter=150))
power.pred=power.func(ndvi.df$DOY,coef(power.nls1))


expon.func=function(dframe,doy)
{
  exp.f=as.formula(ndvi~c+z*log(DOY))
  exp.nls1=nlsLM(exp.f,dframe,start=list(c=.5,z=.5),control=nls.lm.control(maxiter=150))
  print(exp.nls1)
  coefs=coef(exp.nls1)
  c=coefs[1]
  z=coefs[2]
  ndvi=c+z*log(dframe[[doy]])  
  
  plot(dframe[[doy]],dframe$ndvi)
  lines(dframe[[doy]],ndvi,col='blue')
  return(ndvi)
}
exp.pred=expon.func(ndvi.df,'DOY')

negexp.func=function(doy,d,z)
{
  ndvi=d*1-exp(-z*doy)
}

logis.func=function(dframe,doy)
{
  logis=as.formula(ndvi~d/1+exp(-z*DOY+f))
  logis.nls=nlsLM(logis,dframe,start=list(d=.5,z=0.00001,f=9),control=nls.lm.control(maxiter=150))
  print(logis.nls)
  coefs=coef(logis.nls)
  d=coefs[1]
  z=coefs[2]
  f=coefs[3]
  
  ndvi=d/1+exp(-z*dframe$DOY+f)
  dframe=dframe[order(dframe$DOY),]
  plot(dframe$DOY,dframe$ndvi)
  lines(dframe$DOY,ndvi,col='blue')
  return(ndvi)
}

logis.func(ndvi.df,'DOY')

weibull=function(doy,d,z,f)
{
  S=d*(1-exp(-z*doy^f))
}

ndvi.pred=logis.func(ndvi.df$DOY,Asym,xmid,scale)
