l.emo<-ndvi.phenolvals.func(l.ndvi.df,'Yr','NDVI_eMO',qua=F,trd=F,wide.doy=T)

# #coefficients for lonelycamp
#first tried with SR_eMO_NDVI
# maxDer.coefs=coefficients(sr.phenol.maxDer$model)
# halfmax.coefs=coefficients(sr.phenol.50max$model)
# maxVal.coefs=coefficients(sr.phenol.maxval$model)
# firstday.coefs=coefficients(sr.poly.firstday$model)

# #now using SR_rededge
maxDer.coefs=coefficients(sr.re.phenol.maxDer$model)
halfmax.coefs=coefficients(sr.re.phenol.50max$model)
maxVal.coefs=coefficients(sr.re.phenol.maxval$model)
firstday.coefs=coefficients(sr.re.poly.firstday$model)

#now trying it with eMODIS vals
# maxDer.coefs=coefficients(emo.phenol.maxDer$model)
# halfmax.coefs=coefficients(emo.phenol.50max$model)
# maxVal.coefs=coefficients(emo.phenol.maxval$model)
# firstday.coefs=coefficients(emo.poly.firstday$model)

coefs<-rbind(maxDer.coefs,halfmax.coefs,maxVal.coefs,firstday.coefs)
colnames(coefs)=c('Int','Slope')
row.names(coefs)=c('maxDer','halfmax','maxval','firstday')

colville.pred<-function(l.emo,coefs)
{
  maxder.coefs=coefs[1,] #coefficients copied from colville models
  halfmax.coefs=coefs[2,]
  maxval.coefs=coefs[3,]
  
  pred.maxder=maxder.coefs[1]+maxder.coefs[2]*l.emo$NDVI_eMO.doy.maxDeriv
  pred.50max=halfmax.coefs[1]+halfmax.coefs[2]*l.emo$NDVI_eMO.doy.50max
  pred.maxval=maxval.coefs[1]+maxval.coefs[2]*l.emo$NDVI_eMO.doy.maxVal
  
  l.emo<-cbind(l.emo,pred.maxder,pred.50max,pred.maxval)
  return(l.emo)
}

l.emo=colville.pred(l.emo,coefs)

l.emo=l.emo[l.emo$doy.maxN>145,]
#l.emo=l.emo[l.emo$doy.maxN<=188,]


maxder.lm<-lm(doy.maxN~pred.maxder,data=l.emo)
halfmax.lm<-lm(doy.maxN~pred.50max,data=l.emo)
maxval.lm<-lm(doy.maxN~pred.maxval,data=l.emo)


#plotting results
par(mfrow=c(3,1))
par(oma = c(5,4,0,0) + 0.1,
    mar = c(2,0,2,1.5) + 0.1)
plot(doy.maxN~pred.maxder,data=l.emo,xlim=c(160,190),ylim=c(160,200),main='colville maxDer pred')
abline(a=0,b=1,lty=2)
lm_eqn(maxder.lm)

plot(doy.maxN~pred.50max,data=l.emo,xlim=c(160,190),ylim=c(160,200),main='colville 50max pred')
abline(a=0,b=1,lty=2)
lm_eqn(halfmax.lm)

plot(doy.maxN~pred.maxval,data=l.emo,xlim=c(160,190),ylim=c(160,200),main='colville doy.maxval pred')
abline(a=0,b=1,lty=2)
lm_eqn(maxval.lm)

title(xlab = "Obs DOY %N",
      ylab = "Pred DOY %N",
      outer = TRUE, line = 2)