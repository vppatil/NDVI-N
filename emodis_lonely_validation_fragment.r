#lonely.emodis validation
colville.pred2<-function(l.emo,coefs)
{
  maxder.coefs=coefs[1,] #coefficients copied from colville models
  halfmax.coefs=coefs[2,]
  maxval.coefs=coefs[3,]
  
  l.emo$pred.maxder=maxder.coefs[1]+maxder.coefs[2]*l.emo$eMO_NDVI.doy.maxDeriv
  l.emo$pred.50max=halfmax.coefs[1]+halfmax.coefs[2]*l.emo$eMO_NDVI.doy.50max
  l.emo$pred.maxval=maxval.coefs[1]+maxval.coefs[2]*l.emo$eMO_NDVI.doy.maxVal
  
  return(l.emo)
}
lonely.emo<-read.csv('../Data/lonely_emo_transformed.csv')
lonely.emo$Yr=lonely.emo$Year
lonely.emo<-subset(lonely.emo,lonely.emo$Single.Pixel==0,select=c('Year','Yr','DOY','eMO_NDVI','Plot'))
phenol_curve_pts(lonely.emo,'eMO_NDVI',qua=F,trd=F)

#using 2012 only- for predicting 2012 median date of peak N at lonely.
lonely.emo<-subset(lonely.emo,lonely.emo$Year==2012)

#fitting phenol curves for each pixel
l.emo.p1<-phenol_curve_pts(lonely.emo[lonely.emo$Plot=='Pixel1',],'eMO_NDVI',qua=F,trd=F)
l.emo.p2<-phenol_curve_pts(lonely.emo[lonely.emo$Plot=='Pixel2',],'eMO_NDVI',qua=F,trd=F)
l.emo.p3<-phenol_curve_pts(lonely.emo[lonely.emo$Plot=='Pixel3',],'eMO_NDVI',qua=F,trd=F)

l.emo<-rbind(l.emo.p1,l.emo.p2,l.emo.p3)
l.emo=colville.pred2(l.emo,emo.coefs)