# Colville Veg / NDVI data exploration

# Set working directory
setwd("C:/Data/Dan/Changing Arctic Ecosystems/NDVI summaries Kyle 2015")

# Import working file
ColVeg<-read.csv(file = "VegPlots_for R.csv")

# Verify the structure of the database
str(ColVeg);head(ColVeg)
# How does biomass vary by NDVI? Looks at SR_WVnir_NDVI
# Refine dataset for sake of plotting (see later): remove all NA records
NDVIbiomass<-ColVeg[which(ColVeg$SR_WVnir_NDVI>0&ColVeg$Biomass>0),];str(NDVIbiomass)
# Sort on NDVI: helps figure clarity
NDVIbiomass<-NDVIbiomass[order(NDVIbiomass$SR_WVnir_NDVI),];head(NDVIbiomass)

# Raw data, coded by year
par(mar=c(5.1,6.1,4.1,2.1))
plot(ColVeg$SR_WVnir_NDVI,ColVeg$Biomass,col= ifelse(ColVeg$Yr==2011,"black",ifelse(ColVeg$Yr==2012,"red",ifelse(ColVeg$Yr==2013,"blue",ifelse(ColVeg$Yr==2014,"green","orange")))),cex=2,lwd=2,xlab="NDVI",ylab="Biomass (g)",cex.lab=1.75,cex.axis=1.5) 
text(.3,200,"2011",col="black",cex=2);text(.3,180,"2012",col="red",cex=2);text(.3,160,"2013",col="blue",cex=2);text(.3,140,"2014",col="green",cex=2);text(.3,120,"2015",col="orange",cex=2) # NDVI only exceeded .8 in 2011 essentially...biomass >150 in 2013 only essentially....
# coded by treatment
plot(ColVeg$SR_WVnir_NDVI,ColVeg$Biomass,col= ifelse(ColVeg$Ex_Co=="Ex","black","red"),cex=2,lwd=2,xlab="NDVI",ylab="Biomass (g)",cex.lab=1.75,cex.axis=1.5)
text(.3,200,"Exclosed",col="black",cex=2);text(.3,180,"Control",col="red",cex=2)

# Drop 2011: Kyle and David aren't super confident of spectrometer techniques in this first year...
NDVIbiomass1215<-NDVIbiomass[which(NDVIbiomass$Yr!='2011'),];str(NDVIbiomass1215)
# Drop unused levels
drop.levels <- function(dat){
  # Drop unused factor levels from all factors in a data.frame
  # Author: Kevin Wright.  Idea by Brian Ripley. Found on http://rwiki.sciviews.org/doku.php?id=tips:data-manip:drop_unused_levels
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}
NDVIbiomass1215<-drop.levels(NDVIbiomass1215)
levels(NDVIbiomass1215$Yr)
# Plot again by year, with 2011 removed:
par(mar=c(5.1,6.1,4.1,2.1))
plot(NDVIbiomass1215$SR_WVnir_NDVI,NDVIbiomass1215$Biomass,col= ifelse(NDVIbiomass1215$Yr==2011,"black",ifelse(NDVIbiomass1215$Yr==2012,"red",ifelse(NDVIbiomass1215$Yr==2013,"blue",ifelse(NDVIbiomass1215$Yr==2014,"green","orange")))),cex=2,lwd=2,xlab="NDVI",ylab="Biomass (g)",cex.lab=1.75,cex.axis=1.5) 
text(.3,180,"2012",col="red",cex=2);text(.3,160,"2013",col="blue",cex=2);text(.3,140,"2014",col="green",cex=2);text(.3,120,"2015",col="orange",cex=2) 
# Really, to fit a nonlinear curve, one should approach it biologically--what are the underlying phenomena? What is appropraite? See http://stackoverflow.com/questions/21969132/comparing-nonlinear-regression-models. Since this is pretty novel, and there's no obvious model to apply, it seems reasonable to simply look at the data and assess what looks 'good'. The relationship is really ~exponential: at values of NDVI > ~.7, biomass climbs steeply. One can fit an exponential model using nls (y~exp(a+b*x; requires start values for a and b) or log transform the y variable (biomass) to squish down the high values in biomass and high NDVI values...I chose the latter approach due to my familiarity with linear models and interpreting outputs, etc...  
# log transformed
plot(NDVIbiomass1215$SR_WVnir_NDVI,log10(NDVIbiomass1215$Biomass)) # Biomass log transformed to smooth some of the noise...not bad. A linear regression of log(y) is essentially a transformed exponential... 
plot(log10(NDVIbiomass1215$SR_WVnir_NDVI),log10(NDVIbiomass1215$Biomass)) # both log transformed: much better looking, but no gain in r^2 (see below)...

# Compare linear models. Simple linear
a<-lm(NDVIbiomass1215$Biomass~NDVIbiomass1215$SR_WVnir_NDVI);summary(a) # R squared of .64
# log transform biomass (~exponential)
b<-lm(log10(NDVIbiomass1215$Biomass)~NDVIbiomass1215$SR_WVnir_NDVI);summary(b) # R squared of .75
# So, log10(Biomass)=0.7795 + 1.81348(NDVI)
bb<-lm(log10(NDVIbiomass1215$Biomass)~NDVIbiomass1215$eMO_NDVI);summary(bb) #eMODIS is crappy: R squared of 

# Quadratic?
NDVIbiomass1215$NDVI2<-(NDVIbiomass1215$SR_WVnir_NDVI)^2 
quad<-lm(NDVIbiomass1215$Biomass~NDVIbiomass1215$SR_WVnir_NDVI+NDVIbiomass1215$NDVI2);summary(quad) # R squared of .67
# log-log
c<-lm(log10(NDVIbiomass1215$Biomass)~log10(NDVIbiomass1215$SR_WVnir_NDVI));summary(c) # R squared of .74

# Plot model predictions from model b: log10Biomass~NDVI
res<-data.frame(NDVIbiomass1215$SR_WVnir_NDVI,pred=predict(b),na.omit=TRUE)
par(mar=c(5.1,6.1,4.1,2.1))   # adjust margins for extra axis
plot(NDVIbiomass1215$SR_WVnir_NDVI,log10(NDVIbiomass1215$Biomass),axes=FALSE,xlab="NDVI",ylab="Biomass (g)",cex.lab=2,cex.axis=1.75,lwd=2,cex=2,las=1)   # axes=FALSE,
axis(1, at=c(.2,.4,.6,.8), labels=c("0.2","0.4","0.6","0.8"),cex.axis=1.75) 
axis(2, las=1,at=c(1,1.544,2,2.301), labels=c("10","35","100","200"),cex.axis=1.75,mgp=c(.5,.5,0))
points(res$NDVIbiomass1215.SR_WVnir_NDVI,res$pred,type='l',lwd=2) # Plots predicted values
CIs<-predict(b,interval="confidence") # Note, too, that typing 'CIs' will yield the CIs, which could export and back-transform
# Note that this plots 95% confidence interval...."confidence" is for confidence interval
# A bit confusing, but I think worth noting. It seems the latter refers to the mean of the dependent variable, the former the dependent variable itself when predicted at a future time. So, my best est of the slope parameter is shown by the confidence interval, but prediction interval would show my best est of ranges of biomass derived from future NDVI readings...I think! See http://robjhyndman.com/hyndsight/intervals/
lines(NDVIbiomass1215$SR_WVnir_NDVI,CIs[,2],lty=2,lwd=1.5);lines(NDVIbiomass1215$SR_WVnir_NDVI,CIs[,3],lty=2,lwd=1.5) # Adds upper and lower 95% CIs
PIs<-predict(b,interval="prediction")
lines(NDVIbiomass1215$SR_WVnir_NDVI,PIs[,2],lwd=1.5);lines(NDVIbiomass1215$SR_WVnir_NDVI,PIs[,3],lwd=1.5)
text(-.6,2.1,"2012-2015",cex=2)
box(lwd=2)
# OK, what are these predictions like? Back-transform CI and PI...
plot(NDVIbiomass1215$SR_WVnir_NDVI,NDVIbiomass1215$Biomass,xlab="NDVI",ylab="Biomass (g)",cex.axis=1.5,cex.lab=2)
lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(PIs[,1]),lwd=1.5);lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(PIs[,2]),lwd=2,lty=2);lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(PIs[,3]),lwd=2,lty=2) # Prediction interval
lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(CIs[,2]),lwd=1.5,lty=2);lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(CIs[,3]),lwd=1.5,lty=2)# Confidende interval on model... 
text(.25,175,"2012-2015",cex=2);text(.3,150,"Log10 Biomass ~ NDVI",cex=2)

