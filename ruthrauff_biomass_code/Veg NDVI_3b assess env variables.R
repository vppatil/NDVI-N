# Colville Veg / NDVI model
# assuming that log10(Biomass)~NDVI is our best underlying fit, do environmental variables help predict biomass?
# Load packages: AICcmodavg and lme4
library(AICcmodavg)
library(lme4)
library(car)
library(psych)

# Set working directory
setwd("C:/Data/Dan/Changing Arctic Ecosystems/NDVI summaries Kyle 2015")

# Import working file
ColVeg<-read.csv(file = "VegPlots_for R.csv")

# Verify the structure of the database
str(ColVeg);head(ColVeg)
# Create a variable combining plot ID and Exclosure/control....
ColVeg$Union<-paste(ColVeg$SiteID,ColVeg$Ex_Co,sep="_");
# Define Union as factor:
ColVeg$Union<-factor(ColVeg$Union);str(ColVeg)
# Define year as factor:
ColVeg$Yr<-factor(ColVeg$Yr);str(ColVeg);head(ColVeg)
# Create a variable combining Union and year to use as random variable ....
ColVeg$UnionYr<-paste(ColVeg$Union,ColVeg$Yr,sep="_");str(ColVeg)
# Define UnionYr as factor:
ColVeg$UnionYr<-factor(ColVeg$UnionYr);str(ColVeg)

# NDVI as a function of year: plot
plot(ColVeg$DOY,ColVeg$SR_WVnir_NDVI)
# Remove 2011--protocols not standardized, collections not meticulous
ColVeg1215<-ColVeg[which(ColVeg$Yr!='2011'),];str(ColVeg1215)
# Drop unused levels
drop.levels <- function(dat){
  # Drop unused factor levels from all factors in a data.frame
  # Author: Kevin Wright.  Idea by Brian Ripley. Found on http://rwiki.sciviews.org/doku.php?id=tips:data-manip:drop_unused_levels
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}
ColVeg1215<-drop.levels(ColVeg1215)
levels(ColVeg1215$Yr)

# Predicting biomass with NDVI, environmental variables....
# From 'Veg NDVI_1 data exploration.r' I determined which environmental variables were 'best' (accounted for most variance of other variables)
# Simple linear model? logBiomass improves fit...temp and lumin log-transformed b/c on much different scale--returned warnings otherwise...assumes no effect of treatment, so each plot each year (ie, 'UnionYr') are unique and modeled as random effect
m<-list()
(m[[1]]<-lmer(log10(Biomass)~SR_WVnir_NDVI + log10(Soil_TDDcum) + log10(Lumin_Cum) + Precip_CumPeriod + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215)) # Model each plot (Plot + treatment + year; UnionYr) as random effect
(m[[2]]<-lmer(log10(Biomass)~SR_WVnir_NDVI +  log10(Lumin_Cum) + Precip_CumPeriod + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[3]]<-lmer(log10(Biomass)~SR_WVnir_NDVI +  log10(Soil_TDDcum) + Precip_CumPeriod + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[4]]<-lmer(log10(Biomass)~SR_WVnir_NDVI +  log10(Soil_TDDcum) + log10(Lumin_Cum) + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[5]]<-lmer(log10(Biomass)~SR_WVnir_NDVI +  log10(Soil_TDDcum) + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[6]]<-lmer(log10(Biomass)~SR_WVnir_NDVI +  log10(Lumin_Cum) + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[7]]<-lmer(log10(Biomass)~SR_WVnir_NDVI +  Precip_CumPeriod + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[8]]<-lmer(log10(Biomass)~SR_WVnir_NDVI + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[9]]<-lmer(log10(Biomass)~SR_WVnir_NDVI + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[10]]<-lmer(log10(Biomass)~1 + (1|UnionYr), REML=FALSE, data=ColVeg1215))
(m[[11]]<-lmer(log10(Biomass)~log10(Soil_TDDcum) + log10(Lumin_Cum) + Precip_CumPeriod + DOY + (1|UnionYr), REML=FALSE, data=ColVeg1215))

# Assign names to each model....
Modnames<-c("1: Full","2: NDVI + Lumin + Pre","3: NDVI + Temp + Pre","4: NDVI + Temp + Lumin","5: NDVI + Temp",
"6: NDVI + Lumin","7: NDVI + Precip","8: NDVI","9: NDVI + DOY","10: Intercept","11: Environmental")

# create AIC table, sorted by DeltaAIC...
aictab(cand.set = m, modnames = Modnames, sort = TRUE) # Model 2 .99 model wt.; #6 next at deltaAIC=8.4

# Model-averaged parameters...
print(modavg(parm="SR_WVnir_NDVI", cand.set = m, modnames = Modnames), digits =4) # Positive; non-zero
print(modavg(parm="log10(Soil_TDDcum)", cand.set = m, modnames = Modnames), digits =4) # Overlaps zero
print(modavg(parm="log10(Lumin_Cum)", cand.set = m, modnames = Modnames), digits =4) # Positive; non-zero
print(modavg(parm="Precip_CumPeriod", cand.set = m, modnames = Modnames), digits =4) # Positive, non-zero
print(modavg(parm="DOY", cand.set = m, modnames = Modnames), digits =4) # Negative, overlaps zero 
print(modavg(parm="(Intercept)", cand.set = m, modnames = Modnames), digits =4) # Positive, non-zero 

# How do these variables improve point estimates above and beyond the simple linear model including NDVI alone?
new.dat<-data.frame(SR_WVnir_NDVI= .6 ,Soil_TDDcum = 100,Lumin_Cum = 50200,Precip_CumPeriod = 3,DOY = 180)
print(modavgPred(cand.set = m, modnames = Modnames, newdata = new.dat, uncond.se = "revised"),digits=4)
# Yields a model-averaged prediction of 1.825, which is 10^1.825 = 66.83 g. The values for the environmental variables were crudely estimated off plots of those variables across season for DOY=180. Using model from 'Veg NDVI_3 fit linear model to predict biomass.r', a NDVI value of .6 yields a biomass estimate of 1.8676, which is 10^1.8676=73.7 g. Not much different.

# How does temp, say, affect this? Soil TDD peaks at ~320.....
new.dat<-data.frame(SR_WVnir_NDVI= .6 ,Soil_TDDcum = 320,Lumin_Cum = 50200,Precip_CumPeriod = 3,DOY = 180)
print(modavgPred(cand.set = m, modnames = Modnames, newdata = new.dat, uncond.se = "revised"),digits=4) # No influence: 1.825 
# How does luminosity affect this? Lumin_Cum peaks at ~200,000.....
new.dat<-data.frame(SR_WVnir_NDVI= .6 ,Soil_TDDcum = 100,Lumin_Cum = 200000,Precip_CumPeriod = 3,DOY = 180)
print(modavgPred(cand.set = m, modnames = Modnames, newdata = new.dat, uncond.se = "revised"),digits=4) # 1.859, which equals an estimate of 72.3 g. So, really it still is mostly pushed by NDVI alone.
results<-data.frame(modavgPred(cand.set = m, modnames = Modnames, newdata = new.dat, uncond.se = "revised")) 


# Estimate mean values of temp, precip, and lumen for plotting CIs
describe(ColVeg1215$Soil_TDDcum);plot(ColVeg1215$DOY,ColVeg1215$Soil_TDDcum) #198.77; increases linearly across season
describe(ColVeg1215$Precip_CumPeriod);plot(ColVeg1215$DOY,ColVeg1215$Precip_CumPeriod)  #5.02; slight increase across season, but pretty variable
describe(ColVeg1215$Lumin_Cum);plot(ColVeg1215$DOY,ColVeg1215$Lumin_Cum) #112,232.2; increases linearly across season
describe(ColVeg1215$DOY) #192.65...
new.dat<-data.frame(SR_WVnir_NDVI= seq(from=.2, to=.83,by=0.01) ,Soil_TDDcum = seq(from=.5,to=529.7,by=8.4),Lumin_Cum = seq(from=4500,to=226260,by=3520),Precip_CumPeriod = seq(from=0,to=20.664,by=.328),DOY = seq(from=161,to=224,by=1)) # Uses min and max values for environmental variables, increasing linearly across season with NDVI, and DOY increases across season in step with NDVI
avPred<-data.frame(modavgPred(cand.set = m, modnames = Modnames, newdata = new.dat, uncond.se = "revised"))
# merge with new.dat so can plot below; create ID in both data frames
new.dat$ID<-seq(1:64);head(new.dat);avPred$ID<-seq(1:64);head(avPred)
avPred1<-merge(new.dat,avPred,by="ID");head(avPred1) 


######################################################################################################
# compare confidence and prediction intervals between simple logBiomass~NDVI model and model 2 above #
######################################################################################################
# Will need to run 'Veg NDVI_3 etc etc' to get working dataset, NDVIbiomass1215...
b<-lm(log10(NDVIbiomass1215$Biomass)~NDVIbiomass1215$SR_WVnir_NDVI);summary(b)
# Plot model predictions from model b: log10Biomass~NDVI
res<-data.frame(NDVIbiomass1215$SR_WVnir_NDVI,pred=predict(b),na.omit=TRUE)
par(mar=c(5.1,6.1,4.1,2.1))   # adjust margins for extra axis
plot(NDVIbiomass1215$SR_WVnir_NDVI,log10(NDVIbiomass1215$Biomass),axes=FALSE,xlab="NDVI",ylab="Biomass (g)",cex.lab=2,cex.axis=1.75,lwd=2,cex=2,las=1)   # axes=FALSE,
axis(1, at=c(.2,.4,.6,.8), labels=c("0.2","0.4","0.6","0.8"),cex.axis=1.75) 
axis(2, las=1,at=c(1,1.544,2,2.301), labels=c("10","35","100","200"),cex.axis=1.75,mgp=c(.5,.5,0))
points(res$NDVIbiomass1215.SR_WVnir_NDVI,res$pred,type='l',lwd=2) # Plots predicted values
CIs<-predict(b,interval="confidence") # Note, too, that typing 'CIs' will yield the CIs, which could export and back-transform
lines(NDVIbiomass1215$SR_WVnir_NDVI,CIs[,2],lty=2,lwd=1.5);lines(NDVIbiomass1215$SR_WVnir_NDVI,CIs[,3],lty=2,lwd=1.5) # Adds upper and lower 95% CIs
PIs<-predict(b,interval="prediction")
lines(NDVIbiomass1215$SR_WVnir_NDVI,PIs[,2],lwd=1.5);lines(NDVIbiomass1215$SR_WVnir_NDVI,PIs[,3],lwd=1.5)
text(-.6,2.1,"2012-2015",cex=2)
box(lwd=2)
###############################
## Add model-averaged predictions
###############################
# Note: I tried to add prections derived only from model 2 using same approach as above, but the predictions were not smooth b/c for a given NDVI observation, temp, precip, and lumin all varied enough to make the predictions a huge mess...so, I think better to run predictions over a range of environmental values as above to see how sensitive the model is to these variations....this represents output integrating values increasing linearly across season
lines(avPred1$SR_WVnir_NDVI,avPred1$mod.avg.pred,lwd=2,col="red") # Prediction estimate in red
lines(avPred1$SR_WVnir_NDVI,(avPred1$mod.avg.pred+(1.96*avPred1$uncond.se)),lty=2,lwd=1.5,col="red") # Adds upper 95% CI
lines(avPred1$SR_WVnir_NDVI,(avPred1$mod.avg.pred-(1.96*avPred1$uncond.se)),lty=2,lwd=1.5,col="red") # Adds lower 95% CI


# OK, what are these predictions like? Back-transform CI and PI...
par(mar=c(5.1,6.1,4.1,2.1))   # adjust margins 
par(mgp=c(4,1,0)) # Bump 'Biomass' label out a bit...
plot(NDVIbiomass1215$SR_WVnir_NDVI,NDVIbiomass1215$Biomass,xlab="WorldView-2 NDVI",ylab="Biomass (g)",cex.axis=1.75,cex.lab=2,cex=1.5,lwd=1.5,las=1)
lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(PIs[,1]),lwd=2);# lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(PIs[,2]),lwd=2,lty=2);lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(PIs[,3]),lwd=2,lty=2) # Prediction interval
lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(CIs[,2]),lwd=1.5,lty=2);lines(NDVIbiomass1215$SR_WVnir_NDVI,10^(CIs[,3]),lwd=1.5,lty=2)# Confidence interval on model...

lines(avPred1$SR_WVnir_NDVI,10^(avPred1$mod.avg.pred),lwd=2,lty=6) # Back-transformed prediction estimate bold, dot-dash
# lines(avPred1$SR_WVnir_NDVI,10^(avPred1$mod.avg.pred),lwd=2,col="red") # Back-transformed prediction estimate in red
lines(avPred1$SR_WVnir_NDVI,10^(avPred1$mod.avg.pred+(1.96*avPred1$uncond.se)),lty=2,lwd=1.5,col="red") # Adds upper 95% CI
lines(avPred1$SR_WVnir_NDVI,10^(avPred1$mod.avg.pred-(1.96*avPred1$uncond.se)),lty=2,lwd=1.5,col="red") # Adds lower 95% CI
 
text(.25,175,"2012-2015",cex=2)


 

