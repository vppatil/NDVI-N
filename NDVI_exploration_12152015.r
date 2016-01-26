setwd('D:/aksciencecenter files from drive/Data/')
dir()
library(XLConnect)

vegplots<-read.csv('VegPlots_for R.csv')

biomass.ndvi<-subset(vegplots,select=c('Biomass','eMO_NDVI','SR_eMO_NDVI','SR_WVre_NDVI','SR_WVnir_NDVI'))
biomass.ndvi<-biomass.ndvi[!is.na(biomass.ndvi$Biomass),]
biomass.ndvi.impute<-rfImpute(Biomass~.,data=biomass.ndvi)

biomass.ndvi.rf1<-randomForest(Biomass~.,data=biomass.ndvi.impute)
partialPlot(biomass.ndvi.rf1,biomass.ndvi.impute,'SR_eMO_NDVI') #partial dependence plot.

prcntN.ndvi<-subset(vegplots,select=c('PrcntN','eMO_NDVI','SR_eMO_NDVI','SR_WVre_NDVI','SR_WVnir_NDVI'))
prcntN.ndvi<-prcntN.ndvi[!is.na(prcntN.ndvi$PrcntN),]
prcntN.ndvi.impute<-rfImpute(PrcntN~.,data=prcntN.ndvi)

prcntN.ndvi.rf1<-randomForest(PrcntN~.,data=prcntN.ndvi.impute,mtry=2,ntree=500)
varImpPlot(prcntN.ndvi.rf1)
partialPlot(prcntN.ndvi.rf1,prcntN.ndvi.impute,'eMO_NDVI') #partial dependence plot.

plot(PrcntN~eMO_NDVI)

a<-lm(PrcntN~I(eMO_NDVI^2)+eMO_NDVI)

#nls try
nl1<-nls(Biomass~exp(a+b*SR_WVnir_NDVI),data=vegplots)


#building a matrix of 