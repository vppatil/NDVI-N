#emodis_reshape

setwd('~/Data/')
library(reshape2)
emo.wide<-read.csv('emo_wide.csv')
emo.wide=emo.wide[!is.na(emo.wide$Julian.Date),]
emo.names=grep('eMODIS_NDVI',names(emo.wide),value = TRUE)

emo.wide=subset(emo.wide,select=c('Julian.Date',emo.names))

emo.reshape=melt(emo.wide,measure.vars=emo.names)

emo.reshape$variable<-as.character(emo.reshape$variable)
emo.reshape$Year=sapply(emo.reshape$variable,function(x) as.numeric(gsub('X','',strsplit(x,split='_')[[1]][1]))+2000)
emo.reshape=na.omit(emo.reshape)
names(emo.reshape)=c('DOY','variable','eMO_NDVI','Yr')
emo.reshape=subset(emo.reshape,select=c('DOY','eMO_NDVI','Yr'))

write.csv(emo.reshape,'eMODIS_long.csv')

###
lonely.emo<-read.csv('lonely_emodis_long.csv')
lonely.emo$Date<-as.character(lonely.emo$Date)
vars=c('Pixel1','Pixel2','Pixel3')
lonely.emo=melt(lonely.emo,measure.vars=vars)
names(lonely.emo)
lonely.emo.dates=t(sapply(lonely.emo$Date,function(x) strsplit(x,split='/')[[1]]))
lonely.emo$Month=as.numeric(lonely.emo.dates[,1])
lonely.emo$Day=as.numeric(lonely.emo.dates[,2])
lonely.emo$Year<-as.numeric(lonely.emo$Year)
names(lonely.emo)[5:6]=c('Plot','eMO_NDVI')

lonely.emo$DOY<-julian.df(lonely.emo)
write.csv(lonely.emo,'lonely_emo_transformed.csv')
