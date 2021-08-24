setwd('G:/2=Flood seasonality/StatResults/Hotspots')

library('trend')
library("ggplot2")
library("circular")

rm(list=ls())

DateData = read.csv('G:/2=Flood seasonality/BufferAnalysis/DateRecords_1960.csv')

###### flood hotspots #########
for(j in 1:6){
  fn = paste0('G:\\2=Flood seasonality\\StationShp\\HotspotShp6\\HS',j,'.shp')
  hotspot = readOGR(fn)
  hotspot = as.data.frame(hotspot)
  hotspot = hotspot[,c('GaugeID','area','lon','lat','mmkp','mmk','x','y')]
  colnames(hotspot) = c('GaugeID','area','lon','lat','pval','Sen','x','y')
  #hotspot = hotspot[which(as.numeric(as.character(hotspot$area)) <= 50000),]
  #hotspot = merge(hotspot,dams_big,by='GaugeID')
  Dates = merge(hotspot,DateData,by.x = 'GaugeID',by.y='V1')
  
  Smin = min(Dates$Sen)
  Smax = max(Dates$Sen)
  ###### find medium of each year ########
  meds = data.frame(year = NA, med = NA)
  for(i in 2:59){# 1960 - 2016
    rname = paste0('V',as.character(i))
    yrData = na.omit(Dates[,rname])
    meds[i-1,'year'] = 1960+i-2
    meds[i-1,'med'] = round(mean(yrData))
  }
  if(j == 1){
    output = meds
  }
  else{
    colnames(meds) = c('Year',paste0('Hotspot',j))
    output = cbind(output,meds[,2])
  }
}
colnames(output) = c('year','box1','box2','box3','box4','box5','box6')
write.csv(output,'./Hotspot_AMF.csv',row.names = F)