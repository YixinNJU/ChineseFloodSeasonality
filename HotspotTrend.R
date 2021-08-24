setwd('G:/2=Flood seasonality/StatResults/Hotspots')

library('trend')
library("ggplot2")
library("circular")

rm(list=ls())

####### calculate correlation coefs #########
hsVarCoef = data.frame(PRE = NA, SM = NA, SW = NA)
hsVarCoefSig = hsVarCoef
for(i in 2:4){
  var = read.csv('G:/2=Flood seasonality/StatResults/Hotspots/Hotspot_filt_AMF.csv')
  amf = read.csv(paste0('G:/2=Flood seasonality/StatResults/Hotspots/Hotspot_filt_',fns[i],'.csv'))
  for(j in 1:6){
    v1 = circular(var[1:57,j+1]/366*2*pi,units='radians',zero=0,rotation = c("counter"))
    v2 = circular(amf[1:57,j+1]/366*2*pi,units='radians',zero=0,rotation = c("counter"))
    hsVarCoef[j,i-1] = cor.circular(v1,v2,test = T)$cor
    hsVarCoefSig[j,i-1] = cor.circular(v1,v2,test = T)$p.value
  }
}
write.csv(hsVarCoef,'G:/2=Flood seasonality/StatResults/Hotspots/HSFilteredCoef.csv',row.names = F)  
write.csv(hsVarCoefSig,'G:/2=Flood seasonality/StatResults/Hotspots/HSFilteredCoefSig.csv',row.names = F)  

