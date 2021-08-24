# This script is used to derive seasonality statistics 
# and trend in Chinese flood seasonality

setwd("G:\\2=Flood seasonality\\StatResults\\")

library("trend")
library("mgcv")
library("nlme")
library("ismev")
library("ggplot2")
library("pracma")
library("changepoint")
library("zoo")
library("Rcpp")
library("readxl")
library("ecp")
library("mgcv")
library("nlme")
library("ismev")
library("ggplot2")
library("circular")
library("lawstat")
library("diptest")
library("readxl")
library("CircStats")

rm(list = ls())

FILELIST = list.files("G:\\2=Flood seasonality\\DateByStation",'*.dat')
FILELIST_PEAK = list.files("G:\\7=PeakData\\PeakByStation\\",'*.dat')

IDXY = read_xlsx('G:\\Final\\ChinaWatersheds1617.xlsx')
IDXY = IDXY[,c('GaugeID','lon','lat','area')]

### statistics for trend and seasonality
Data=data.frame(GaugeID = NA,len1=NA,len2=NA,
                mmk=NA,mmkp=NA, # trend in mean
                md=NA, mrl=NA, kuiper=NA# seasonality  
)

for(i in (1:length(FILELIST))){
  GaugeID = as.numeric(substr(x = FILELIST[i], start = 1, stop = 9))
  if(length(which(IDXY$GaugeID == GaugeID)) == 0){next()} # the station does not have a watershed
  
  file1=paste("G:\\2=Flood seasonality\\DateByStation\\",FILELIST[i],sep='')
  dates=read.table(file1, header = FALSE, sep = ",")
  dates = dates[which(dates$V1 >= 1960),]
  NewDates = dates[which(dates$V2 != 1 & dates$V2 != 365 & dates$V2 != 366),]
  
  # make sure all dates left have real AMFs
  if(length(intersect(FILELIST[i],FILELIST_PEAK)) != 0){ 
    file2=paste("G:\\7=PeakData\\PeakByStation\\",FILELIST[i],sep='')
    peaks =read.table(file2, header = FALSE, sep = ",")
    peaks = peaks[which(peaks$V1 >= 1960),]
    NewPeaks = peaks[!is.na(peaks$V2),]
    Year = intersect(NewDates$V1,NewPeaks$V1)
    NewDates = NewDates[NewDates$V1 %in% Year,]
  }
  
  if(nrow(NewDates) < 30){next()} # record length must more than 50 years
  
  Data[i,'GaugeID'] <- as.integer(GaugeID)
  Data[i,'len1'] <- nrow(dates) # how many years after 1960
  Data[i,'len2'] <- nrow(NewDates) # real length
  
  Data[i,'mmk'] <- sens.slope(NewDates$V2)$estimates*10  # sen's slope
  Data[i,'mmkp'] <- mk.test(NewDates$V2)$p.value
  
  ## seasonality
  cc <- circular(NewDates$V2/366*2*pi, units='radians',zero=0,rotation = c("counter"))
  
  # mean date
  Data[i,'md'] <- mean.circular(cc)/2/pi*366
  if (Data[i,'md']<0){
    Data[i,'md']=Data[i,'md']+366
  }
  # mean resultant length
  Data[i,'mrl'] <- rho.circular(cc)
  
  if(Kuipertest$statistic<=1.747){
    Data[i,'kuiper']= 1 # this means the data series is uniform, can not reject null hypo
  }
  else{Data[i,'kuiper']= 0}
}
Data = na.omit(Data)

Data = as.data.frame(Data)
Data = merge(Data,IDXY,by='GaugeID')

tmp40 = Data[Data$len2 > 40,]
#tmp45 = Data[Data$len2 >= 45,]
#tmp50 = Data[Data$len2 >= 50,]

write.csv(tmp40,paste0('./Seasonality',as.character(nrow(tmp40)),'.csv'),row.names = FALSE)
#write.csv(tmp50,paste0('./Seasonality',as.character(nrow(tmp50)),'.csv'),row.names = FALSE)
