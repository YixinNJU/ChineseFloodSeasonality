# This script is used to do DBSCAN for flood and precipitation (i.e., rainfall)

setwd('G:/2=Flood seasonality/Scale_DBSCAN/')
Sys.setenv(LANGUAGE = 'en')

library('fpc')
library('ggplot2')
library('cluster')
library('readxl')
library('geometry')
library('rgdal')
library('sp')
library('sf')
library('spatialEco')
library('trend')
library('raster')
library("dbscan")
library(segmented)
library(changepoint)
library(concaveman)

########### Flood DBSCAN ###########
rm(list=ls())

wd = 7 # time window

# prepare data
Dates = read.csv('G:/2=Flood seasonality/BufferAnalysis/DateRecords_1960.csv')
Dates = Dates[,1:59]
colnames(Dates) = c('GaugeID',as.character(1960:2017))

Data = readOGR('G:/2=Flood seasonality/StationShp/S1257_EASM.shp')
crs = Data@proj4string

Data = Data@data
id = read.csv('G:/2=Flood seasonality/BufferAnalysis/Sync1967_2017_part/ID596.csv')
Data = merge(id,Data,by=1)
df = Data[,c('GaugeID','x','y','lon','lat')]
df = merge(df,Dates,by=1)

# DBSCAN

set.seed(10000)

yrScale = data.frame(year = 1960:2017, scale = NA,tw = NA)
for(year in 1960:2017){
  print(year)
  colName = as.character(year)
  yrLayer = df[,c('x','y',colName)]
  yrLayer = na.omit(yrLayer) # using {spatialEco} package
  yrLayer[yrLayer[colName] == 1,colName] = NA
  yrLayer = na.omit(yrLayer)
  
  startDay = min(yrLayer[colName])-1
  endDay = max(yrLayer[colName])
  
  tw_Area = vector('double',length(startDay:endDay))
  tw_k = vector('double',length(startDay:endDay))
  tw_staNum = vector('double',length(startDay:endDay))
  tw_eps = tw_k
  
  for(day in startDay:endDay){
    tw = c(day-wd,day+wd) # time window
    tmp = yrLayer[(yrLayer[,colName]<= tw[2] & yrLayer[,colName]> tw[1]),]
    
    xy = tmp[,c('x','y')] 
    
    if(nrow(xy) < 20){ # at least 30 stations
      next()
    }
    
    minPts = 10
    
    #png(paste0('./tmpFigs/',day,'.png'))
    kdis = kNNdist(xy,k = minPts)
    kdis_sort = sort(kdis,decreasing = F)
    eps = 2e5
    #plot(kdis_sort)
    #abline(h = eps, lty = 2)
    dbscan::kNNdistplot(xy, k =  minPts)
    abline(h = eps, lty = 2)
    #dev.off()
    
    clusResult = dbscan::dbscan(xy, eps=eps, minPts = minPts)
    
    xy$label = clusResult$cluster # give each point its label
    
    k = max(unique(clusResult$cluster))
      
    # get convex hull for each cluster
    plys = vector("list", k)
    count = 1
    plyArea = vector('double',k)
    for(i in 1:k){ # loop for each cluster
      pts = xy[xy$label == i,c('x','y')]
      if(nrow(pts) < 10){ # there are less than 10 points in this cluster
        next()
      }
      ## convex polygons
      hpts <- chull(pts)
      hpts <- c(hpts, hpts[1])
      p = Polygon(pts[hpts, ], hole=F)
      plys[count] = p
      count = count+1
      
      plyArea[i] = p@area/10^6 # in sqkm
    }

    #pu = st_union(p)
    tw_Area[day-startDay+1] = sum(na.omit(plyArea))
    tw_k[day-startDay+1] = k
    tw_staNum[day-startDay+1] = nrow(xy)
    tw_eps[day-startDay+1] = eps
    }
  
  twYrArea = data.frame(day = startDay:endDay,twArea = tw_Area,knum = tw_k,snum = tw_staNum,eps = tw_eps)
  twYrArea = na.omit(twYrArea)
  area = tw_Area[tw_Area!=0]
  yrScale[year-1959,'scale'] = max(area)
  yrScale[year-1959,'tw'] = twYrArea[which.max(twYrArea$twArea),'day']
  yrScale[year-1959,'knum'] = twYrArea[which.max(twYrArea$twArea),'knum']
  yrScale[year-1959,'snum'] = twYrArea[which.max(twYrArea$twArea),'snum']
  yrScale[year-1959,'eps'] = twYrArea[which.max(twYrArea$twArea),'eps']
}

sens.slope(yrScale$scale)
hist(yrScale$eps)

write.csv(yrScale,'G:/2=Flood seasonality/Scale_DBSCAN_3d/maxYear_10nn_2e5/SumMaxScale3d_dbscan_10nn_2e5.csv',row.names = F)

############ Rain DBSCAN ##############

rm(list=ls())

df = read.csv('G:/2=Flood seasonality/Scale_Kmeans/PR_NC/PR_Dates.csv') 
colnames(df) = c('x','y','lon','lat',as.character(1961:2015))

df = df[df$lat <= 40 & df$lat >= 25,]

## get each cluster polygon area
set.seed(10000)

## stat for the maximum cluster number
Knum = matrix(data = NA, nrow = length(1961:2015), 
              ncol = length(startDay:endDay))
yearK = vector(mode = 'double',length = length(1961:2015))

for(year in 1961:2015){
  print(year)
  colName = as.character(year)
  yrLayer = df[,c('x','y',colName)]
  yrLayer = na.omit(yrLayer) # using {spatialEco} package
  
  tw_k = vector('double',length(startDay:endDay))
  startDay = min(yrLayer[colName]) -1
  endDay = max(yrLayer[colName])
  for(day in startDay:endDay){
    tw = c(day-wd,day+wd) # time window
    tmp = yrLayer[(yrLayer[,colName]<= tw[2] & yrLayer[,colName]> tw[1]),]
    
    xy = tmp[,c('x','y')] 
    
    if(nrow(xy) < 20){ # at least 30 stations
      next()
    }
    
    minPts = 10
    
    
    kdis = kNNdist(xy,k = minPts)
    kdis_sort = sort(kdis,decreasing = F)
    eps = 2e5
    
    #png(paste0('G:/2=Flood seasonality/Scale_DBSCAN/PR/10nn_plot/',day,'.png'))
    #dbscan::kNNdistplot(xy, k =  minPts)
    #abline(h = eps, lty = 2)
    #dev.off()
    
    clusResult = dbscan::dbscan(xy, eps=eps, minPts = minPts)
    
    xy$label = clusResult$cluster # give each point its label
    
    k = max(unique(clusResult$cluster))
    
    Knum[year-1960,day-startDay+1] = k
  }
  yearK[year-1960]= sum(na.omit(Knum[year-1960,]))
}

############### Export shapefiles and points for plot##################
getYearTw = function(df,year,tw){
  print(year)
  colName = as.character(year)
  yrLayer = df[,c('x','y',colName)]
  yrLayer = na.omit(yrLayer) # using {spatialEco} package
  
  tw = c(day-wd,day+wd) # time window
  tmp = yrLayer[(yrLayer[,colName]<= tw[2] & yrLayer[,colName]> tw[1]),]
  
  xy = tmp[,c('x','y')]
  return(xy)
}

sens.slope(yrScale$scale)
## get save xy data for specific maximum(median) years
for(i in 1:nrow(yrScale)){
  year = yrScale$year[i]
  scale = yrScale$scale[i]
  day = yrScale$tw[i]
  kn = yrScale$knum[i]
  
  ## write xy data to output
  xy = getYearTw(df,year = year,tw = day)
  
  ## do cluster for this year
  tw = c(day-wd,day+wd)
  
  minPts = 4
  
  ### plot knn
  png(paste0('./test//4nn_',year,'_',day,'.png'))
  
  kdis = kNNdist(xy,k = minPts)
  kdis_sort = sort(kdis,decreasing = F)
  eps = 2e5
  dbscan::kNNdistplot(xy, k =  minPts)
  abline(h = eps, lty = 2)
  dev.off()
  # end plot knn
  
  # do cluster 
  clusResult = dbscan::dbscan(xy, eps=eps, minPts = minPts)
  
  xy$label = clusResult$cluster # give each point its label
  write.csv(xy,paste0('./test/',year,'_',day,'.csv'),row.names = F)
  
  k = max(unique(clusResult$cluster))
  
  # get convex hull for each cluster
  hnum = 1
  for(i in 1:k){ # loop for each cluster
    pts = xy[xy$label == i,c('x','y')]
    if(nrow(pts) < 4){ # there are less than 3 points in this cluster
      next()
    }
    
    #hpts <- chull(pts)
    #hpts <- c(hpts, hpts[1])
    #p = Polygon(pts[hpts, ], hole=F)
    
    # concave poly
    ps = st_as_sf(SpatialPoints(pts))
    ccv = concaveman(ps)
    p = Polygon(ccv[[1]][[1]][[1]])
    
    ply = SpatialPolygons(list(Polygons(list(p),ID='A')))
    shapefile(x = ply, file = paste0("./test//AMF_",year,'_',day,'_',hnum,'.shp'),overwrite=TRUE)
    hnum = hnum+1
  }
}