# This script is used for trend interpolation

setwd('G:/2=Flood seasonality/StatResults/')

library(readxl)
library(sp)
library(automap)
library(raster)
library(rgdal)
library(ggplot2)
library(latex2exp)

rm(list = ls())

# prepare data
shp = readOGR("G:\\Automap\\WSDIV1\\WSDIV1.shp")
prj = shp@proj4string
WGS84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

Data = read.csv('./Seasonality1257.csv')
input = Data[,c('lon','lat','mmkp','mmk')]
colnames(input) = c("lon","lat","pval",'senMean')

# make sure data is numeric
df = as.data.frame(lapply(input, as.numeric))

#############################################################
## make a SpatialPointDataFrame for data that is to be interpolated.
df = as.data.frame(lapply(df, as.numeric))
rawData = SpatialPointsDataFrame(coords=df[,c("lon", "lat")], 
                                 data=data.frame(df[,"senMean"]),
                                 proj4string=CRS(WGS84))
zd = zerodist(rawData) # find those with the same locations
inData = rawData[-zd[,1],] # remove replicas
inData = rawData
remove(rawData,zd)

# Transform the cooords of inData
inData = spTransform(inData,prj)

# square bound is from shp@bbox
shp@bbox[1,1]
cover = data.frame(x=c(shp@bbox[1,1],shp@bbox[1,1],shp@bbox[1,2],shp@bbox[1,2]),y=c(shp@bbox[2,1],shp@bbox[2,2],shp@bbox[2,2],shp@bbox[2,1]),data=rep("data",4))
coordinates(cover) =~x+y
grid = makegrid(cover,cellsize=20000)
#plot(grid)

# make a country boundary dataframe
krige_grid = data.frame(lon=NA,lat=NA)
for(j in c(1:nrow(shp@data))){
  boundf = data.frame(id =NA, lon=NA,lat=NA)
  n_polys = length(shp@polygons[[j]]@Polygons)
  for(i in c(1:n_polys)){
    shpoly = shp@polygons[[j]]
    tmp_frame = data.frame(id=rep(shpoly@ID,nrow(shpoly@Polygons[[i]]@coords)),
                           lon=shpoly@Polygons[[i]]@coords[,1],
                           lat=shpoly@Polygons[[i]]@coords[,2])
    boundf = rbind(boundf,tmp_frame)
    boundf = na.omit(boundf)
  }
  # only keep points that are inside mainland China and make them a grid
  keep_point = rep(0,nrow(grid))
  ptinChina = point.in.polygon(point.x=grid[,1], point.y=grid[,2], pol.x=boundf[,"lon"], pol.y=boundf[,"lat"], mode.checked=FALSE)
  keep_point[which(ptinChina!=0)] = 1
  gridtmp = grid[which(keep_point==1),]
  names(gridtmp) = c("lon","lat")
  krige_grid = rbind(krige_grid,gridtmp)
  print(j)
}

krige_grid = na.omit(krige_grid)

# transfrom the grid dataframe to a SpatialPointsDataFrame
grid_krig = SpatialPointsDataFrame(coords=krige_grid[,c(1,2)], data=krige_grid[,c(1,2)],proj4string=prj)

# do kriging interpolation
kriging_result <- autoKrige(df....senMean..~1, input_data=inData, new_data = grid_krig)
vgrm = autofitVariogram(df....senMean..~1, input_data=inData)
plot(vgrm)
vgrm
# plot them
#automapPlot(kriging_result$krige_output,"var1.pred")

# get results in dataframe
meanSendf = data.frame(kriging_result$krige_output@coords, 
                       pred = kriging_result$krige_output@data$var1.pred) 
MSvar = data.frame(kriging_result$krige_output@coords, 
                   var = kriging_result$krige_output@data$var1.var)
MSstdev = data.frame(kriging_result$krige_output@coords, 
                     stdev = kriging_result$krige_output@data$var1.stdev) 

# transfrom the projs
DFtoWGS = function(df,val){
  df = SpatialPointsDataFrame(coords=df[,c("lon", "lat")], 
                              data=data.frame(df[,c(val)]),
                              proj4string=prj)
  newprj = WGS84
  ndf = spTransform(df,CRS(newprj))
  ndf = data.frame(ndf@coords, 
                   val = ndf@data[[1]]) 
  names(ndf) = c("lon","lat",val)
  return(ndf)
}

MSpred = DFtoWGS(meanSendf,"pred")
MSvar = DFtoWGS(MSvar,"var")
MSstdev = DFtoWGS(MSstdev,"stdev")
meanSensOK = cbind(MSpred,MSvar$var,MSstdev$stdev)
names(meanSensOK) = c("lon","lat","pred","var","stdev")

write.csv(meanSensOK,'G:/2=Flood seasonality/KrigeResults/OKcsv/OK1257.csv',row.names = F)
