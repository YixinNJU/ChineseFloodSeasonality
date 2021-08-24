%% cut zonal bands from statistics (mean occurrence date and mean trend) of AMPD & MSMD

clc;
clear;

cn = shaperead('G:/2=Flood seasonality\EASMshp/EASMDIS_clip_wgs.shp');

lon  = ncread('G:/CPCSoilMoisture\SoilMoisture/MSMD.nc','lon');
lat  = ncread('G:/CPCSoilMoisture\SoilMoisture/MSMD.nc','lat');
Dates = ncread('G:/CPCSoilMoisture\SoilMoisture/MSMD.nc','MSMD');

%% create EASM region mask
xx = lon;
yy = lat;

[X,Y] = meshgrid(xx,yy);

rx = cn.X(1:end-1);
ry = cn.Y(1:end-1);
     
mask = inpolygon(X,Y,rx,ry);

writematrix(mask,'G:\2=Flood seasonality\StationShp\SM_Pre_Bands/SMask.txt')
%%  check mask
mask = readmatrix('G:\2=Flood seasonality\StationShp\SM_Pre_Bands/SMask.txt');
figure;
imagesc((Dates(:,:,1)));
figure;
imagesc(flipud(rot90(Dates(:,:,1))));
figure;
mask = flipud(rot90(mask));
%mask  = mask';
imagesc(mask);
%DatesInBand = Dates(mask);

%% use mask to clip MSMD (58 year soil moisture timing)
mask = readmatrix('G:\2=Flood seasonality\StationShp\SM_Pre_Bands/SMask.txt');
matsize = size(Dates);
medDate = ones(matsize(1),matsize(2))*nan;
meanDate = ones(matsize(1),matsize(2))*nan;
dateTrend = ones(matsize(1),matsize(2))*nan;

for i = 1:matsize(1)
    for j = 1:matsize(2)
        yr = 1960:2017;
        ts = reshape(Dates(i,j,:),matsize(3),1)/366*2*pi; % change dim of data
        ts = ts(~isnan(ts));
        goodYr = yr(~isnan(ts)); 
        if ~isempty(ts)
            ts = movmean(ts,10);
            medDate(i,j) = circ_median(ts,1)/2/pi*366; % 1 means along this dimension
            meanDate(i,j) = circ_mean(ts,[],1)/2/pi*366;  % do not forget to convert into DOY data
            dateTrend(i,j) = Theil_Sen_Regress(goodYr',ts)*10;
        end
    end
end

% incase there are negative circular data
medDate(medDate<0) = medDate(medDate<0)+366;
meanDate(meanDate<0) = meanDate(meanDate<0)+366;
figure;
imagesc(dateTrend);
colorbar

figure;
imagesc((mask));
colorbar

% set those outside EASM and BASIN region to NAN
medDate(~(mask)) = nan;
meanDate(~(mask)) = nan;
dateTrend(~(mask)) = nan;


figure;
imagesc(dateTrend);
colorbar
%% write annual maximum precipitation date data to output
%AMPD = netcdf.create('./SoilMoisture/MSMDMeanMed.nc','NC_WRITE');
AMPD = netcdf.create('G:/2=Flood seasonality\StationShp\SM_Pre_Bands/SM_stats_EASM.nc','NC_WRITE');

dimlon = netcdf.defDim(AMPD,'lon',length(lon));
dimlat = netcdf.defDim(AMPD,'lat',length(lat));
var1 = netcdf.defVar(AMPD,'lon','NC_DOUBLE',dimlon);
var2 = netcdf.defVar(AMPD,'lat','NC_DOUBLE',dimlat);
var3 = netcdf.defVar(AMPD,'mean','NC_DOUBLE',[dimlon,dimlat]);
var4 = netcdf.defVar(AMPD,'median','NC_DOUBLE',[dimlon,dimlat]);
var5 = netcdf.defVar(AMPD,'trend','NC_DOUBLE',[dimlon,dimlat]);
netcdf.endDef(AMPD);

netcdf.putVar(AMPD,var1,lon);
netcdf.putVar(AMPD,var2,lat);
netcdf.putVar(AMPD,var3,meanDate);
netcdf.putVar(AMPD,var4,medDate);
netcdf.putVar(AMPD,var5,dateTrend);
netcdf.close(AMPD);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar for precipitation
cn = shaperead('G:/2=Flood seasonality\EASMshp/EASMDIS_clip_wgs.shp');

lon  = ncread('G:/PrecipGrid/AMPD.nc','lon');
lat  = ncread('G:/PrecipGrid/AMPD.nc','lat');
Dates = ncread('G:/PrecipGrid/AMPD.nc','AMPD');

xx = lon;
yy = lat;

[X,Y] = meshgrid(xx,yy);

rx = cn.X(1:end-1);
ry = cn.Y(1:end-1);
     
mask = inpolygon(X,Y,rx,ry);

writematrix(mask,'G:\2=Flood seasonality\StationShp\SM_Pre_Bands/RainMask.txt')
%DatesInBand = Dates(mask);

%% check mask
figure;
imagesc((Dates(:,:,1)));
figure;
mask = flipud(rot90(mask));
imagesc(mask);

%%
clc;clear;
cn = shaperead('G:/2=Flood seasonality\EASMshp/EASMDIS_clip_wgs.shp');

lon  = ncread('G:/PrecipGrid/AMPD.nc','lon');
lat  = ncread('G:/PrecipGrid/AMPD.nc','lat');
Dates = ncread('G:/PrecipGrid/AMPD.nc','AMPD');

mask = readmatrix('G:\2=Flood seasonality\StationShp\SM_Pre_Bands/RainMask.txt');

yr = ncread('G:/PrecipGrid/AMPD.nc','year');
matsize = size(Dates);
medDate = ones(matsize(1),matsize(2))*nan;
meanDate = ones(matsize(1),matsize(2))*nan;
dateTrend = ones(matsize(1),matsize(2))*nan;

for i = 1:matsize(1)
    for j = 1:matsize(2)
        ts = reshape(Dates(i,j,:),matsize(3),1)/366*2*pi; % change dim of data 
        ts0 = reshape(Dates(i,j,:),matsize(3),1);
        ts0 = ts0(~isnan(ts0));
        goodYr = yr(~isnan(ts0)); 
        ts = ts(~isnan(ts));
        if length(ts) > 40
            medDate(i,j) = circ_median(ts,1)/2/pi*366; % 1 means along this dimension
            meanDate(i,j) = circ_mean(ts,[],1)/2/pi*366;  % do not forget to convert into DOY data
            dateTrend(i,j) = Theil_Sen_Regress(goodYr,ts0')*10;
        end
    end
end

% incase there are negative circular data
medDate(medDate<0) = medDate(medDate<0)+366;
meanDate(meanDate<0) = meanDate(meanDate<0)+366;


% set those outside EASM and BASIN region to NAN

medDate(~mask') = nan;
meanDate(~mask') = nan;
dateTrend(~mask') = nan;
figure;
imagesc(dateTrend);
colorbar
%% write annual maximum precipitation date data to output
%AMPD = netcdf.create('G:/2=Flood seasonality/StationShp/SM_Pre_Bands/SMdateTrendBand.nc','NC_WRITE');
AMPD = netcdf.create('G:/2=Flood seasonality\StationShp\SM_Pre_Bands/PRE_stats_EASM.nc','NC_WRITE');
dimlon = netcdf.defDim(AMPD,'lon',length(lon));
dimlat = netcdf.defDim(AMPD,'lat',length(lat));
var1 = netcdf.defVar(AMPD,'lon','NC_DOUBLE',dimlon);
var2 = netcdf.defVar(AMPD,'lat','NC_DOUBLE',dimlat);
var3 = netcdf.defVar(AMPD,'mean','NC_DOUBLE',[dimlon,dimlat]);
var4 = netcdf.defVar(AMPD,'median','NC_DOUBLE',[dimlon,dimlat]);
var5 = netcdf.defVar(AMPD,'trend','NC_DOUBLE',[dimlon,dimlat]);
netcdf.endDef(AMPD);
netcdf.putVar(AMPD,var1,lon);
netcdf.putVar(AMPD,var2,lat);
netcdf.putVar(AMPD,var3,meanDate);
netcdf.putVar(AMPD,var4,medDate);
netcdf.putVar(AMPD,var5,dateTrend);
netcdf.close(AMPD);
