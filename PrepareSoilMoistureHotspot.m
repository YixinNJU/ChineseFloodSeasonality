clear;
clc;

% ncdisp('./soilw.mon.mean.v2.nc');

% load data and read variables
% time span: 1948 - 2020
time = ncread('./soilw.mon.mean.v2.nc','time') + datenum(1800,1,1);
lon = ncread('./soilw.mon.mean.v2.nc','lon');
lat = ncread('./soilw.mon.mean.v2.nc','lat');
soilw = ncread('./soilw.mon.mean.v2.nc','soilw');

time = datetime(time,'ConvertFrom','datenum','Format','y-M-d'); 

%% SOIL MOISTURE MAGNITUDE 
% loop each grid to find the maximum Soil Moisture
msm = zeros(length(lon),length(lat), 2017-1959)*nan; % maximum soil moisture dates 

lon_ext = find(lon > 71 & lon < 150);
lat_ext = find(lat > 3 & lat < 60);
count = 1;
for t = 1960:2017
    if time.Year < 1960
        continue
    end
    disp(t);
    yrMat = soilw(lon_ext,lat_ext,(time.Year == t));
    yrMat(yrMat<=0) = nan;
    msm(lon_ext, lat_ext,count) = max(yrMat,[],3,'omitnan');
    count = count +1;
end

msm = msm(lon_ext, lat_ext,:); % maximum soil moisture dates
mmsm = mean(msm,3,'omitnan');
stdmsm = std(msm,[],3,'omitnan'); 
%values


%% Trend for grid data
% coherent trend for mainland China. Do not extract by basin.
 
[xlen, ylen, yrlen] = size(msm);
trend = zeros(xlen,ylen)*nan;
%md = trend;
yr = 1960:2017;

for i = 1:xlen
    for j = 1:ylen
        ts = msm(i,j,:);
        mm = mmsm(i,j);
        ts = ts/mm * 1000; % to /decades %:  /mm*10
        ts = ts(~isnan(ts));
        goodYr = yr(~isnan(ts));
        if isempty(ts) == 0
            trend(i,j,1) = Theil_Sen_Regress(goodYr',ts); 
            trend(i,j,2) = Mann_Kendall(ts,0.05) ; 
            %ts_deg = ts/366*360;
            %md(i,j) = cir_median(ts_deg,1);
        end
    end
end


%% TIMING ANALYSIS 
% loop for each grid to find the median DOY of the month with highest SM
msmd = zeros(length(lon),length(lat), 2017-1959)*nan; % maximum soil moisture dates 

count = 1;
%lon_ext = find(lon > 71 & lon < 150); % 137
%lat_ext = find(lat > 3 & lat < 60); % 55

for t = 1960:2017
    if time.Year < 1960
        continue
    end
    disp(t);
    
    for i = lon_ext(1):lon_ext(end)
        for j = lat_ext(1):lat_ext(end)
            yearData = soilw(i,j,(time.Year == t));
            yearData(yearData <= 0) = nan;
            
            yrDates = time(time.Year == t);
            
            if length(yrDates(isnan(yearData))) == 12
                %disp('No data this year');
                continue
            end
            
            maxMonth = find(yearData == max(yearData)); 
            
            % meadian date of that month
            %if maxMonth == 12 % December
            %   msmd(i,j,count) = day(yrDates(maxMonth),'dayofyear')+(366 - day(yrDates(maxMonth),'dayofyear')+1)/2-1;
            %    continue
            %end
            
            %msmd(i,j,count) = day(yrDates(maxMonth),'dayofyear')+(day(yrDates(maxMonth+1),'dayofyear') - day(yrDates(maxMonth),'dayofyear')+1)/2-1;
            msmd(i,j,count) = day(yrDates(maxMonth),'dayofyear')+14;
            
        end
    end
    count = count + 1;
end

msmd = msmd(lon_ext, lat_ext,:); % maximum soil moisture dates

%% write annual maximum SM date data to output
%lon = ncread('G:/CPCSoilMoisture\SoilMoisture/soilw.mon.mean.v2.nc','lon');
%lat = ncread('G:/CPCSoilMoisture\SoilMoisture/soilw.mon.mean.v2.nc','lat');
MSMD = netcdf.create('G:/CPCSoilMoisture\SoilMoisture/MSMD.nc','NC_WRITE');
dimlon = netcdf.defDim(MSMD,'lon',length(lon_ext));
dimlat = netcdf.defDim(MSMD,'lat',length(lat_ext));
dimtime = netcdf.defDim(MSMD,'year',length(1960:2017));
var1 = netcdf.defVar(MSMD,'lon','NC_DOUBLE',dimlon);
var2 = netcdf.defVar(MSMD,'lat','NC_DOUBLE',dimlat);
var3 = netcdf.defVar(MSMD,'year','NC_DOUBLE',dimtime);
var4 = netcdf.defVar(MSMD,'MSMD','NC_DOUBLE',[dimlon,dimlat,dimtime]);
netcdf.endDef(MSMD);
netcdf.putVar(MSMD,var1,lon(lon_ext));
netcdf.putVar(MSMD,var2,lat(lat_ext));
netcdf.putVar(MSMD,var3,1960:2017);
netcdf.putVar(MSMD,var4,msmd);
netcdf.close(MSMD);

 
 %% Hotspot analysis
 
% Hydro watersheds within hotspot
boxfn ='G:\\2=Flood seasonality\\StationShp\\HotspotShp6\\AllBoxes6_wgs.shp'; % NEwshd42.shp
box = shaperead(boxfn);

% gridded data to search
ncfn = 'G:\CPCSoilMoisture\SoilMoisture/MSMD.nc'; % G:/CPCSoilMoisture\SoilMoisture/MSMD.nc
ncdisp(ncfn);
yrs = ncread(ncfn,'year');
%lon_ext = ncread(ncfn,'lon');
%lons = ncread('./soilw.mon.mean.v2.nc','lon');
%lon = lons(lon_ext);
lon=ncread(ncfn,'lon');

%lat_ext = ncread(ncfn,'lat');
lats = ncread('G:\CPCSoilMoisture\SoilMoisture/soilw.mon.mean.v2.nc','lat');
%lat = lats(lat_ext);
lat=ncread(ncfn,'lat');

Data = ncread(ncfn,'MSMD'); % MSMD

stats = zeros(length(yrs),7)*nan; %  year; Precip_box1;Precip_box2;Precip_box3;Precip_box4;Precip_box5;Precip_box6;
stats(:,1) = yrs;
[X,Y] = meshgrid(lon,lat);
for i = 1:6
    rx = box(i).X(1:end-1);
    ry = box(i).Y(1:end-1);
    mask = inpolygon(X,Y,rx,ry);
    for k = 1:length(yrs)
        yrData = Data(:,:,k)';
        if sum(sum(isnan(mask)))/numel(mask) == 1
            disp('No data')
            continue
        end
        stats(k,i+1) = round(mean(yrData(mask),'omitnan'));
    end
end

 writematrix(stats, 'G:/2=Flood seasonality/StatResults/Hotspots/Hotspot_SM.csv');