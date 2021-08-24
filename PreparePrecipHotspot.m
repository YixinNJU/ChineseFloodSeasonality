%clc;
%clear;

ncdisp('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc');

%  read variables
time = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','time')/24+datenum(1961,1,1);
Data = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','pre');
lon = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','lon');
lat = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','lat');
disp('Data loaded.');

time = datetime(time,'ConvertFrom','datenum','Format','y-M-d'); % no need

%%  process data for 3-day and 7-day window
n = 7;
% calculate moving sum 7-day max precip for each grid
%newData = movsum(Data,n,3,'omitnan'); % along time dimension
newData = zeros(size(Data));
for i = 1:length(lon)
   for j = 1:length(lat)
        tmpData = Data(i,j,:);
       newData(i,j,:) = movsum(tmpData, n,'omitnan');
    end
end
disp('Moving n-day window for precipitation amout.')

%% Annual maximum precipitation magnitude
amp = zeros(length(lon), length(lat), 2015-1960)*nan;
yrs = 1961:2015;
for k = 1:length(yrs)
    tindx = find(time.Year == yrs(k));
    amp(:,:,k) = max(newData(:,:,tindx),[],3,'omitnan');
    nanid = find(sum(isnan(newData(:,:,tindx)),3) > 60); % find locations with more than two month missing value
    nanid(:,:,1:55) = k; % the k th year
    amp(nanid) = nan; % give nan to locations mentioned above
end
disp('Now we have each year max-n-d precipitation')
%% save amp data
lon = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','lon');
lat = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','lat');
AMP = netcdf.create('G:/PrecipGrid/AMP3d.nc','NC_WRITE');
dimlon = netcdf.defDim(AMP,'lon',length(lon));
dimlat = netcdf.defDim(AMP,'lat',length(lat));
dimtime = netcdf.defDim(AMP,'year',length(1961:2015));
var1 = netcdf.defVar(AMP,'lon','NC_DOUBLE',dimlon);
var2 = netcdf.defVar(AMP,'lat','NC_DOUBLE',dimlat);
var3 = netcdf.defVar(AMP,'year','NC_DOUBLE',dimtime);
var4 = netcdf.defVar(AMP,'AMP','NC_DOUBLE',[dimlon,dimlat,dimtime]);
netcdf.endDef(AMP);
netcdf.putVar(AMP,var1,lon);
netcdf.putVar(AMP,var2,lat);
netcdf.putVar(AMP,var3,1961:2015);
netcdf.putVar(AMP,var4,amp);

%% TIMING ANALYSIS
% squeeze to 55 year max date data
%find the date of max n day precip AMPD (annual maximum 7-day precip) 
%newData = data; % for 1-day time window
ampd = zeros(length(lon), length(lat), 2015-1960)*nan;
yrs = 1961:2015;

for i = 1:length(lon)
    %disp(i);   
    for j = 1:length(lat)
        for k = 1:length(yrs)
            tindx = find(time.Year == yrs(k));
            yrData = newData(i,j,tindx);
            yrData = reshape(yrData,length(yrData),1);
            
            if sum(isnan(yrData)) > 105 % this grid this year has less than 260 days  
                ampd(i,j,k) = nan; % this means there is no date for this year
                continue
            end 
            % if this year has precip
            tmpYear = time(tindx);
            tmpYear = reshape(tmpYear,length(tmpYear),1);
            pmax = find(yrData == max(yrData)); 
            tmax = round(median(day(tmpYear(pmax), 'dayofyear'),'omitnan'));
            ampd(i,j,k) = tmax;
        end
    end
end
    
disp('Max-n-day precipitation date each year each basin.');

%% write annual maximum precipitation date data to output
lon = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','lon');
lat = ncread('./CN05.1_Pre_1961_2015_daily_025x025/CN05.1_Pre_1961_2015_daily_025x025.nc','lat');
AMPD = netcdf.create('G:/PrecipGrid/AMPD3d.nc','NC_WRITE');
dimlon = netcdf.defDim(AMPD,'lon',length(lon));
dimlat = netcdf.defDim(AMPD,'lat',length(lat));
dimtime = netcdf.defDim(AMPD,'year',length(1961:2015));
var1 = netcdf.defVar(AMPD,'lon','NC_DOUBLE',dimlon);
var2 = netcdf.defVar(AMPD,'lat','NC_DOUBLE',dimlat);
var3 = netcdf.defVar(AMPD,'year','NC_DOUBLE',dimtime);
var4 = netcdf.defVar(AMPD,'AMPD','NC_DOUBLE',[dimlon,dimlat,dimtime]);
netcdf.endDef(AMPD);
netcdf.putVar(AMPD,var1,lon);
netcdf.putVar(AMPD,var2,lat);
netcdf.putVar(AMPD,var3,1961:2015);
netcdf.putVar(AMPD,var4,ampd);
netcdf.close(AMPD);

%% Hotspot analysis

% Hydro watersheds within hotspot
boxfn ='G:\\2=Flood seasonality\\StationShp\\HotspotShp6\\AllBoxes6_wgs.shp'; % NEwshd42.shp
box = shaperead(boxfn);

% gridded data to search
ncfn = 'G:/PrecipGrid/AMPD.nc'; % G:/CPCSoilMoisture\SoilMoisture/MSMD.nc
ncdisp(ncfn);
yrs = ncread(ncfn,'year');
lon = ncread(ncfn,'lon');
lat = ncread(ncfn,'lat');
Data = ncread(ncfn,'AMPD'); % MSMD

stats = zeros(length(yrs),7)*nan; %  year; Precip_box1;Precip_box2;Precip_box3;Precip_box4;Precip_box5;Precip_box6;
stats(:,1) = yrs;
[X,Y] = meshgrid(lon,lat);
for i = 1:6
    rx = box(i).X(1:end-1);
    ry = box(i).Y(1:end-1);
    mask = inpolygon(X,Y,rx,ry);
    for k = 1:length(yrs)
        yrData = Data(:,:,k)';
        %if k == 1 
        %    figure;
        %    imagesc(yrData);
        %    break
        %end
        if sum(sum(isnan(mask)))/numel(mask) == 1
            disp('No data')
            continue
        end
        stats(k,i+1) = round(mean(yrData(mask),'omitnan'));
    end
end
%
writematrix(stats, 'G:/2=Flood seasonality/StatResults/Hotspots/Hotspot_PRE.csv');