# ChineseFloodSeasonality
Descriptions of the scripts used in Yang et al. (2021): Climate Leads to Reversed Latitudinal Changes in Chinese Flood Peak Timing.
1. *FloodDateSeasonalityTrend.R*: Calculate mean **seasonality (i.e., occurrence date)** and **timing trend** based on circular statistics and Sen's slope estimator. 
2. *FloodTrendInterpolation.R*: Interpolate flood timing trend using **Ordinary Kriging** automatically.
3. *CutZonalBelts.m*: Calculate **statistics** (mean seasonality and timing trend) for precipitation and soil moisture grids masked by **EASM region**.
4. *PrepareFloodHotspot.R*: Mask **flood timing** time series using extents of 6 hotspots.
5. *PreparePrecipHotspot.m*: Calculate timing of **annual maximum 7-day precipitation** during 1961-2015 and clip time series by location of 6 hotspots.
6. *PrepareSoilMoistureHotspot.m*: Calculate timing of **annual maximum monthly soil moisture** during 1960-2017 and clip time series by location of 6 hotspots.
7. *HotspotTrend.R*: Calculate **trends** and do **circular correlation** based on time series obtained from 4-6.
8. *FloodRainDBSCAN.R*: Use **DBSCAN** to do **spatial cluster** based on flood and rainfall timing data.
