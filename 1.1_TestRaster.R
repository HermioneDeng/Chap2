library(raster)

wind = brick("./data/Rasters/4_ClimateDatabase_2000to2020_Wind_AvDay_uwnd.tif")
plot(wind[[3]])

wind = brick("./data/Rasters/4_ClimateDatabase_2000to2020_Wind_AvDay_vwnd.tif")
plot(wind[[3]])

wind = brick("./data/Rasters/4_ClimateDatabase_2000to2020_Wind_AvDay_WindDirection.tif")
plot(wind[[3]])

wind_rds = readRDS("./data/4_ClimateDatabase_2000to2020_Wind_AvDay.rds")
str(wind_rds)

#########
precip = brick("./data/Rasters/3_ClimateDatabase_2000to2020_Prec_DaySum_shum.2m.tif")
plot(precip[[1]])

############
temp = brick("./data/Rasters/2_ClimateDatabase_2000to2020_Temp_AvDay_air.2m.tif")
plot(temp[[100]])

temp = brick("./data/Rasters/2_ClimateDatabase_2000to2020_Temp_AvDay_tmin.2m.tif")
plot(temp[[1]])

writeRaster(temp,'test.tif', options=c('TFW=YES'))

###############
evi = brick("./data/Rasters/6_MOD13A2.tif")
plot(evi[[1]])

#####
#check the masks
#####

mask = raster("./data/Masks/4_ClimateDatabase_2000to2020_Wind_AvDay_WindDirection_OceanMask.tif")
plot(mask)


library(rgdal)
shp <- readOGR(ShapeMaskPath)
plot(shp)
