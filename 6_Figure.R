ReqPackages = c("terra", "tidyterra", "gdalUtils", "sf", "rts", "stringi", "lubridate", "mgcv", "data.table", 
                "ggplot2", "dplyr", "tidyr", "geodata", "viridis", "ggmap", "maptools", "raster", "ggpubr")
#install.packages(ReqPackages)
invisible(lapply(ReqPackages, library, character.only = T, quietly = T))

####################################################################
# Figure 1: showing all the pixels went into the analysis (NCEP + MODIS)
evi = rast("./data/Rasters/6_MOD13A2_EVI_interpolated_rescaled.tif") #2000-01-01 to 2020-12-31
aqua = rast("./data/Rasters/5_MYD21A2_GreatLakesTemp_interpolated_rescaled.tif") #2003-01-01 to 2020-12-31
prec = rast("./data/Rasters/3_ClimateDatabase_2000to2020_Prec_DaySum_shum.2m.tif")
winddir = rast("./data/Rasters/4_ClimateDatabase_2000to2020_Wind_1000p_Daytime_WindDirection.tif")
temp = rast("./data/Rasters/2_ClimateDatabase_2000to2020_Temp_AvDay_air.2m.tif")
maskwind = rast("./data/Masks/4_ClimateDatabase_2000to2020_Wind_1000p_Daytime_WindDirection_OceanMask.tif")
maskprec = rast("./data/Masks/3_ClimateDatabase_2000to2020_Prec_DaySum_shum.2m_OceanMask.tif")
masktemp = rast("./data/Masks/2_ClimateDatabase_2000to2020_Temp_AvDay_air.2m_OceanMask.tif")

radar = read.csv("./data/maps/radars_great lake.csv", header = TRUE, sep = ",")
radar_points = st_as_sf(radar, coords = c("longitude", "latitude"), crs = canlam)
radar_points = vect(radar, geom=c("longitude", "latitude"), crs="WGS84", keepgeom=FALSE)


# Print information about the SpatVector object
print(points)

# Get the US + Canada + Great Lakes basemap

# lake = st_read("./data/maps/layers/greatlakes.shp")
# canada = readRDS(("./data/maps/gadm/gadm41_CAN_1_pk_low.rds"))
# canada = st_as_sf(canada)
# us = readRDS(("./data/maps/gadm/gadm41_USA_1_pk_low.rds"))
# us = st_as_sf(us)

us <- gadm(country = "USA", level = 1, resolution = 2,
           path = "./data/maps/")
canada <- gadm(country = "CAN", level = 1, resolution = 2,
               path = "./data/maps")
lake = terra::vect("./data/maps/layers/greatlakes.shp")
CanUS <- rbind(us, canada, lake)

plot(CanUS, xlim = c(-180, -50))

###############
canlam <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

## project our vector data:
CanUS.lcc <- project(CanUS, canlam)

plot(evi[[1277]])
plot(aqua[[182]])
plot(prec[[1277]])
plot(winddir[[1227]])

evi.lcc <- project(evi[[1277]], canlam) #2003-07-01
aqua.lcc <- project(aqua[[182]], canlam) #2003-07-01
prec.lcc <- project(prec[[1277]], canlam) #2003-07-01
temp.lcc <- project(temp[[1277]], canlam) #2003-07-01
winddir.lcc <- project(winddir[[1277]], canlam) #2003-07-01
maskwind = project(maskwind, canlam)
maskprec = project(maskprec, canlam)
masktemp = project(masktemp, canlam)
radar_points.lcc = project(radar_points, canlam)

prec = mask(prec.lcc, maskprec)
wind = mask(winddir.lcc, maskwind)
temp = mask(temp.lcc, masktemp)

####################

pdf("/Volumes/T7/Code_changed_2020_HaestB/plot/1_aqua_data.pdf", width=10, height=8, pointsize = 20)
plot(aqua.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000))
plot(CanUS.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000), border = "grey", add = TRUE)
plot(radar_points.lcc, add=TRUE)
dev.off()

pdf("/Volumes/T7/Code_changed_2020_HaestB/plot/1_evi_data.pdf", width=10, height=8, pointsize = 20)
plot(evi.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000))
plot(CanUS.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000), border = "grey", add = TRUE)
plot(radar_points.lcc, add=TRUE)
dev.off()

pdf("/Volumes/T7/Code_changed_2020_HaestB/plot/1_prec_data.pdf", width=10, height=8, pointsize = 20)
plot(prec, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000))
plot(CanUS.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000), border = "grey", add = TRUE)
plot(radar_points.lcc, add=TRUE)
dev.off()

pdf("/Volumes/T7/Code_changed_2020_HaestB/plot/1_wind_data.pdf", width=10, height=8, pointsize = 20)
plot(wind, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000))
plot(CanUS.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000), border = "grey", add = TRUE)
plot(radar_points.lcc, add=TRUE)
dev.off()

pdf("/Volumes/T7/Code_changed_2020_HaestB/plot/1_temp_data.pdf", width=10, height=8, pointsize = 20)
plot(temp, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000))
plot(CanUS.lcc, xlim = c(-3000000, 3000000), ylim = c(-2000000, 3000000), border = "grey", add = TRUE)
plot(radar_points.lcc, add=TRUE)
dev.off()

pdf("/Volumes/T7/Code_changed_2020_HaestB/plot/1_map.pdf", width=10, height=8, pointsize = 20)
plot(temp, xlim = c(100000, 1600000), ylim = c(-1200000, 200000))
plot(CanUS.lcc, border = "white", xlim = c(100000, 1600000), ylim = c(-1000000, 0), add = TRUE)
plot(radar_points.lcc, add=TRUE)
text(radar_points.lcc, labels = c("KAPX", "KBUF", "KCLE", "KDLH", "KDTX", "KGRB", "KGRR", "KIWX", "KLOT", "KMQT", "KMKX", "KTYX"), pos = 1, offset = 0.7)
dev.off()

####################################################################
####################################################################
####################################################################

# Figure 2: Final results in the top signals 

us <- gadm(country = "USA", level = 1, resolution = 2,
           path = "./data/maps/")
canada <- gadm(country = "CAN", level = 1, resolution = 2,
               path = "./data/maps")
lake = terra::vect("./data/maps/layers/greatlakes.shp")
CanUS <- rbind(us, canada, lake)
CanUS = st_as_sf(CanUS)

radar = read.csv("./data/maps/radars_great lake.csv", header = TRUE, sep = ",")
radar_points = st_as_sf(radar, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")

final_infl = read.csv("./results/5_FinalSignals/top5.csv") %>% 
  mutate(median =  as.Date(258-(WOpen+WClose)/2, origin="0999-01-01"))
final_infl_points = st_as_sf(final_infl, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")

puma = terra::rast("./data/maps/puma/purmar_abundance_seasonal_breeding_mean_2021.tif") 
  puma_latlon = project(puma, "EPSG:4326") 
  puma_crop = crop(puma_latlon, e <- ext(-140, -50, 20, 70)) 
  puma_df = as.data.frame(puma_crop, xy=T)
  puma_df$breeding = na_if(puma_df$breeding, 0)
  puma_df = puma_df %>% na.omit()

trsw = terra::rast("./data/maps/trsw/treswa_abundance_seasonal_breeding_mean_2021.tif") 
  trsw_latlon = project(trsw, "EPSG:4326") 
  trsw_crop = crop(trsw_latlon, e <- ext(-140, -50, 20, 70)) 
  trsw_df = as.data.frame(trsw_crop, xy=T)
  trsw_df$breeding = na_if(trsw_df$breeding, 0)
  trsw_df = trsw_df %>% na.omit()
  
# trsw = terra::rast("./data/maps/trsw/treswa_abundance_seasonal_breeding_mean_2021.tif") %>% 
#   project("EPSG:4326") %>% 
#   crop(ext(-140, -50, 20, 70)) %>% 
#   as.data.frame(xy=T) %>% 
#   mutate(breeding = na_if(breeding, 0)) %>% 
#   na.omit() 

overlap = merge(puma_df, trsw_df, by = c("x", "y"), all = FALSE)

# final_infl = read.csv("./results/4_MultiVar/7_C100_Confirmation/yday_50_MADJDAYSWS_RemainingSignals.csv") %>% 
#   filter(is.na(RemovedPriorToVarImp)) %>% 
#   separate(SourceCoords, c("latitude", "longitude"), ",") %>% 
#   mutate(longitude = as.numeric(longitude)-360,
#          latitude = as.numeric(latitude))

# final_infl = st_as_sf(final_infl, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")

pdf(file = "./plot/2_Top5_range_updated.pdf",   # The directory you want to save the file in
    width = 6.98, # The width of the plot in inches
    height = 5.82) # The height of the plot in inches

ggplot() +
  geom_raster(data = trsw_df, aes(x=x, y=y, fill=breeding), fill="#26ABD9", alpha=0.5)+
  geom_raster(data = puma_df, aes(x=x, y=y, fill=breeding), fill="#4914DD", alpha=0.5)+
  geom_raster(data = overlap, aes(x=x, y=y, fill=breeding), fill="#F8766D", alpha=0.5)+
  geom_sf(data = CanUS, fill = NA, color = "grey")+
  geom_sf(data = final_infl_points, aes(color = ClimateVar, size = Mean_Rank)) +
  scale_size(range = c(4, 10))+
  #geom_sf(data = radar_points, color = "black", size=2, shape=13)+
  coord_sf(xlim = c(-140, -50), ylim = c(20, 70), expand = FALSE)+
  scale_colour_viridis_d()+
  labs(fill = "Data Value")+
  theme_classic()

dev.off()

# time-window for the final signals
daily.df = read.csv("./data/Chap1_pheno_data/df_0fillted.csv") %>% 
  aggregate(bird_count_polar ~ date,sum) %>% 
  mutate(fake_date = as.Date(date, format="%Y-%m-%d"))
year(daily.df$fake_date) = 999

b = ggplot() +
  geom_violin(daily.df, mapping=aes(x=fake_date, y=bird_count_polar))+
  scale_x_date(limits = as.Date(c("0999-01-01", "0999-10-01")), date_breaks = "1 month", date_labels = "%b")+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  labs(y="Daily number of swallows")+
  theme(axis.title.x = element_blank())

c = ggplot() +
  geom_point(final_infl, mapping=aes(x=median, y=ClimateVar)) +
  geom_errorbarh(final_infl, mapping=aes(xmin = as.Date(258-WOpen, origin="0999-01-01"), 
                                       xmax = as.Date(258-WClose, origin="0999-01-01"),
                                       y = ClimateVar,
                                       height = .2))+
  scale_x_date(limits = as.Date(c("0999-01-01", "0999-10-01")), date_breaks = "1 month", date_labels = "%b")+
  labs(y="Climate variable")+
  theme(axis.title.x = element_blank())


pdf(file = "./plot/2_Time_window.pdf",   # The directory you want to save the file in
    width = 6.98, # The width of the plot in inches
    height = 5.82) # The height of the plot in inches

ggarrange(b, c, nrow = 2, labels = c("B", "C"), align = "v")
dev.off()

ggarrange(a,                                                 # First row with map
          ggarrange(b, c, nrow = 2, labels = c("B", "C"), align = "v"), # Second row with time-window plots
          ncol = 2, 
          labels = "A"                                        # Labels of the map
) 


####################################################################
all_influence = st_read("./results/5_FinalSignals/2_Shapes/1_All/yday_50_MADJDAYSWS_FinalOverview.shp")
df_of_vars=as.data.frame(all_influence)
df_of_vars$XLoc= df_of_vars$XLoc-360 # all_influence shapefile is on a scale of 0-360, its longitude needs to be converted into 180 to 180
df_of_vars= st_as_sf(df_of_vars, coords=c("XLoc","YLoc"), crs= "+proj=longlat +datum=WGS84")
df_of_vars=as.data.frame(all_influence)
df_of_vars$XLoc= df_of_vars$XLoc-360
df_of_vars= st_as_sf(df_of_vars, coords=c("XLoc","YLoc"), crs= "+proj=longlat +datum=WGS84")

three_infl = read.csv("./results/5_FinalSignals/1_OverviewFiles/yday_50_MADJDAYSWS_FinalOverview.csv") %>% 
  filter(Mean_Rank > 9) %>% # the top 3
  separate(SourceCoords, c("latitude", "longitude"), ",") %>% 
  mutate(longitude = as.numeric(longitude)-360,
         latitude = as.numeric(latitude))
three_infl = st_as_sf(three_infl, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")


tmin = rast("./results/3_SpaWin_Adj/MADJDAYSWS/BestModelSumm/air.min_yday50_Climwin/AICc/yday_50.tif")
df_of_vars = cbind(as.data.frame(tmin), as.data.frame(crds(tmin)))
df_of_vars$x <- df_of_vars$x - 360
df_of_vars <- df_of_vars[, c(2,3,1)]
tmin = rast(df_of_vars, type="xyz", crs="+proj=longlat +datum=WGS84", digits=6, extent=NULL)
plot(tmin)

ggplot() +
  geom_sf(data = CanUS, fill = NA, color = "grey")+
  tidyterra::geom_spatraster(data = tmin)+
  scale_colour_viridis_c() +
  labs(fill = "AICc")+
  theme_classic()

ggplot() +
  geom_sf(data = CanUS, fill = NA, color = "grey")+
  geom_sf(data = df_of_vars, aes(color = ClimtVr, size = R_squrs)) +
  geom_sf(data = radar_points, color = "black", size=2, shape=18)+
  coord_sf(xlim = c(-180, -50), ylim = c(20, 70), expand = FALSE)+
  scale_colour_viridis_d()+
  labs(title = "Top 14 Influential Environmental Virables",
       fill = "Data Value")+
  theme_classic()




####################################################################
####################################################################
####################################################################

#Figure 3: Overview of the effects of the identified weather signals on phenology (chain rules)

pheno_top5_df = read.csv("./results/5_FinalSignals/pheno_top5_df.csv") %>% 
  gather(key="variable", value="value", 4:8) 

#dClimate/dTime
ggplot(pheno_top5_df, aes(x=Year, y=value)) +
  geom_point()+
  geom_smooth(method="lm", formula = y~x, color="black")+
  facet_wrap(~variable, ncol = 3, scales = "free_y")

for (i in unique(pheno_top5_df$variable)){
  model = lm(value ~ Year,
             data = pheno_top5_df %>% filter(variable==i))
  print(summary(model))
}

#dPheno/dClimate
ggplot(pheno_top5_df, aes(x=value, y=MADJDAYSWS)) +
  geom_point()+
  geom_smooth(method="lm", formula = y~x, color="black")+
  facet_wrap(~variable, ncol = 3, scales = "free_x")

for (i in unique(pheno_top3_df$variable)){
  model = lm(MADJDAYSWS ~ value,
           data = pheno_top3_df %>% filter(variable==i))
  print(summary(model))
}

#dPheno/dTime
ggplot(pheno_top3_df, aes(x=Year, y=MADJDAYSWS)) +
  geom_point()+
  geom_smooth(method="lm", formula = y~x, color="black")+
  facet_wrap(~variable, ncol = 3)


####################################################################
####################################################################
####################################################################

#Figure: Time series of a single pixel in an example
library(rts)
temp = rast("./data/Rasters/2_ClimateDatabase_2000to2020_Temp_AvDay_air.2m.tif")
d <- seq(as.Date('2000-01-01'), as.Date('2020-12-31'), by = "1 days") # corresponding dates to 4 rasters
d <- as.Date(d)
timeseries <- rts(temp,d) # we creater a SpatRasterTS (a raster time series)
timeseries = timeseries[[1:365]]
plot(timeseries[[200]])
c <- cellFromXY(timeseries, matrix(c(280, 50),nrow=1))
plot(timeseries[c])

####################################################################
####################################################################
####################################################################
#### Input the tif files from k-fold results ####

inputDir  = "./results/8_Cross_Validation/FinalVariable_allvariables/shum.2m/"
rast = raster(paste0(inputDir,"/AICc/yday_50.tif"))
## S4 method for signature 'RasterLayer'
KML(rast, paste0(inputDir,"/AICc/shum.2m.kml"), col=rev(terrain.colors(255)), 
    colNA=NA, maxpixels=100000, blur=1, zip='', overwrite=FALSE)
