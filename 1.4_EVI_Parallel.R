ReqPackages = c("lattice", "terra", "gdalUtils", "fs", "rts", "stringi", "lubridate", "mgcv", "data.table", "ggplot2", "dplyr", "gridExtra", "tidyterra")
#install.packages(ReqPackages)
invisible(lapply(ReqPackages, library, character.only = T, quietly = T))


## Read in the files and create raster stack

input_EVI = dir_ls("/Volumes/T7/MOD13A2_EVI_tif_merged") #16-day resolution
output = "./data"

# # SpatRaster
# for (i in 1:471){
#   r <- rast(input_EVI[i])
#   print(origin(r))
# }
# 
# # not all rasters' origins are aligned
# terraOptions(tolerance=.0000001)
r <- rast(input_EVI)
writeRaster(r, filename = "./data/Rasters/6_MOD13A2_EVI.tif", overwrite = TRUE)

# files <- list.files(path = "/Volumes/T7/MOD13A2_EVI_tif_merged", pattern = "*.tif", full.names = TRUE)
# rasterList <- lapply(files, function(x) rast(x))
# #rasterBrick <- terra::merge(rasterList)
# rasterBrick <- rasterList[[1]]
# for (i in 2:length(rasterList)) {
#   rasterBrick <- terra::merge(rasterBrick, rasterList[[i]], overwrite = TRUE)
# }

r_Rescaled = aggregate(r, fact=16, fun=mean) #1/4 of the orginal size - 4 x 4 km
writeRaster(r_Rescaled, filename = "./data/Rasters/6_MOD13A2_EVI_rescaled.tif", overwrite = TRUE)
r_Rescaled = rast("./data/Rasters/6_MOD13A2_EVI_rescaled.tif")
plot(r_Rescaled[[101]])

## Create raster stack time series and extract values for each pixel
year = stri_sub(input_EVI, -11, -8)
yday = stri_sub(input_EVI, -7, -5)
yday <- as.numeric(sub("^0+", "", yday)) #Delete Leading Zeros Using sub() Function
d = as.Date(yday-1, origin = paste0(year, "-01-01"))

# rt <- rts(r_Rescaled,d) # creating a RasterStackTS object
# plot(rt[[c(200)]])
# 
# # Extract time series from a pixel
# # you may use cellFromXY to get the cell number of a location given its coordinate
# c <- cellFromXY(rt, matrix(c(-80, 41), nrow=1))
# plot(rt[c])
# extract(rt, c)
# This is how to get the values from a single pixel
cellFromXY(r_Rescaled, matrix(c(-80, 41), nrow=1))

## Analyzing seasonal time series with GAM
d_all = seq(as.Date('2000-01-01'), as.Date('2020-12-31'), by = "1 days")

Database_Arr = array(NA_real_,
                     dim = c(dim(r_Rescaled)[1],
                             dim(r_Rescaled)[2],
                             length(d_all)))

dev_df = data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("pixel", "dev"))))
start = Sys.time()
for (i in 1:dim(r_Rescaled)[1]){
  for (j in 1:dim(r_Rescaled)[2]){
    
    c = cellFromRowCol(r_Rescaled, i, j)
    cat('Index number:', c, '\n')
    
    df = data.frame(doy = yday(d), year = year(d), date = d, EVI = t(extract(r_Rescaled, c)))
    df = na.omit(df)
    
    # ggplot(df, aes(date, EVI))+
    #   geom_line()
    
    if (nrow(df)<21 | mean(table(df$year))<12){
      Database_Arr[i,j,] = NA
    }
    else{
      # gam.model <- gam(EVI ~ s(year, bs = "tp", k = length(unique(df$year))) +
      #                    s(doy, bs = 'cc', k = length(unique(df$doy))),
      #                  data = df,
      #                  family = gaussian)
      
      gam.model <- gam(EVI ~ te(year, doy, 
                                 k = mean(table(df$year))-1,
                                 bs = c("ps", "cr")),
                       data = df,
                       family = gaussian)
      
      # gam.model <- gam(EVI ~ t2(doy, year,
      #                        k = c(20, 17),
      #                        bs = c("cr", "ps"),
      #                        full = TRUE),
      #              data = df,
      #              family = gaussian)
      
      dev_df[nrow(dev_df)+1, ] = data.frame(pixel = c, dev = summary(gam.model)$dev.expl)

      # datas <- rbindlist(list(subset(df, select=c(EVI, date)),
      #                         data.table(EVI = gam.model$fitted.values,
      #                                    date = df$date)))
      # datas$type = c(rep("Real", nrow(df)), rep("Fitted", nrow(df)))
      # 
      # ggplot(data = datas, aes(date, EVI, group = type, colour = type)) +
      #   geom_line(linewidth = 1) +
      #   theme_bw() +
      #   labs(x = "Time", y = "EVI")


      df_daily = data.frame(date = d_all) %>%
        mutate(year=year(date),
               doy=yday(date))
      df_daily = df_daily %>% mutate(EVI = predict.gam(gam.model, newdata=df_daily))
      
      # fill in the temporal gaps (keep the original measurements, only add estimates for the days without measurements)
      df_filled = rbind(df, df_daily[!df_daily$date %in% df$date,])
      df_filled = df_filled[order(as.Date(df_filled$date, format="%Y-%m-%d")),]
      
      # ggplot() +
      #   geom_line(data = datas %>% filter(type=="Real"), aes(date, EVI), size = 0.5, color="red", alpha=1) +
      #   geom_line(data=df_filled, aes(date, EVI), alpha=0.5)+
      #   theme_bw() +
      #   labs(x = "Time", y = "EVI")
      
      Database_Arr[i,j,] = df_filled$EVI
      
    }
  }
}
end = Sys.time()
end-start



## Save the data into 3D array and a raster brick


# Convert the array to a rasterbrick object
# ===============================================================
Outputbrick = rast(Database_Arr)
ext(Outputbrick) = ext(r_Rescaled)
crs(Outputbrick) = crs(r_Rescaled)
# creating a RasterStackTS object:
Outputbrick_rts <- rts(Outputbrick, d_all)

# Agregate again by a factor of 4 (64 x 64 km)
# ===============================================================
Outputbrick_Rescaled = aggregate(Outputbrick, fact=4, fun=mean)

# Write the raster to a file
# ===============================================================
writeRaster(Outputbrick, "./data/Rasters/6_MOD13A2_EVI_interpolated.tif", overwrite=TRUE)
writeRaster(Outputbrick_Rescaled, "./data/Rasters/6_MOD13A2_EVI_interpolated_rescaled.tif", overwrite=TRUE)
write.csv(dev_df, "./EVI_deviance_explained.csv")

# Convert the current database to an array of the desired size
# ===============================================================
Outputbrick_Rescaled = rast("./data/Rasters/6_MOD13A2_EVI_interpolated_rescaled.tif")
Outputbrick_Rescaled_array = as.array(Outputbrick_Rescaled)

Database_Arr = array(NA_real_,
                     dim = c(dim(Outputbrick_Rescaled)[1]*dim(Outputbrick_Rescaled)[2],
                             dim(Outputbrick_Rescaled)[3], 1))

# j represents each day of the data (7671 days in total)
for (j in 1:dim(Outputbrick_Rescaled)[3]){
  Database_Arr[,j,1] = Outputbrick_Rescaled_array[,,j]
}

dimnames(Database_Arr) = list(NULL, as.character(d_all), "EVI")
saveRDS(Database_Arr, file = "./data/6_MOD13A2_EVI_interpolated_rescaled.rds")
Database_Clim = readRDS(file = "./data/6_MOD13A2_EVI_interpolated_rescaled.rds")  

## Checking the results
# ===============================================================
pdf("EVI_pixel_examine.pdf", width = 12, height = 7)

for (i in 1:dim(Outputbrick)[1]){
  for (j in 1:dim(Outputbrick)[2]){
    
c = cellFromRowCol(Outputbrick, i, j)
    
out_df = data.frame(date = d_all, value = t(extract(Outputbrick, c)))
og_df = data.frame(date = d, value = t(extract(r_Rescaled, c)))

if (is.na(out_df$value[100])==TRUE){
  next
}

else{
  cat('Index number:', c, '\n')
  
  print(
  ggplot()+
  geom_line(og_df, mapping=aes(x=date, y=value), color="red", alpha=0.8)+
  geom_line(out_df, mapping=aes(x=date, y=value), color="black", alpha=0.6)+
  labs(title=paste0("row=", i, " ", "col=", j, " ", "index=", c))+
  theme_classic()
)
}

  }
}

dev.off()

plot(Outputbrick_Rescaled_rts[[1:10]])

for(c in 1:56544) {       # Printing ggplot within for-loop
  print(
    plot(Outputbrick_Rescaled_rts[c])
  )
  Sys.sleep(2)
}


for (c in 1:56544) {
  
  skip_to_next <- FALSE
  
  tryCatch(print(plot(Outputbrick_Rescaled_rts[c])), 
           Sys.sleep(2),
           error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}

# ===============================================================
pdf("EVI_interpolated_rescaled.pdf")
for (i in seq(1,7671,3)) {
    p <- list()
    for(j in 1:3){
      p[[j]] <-  ggplot() +
        geom_spatraster(data = Outputbrick_Rescaled_rts[[i+j-1]]) +
        scale_fill_whitebox_c(
          palette = "muted",
          labels = scales::label_number(suffix = ""),
          limits=c(-2000, 10000)) +
        labs(fill = "EVI", title= paste0(Outputbrick_Rescaled_rts[[i+j-1]]@ptr$names, "  index = ", i+j-1))+
        theme_classic()
    }
    do.call(grid.arrange,p)
}
dev.off()


pdf("EVI_pixel_resclaed_examine.pdf", width = 12, height = 7)

for (i in 1:dim(Outputbrick_Rescaled)[1]){
  for (j in 1:dim(Outputbrick_Rescaled)[2]){
    
    c = cellFromRowCol(Outputbrick_Rescaled, i, j)
    
    out_df = data.frame(date = d_all, value = t(extract(Outputbrick_Rescaled, c)))
    
    if (is.na(out_df$value[100])==TRUE){
      next
    }
    
    else{
      cat('Index number:', c, '\n')
      
      print(
        ggplot()+
          geom_line(out_df, mapping=aes(x=date, y=value), color="black", alpha=0.6)+
          labs(title=paste0("row=", i, " ", "col=", j, " ", "index=", c))+
          theme_classic()
      )
    }
    
  }
}
dev.off()
