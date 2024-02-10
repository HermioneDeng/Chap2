library(sf)
library(raster)
library(rgdal)
library(ggplot2)

map = readOGR("./data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
map = readOGR("./data/North_American_Atlas_Bathymetry/bathy_p.shp")
map = st_read("./data/shapefile/bathy_p.shp")
map = st_read("./data/shapefile/ne_10m_admin_0_countries.shp")
map = map$geometry

ncep = raster("./test.tif")
ncep <- raster::rotate(ncep)


st_crs(map)
identical(st_crs(map), crs(ncep))
map <- spTransform(x = map, CRSobj = crs(ncep))
map <- st_transform(x = map, CRSobj = crs(ncep))

ncep_df = as.data.frame(ncep, xy=TRUE)

plot(ncep)
plot(map, add=TRUE)

ggplot() + 
  geom_sf(data = map, size = 0.5, color = "black", fill = "cyan1") + 
  #geom_raster(data = ncep_df , aes(x = x, y = y, fill=layer)) + 
  coord_sf(crs = st_crs(4326))

masked <- mask(x = ncep, mask = map)
masked = st_crop(x=ncep_df, y=map)

plot(masked)

  
