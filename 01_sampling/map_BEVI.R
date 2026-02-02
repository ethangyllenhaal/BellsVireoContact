# code for making sampling maps used in manuscript

# load libraries
library(tidyverse)
library(sf)

### River code based on https://milospopovic.net/map-rivers-with-sf-and-ggplot2-in-r/

# set working dir, read sampling files
setwd("C:/Documents/Projects/BEVI/sampling")
samples <- read.csv("bells_sampling_coords.csv")

# download USA map
us <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

# set sf option
sf::sf_use_s2(F)

# downloaded from https://www.hydrosheds.org/products/hydrorivers
# in working dir
raw_rivers <- sf::st_read("HydroRIVERS_v10_na_shp")

# set projection, subset rivers to focal area
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"
river_subset <- st_sfc(
  st_polygon(list(cbind(
    c(-113, -102, -102, -113, -113),
    c(29, 29, 37, 37, 29)
  ))), crs = crsLONGLAT
)
subset_rivers <- st_crop(raw_rivers, river_subset)

# subset by size for different plotting parameters
big_rivers <- subset_rivers %>%
  filter(ORD_FLOW<=5)  |> st_cast("MULTILINESTRING")
med_rivers <- subset_rivers %>%
  filter(ORD_FLOW==6) |> st_cast("MULTILINESTRING")
small_rivers <- subset_rivers %>%
  filter(ORD_FLOW==7) |> st_cast("MULTILINESTRING")

# river basins (not used in final version)
raw_basins <- sf::st_read("hybas_lake_na_lev05_v1c")
na_basins <- raw_basins |> 
  sf::st_cast("MULTILINESTRING")

# set bounding box
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"
box <- st_sfc(
  st_polygon(list(cbind(
    c(-112, -103, -103, -112, -112),
    c(30, 30, 36, 36, 30)
  ))), crs = crsLONGLAT
)
proj <- st_transform(box, crs=4087)
bounds <- st_bbox(proj)

# import coordinates from sample dataframe
points <- data.frame(samples$Lon, samples$Lat)
coords <- st_as_sf(points, coords = c(1:2), crs = crsLONGLAT)
coords_transform <- st_transform(coords, crs=4087) |> st_jitter(amount=10000)
coords_nojitter <- st_transform(coords, crs=4087)
coords_info <- cbind(coords_transform, data.frame(samples))

# Map with points and rivers for Figure 1
ggplot() +
  geom_sf(data = na_basins, color = "gray85") + # basins (removed for most figures)
  geom_sf(data = med_rivers, color = "blue", lwd=1) + # medium rivers
  geom_sf(data = big_rivers, color = "darkblue", lwd=1.1) +  # big rivers
  geom_sf(data = us, fill = NA) + # US background
  geom_sf(data=coords_info, size=5, aes(geometry=geometry, color=mtDNA)) + # coordinates
  scale_color_manual(values=c("#32DEDE", "#c83200")) +
  coord_sf(crs=4087, xlim=c(bounds["xmin"], bounds["xmax"]), ylim = c(bounds["ymin"], bounds["ymax"])) + # bounding box
  theme_void()

# Map with small rivers for Figure 6
ggplot() +
  geom_sf(data = na_basins, color = "gray85") +
  geom_sf(data = small_rivers, color = "lightblue", lwd=0.8) +
  geom_sf(data = med_rivers, color = "blue", lwd=1) +
  geom_sf(data = big_rivers, color = "darkblue", lwd=1.1) +  
  geom_sf(data = us, fill = NA) +
  coord_sf(crs=4087, xlim=c(bounds["xmin"], bounds["xmax"]), ylim = c(bounds["ymin"], bounds["ymax"])) +
  theme_void()

# Figure 1 inset map
ggplot() +
  geom_sf(data = box, fill=NA, color="black", lwd=1) +
  geom_sf(data = us, fill = NA) +
  theme_void()

# Figure 2 (PCA) map, same idea as above but zoomed in
box_pca <- st_sfc(
  st_polygon(list(cbind(
    c(-111, -103.5, -103.5, -111, -111),
    c(30, 30, 36, 36, 30)
  ))), crs = crsLONGLAT
)
proj_pca <- st_transform(box_pca, crs=4087)
bounds_pca <- st_bbox(proj_pca)

ggplot() +
  geom_sf(data = med_rivers, color = "blue", lwd=1) +
  geom_sf(data = big_rivers, color = "darkblue", lwd=1.1) +  
  geom_sf(data = us, fill = NA) +
  geom_sf(data=coords_nojitter, size=5, aes(geometry=geometry)) +
  coord_sf(crs=4087, xlim=c(bounds_pca["xmin"], bounds_pca["xmax"]), ylim = c(bounds_pca["ymin"], bounds_pca["ymax"])) +
  theme_void()



  
