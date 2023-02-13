# Determine lake sites

# Set up work space ----
# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(sf)
library(dbhydroR)
library(sf)
library(janitor)
# library(raster)
# library(exactextractr)
# library(colormap)
library(lubridate)
# library(quantreg)
library(beepr)
library(wesanderson)

# scales::show_col(
#   wes_palette(
#     name = 'Royal1',
#     n = 4,
#     type = c('continuous')
#   ) #%>% as.vector()
# )


# Set working directory
setwd('/Users/natalie/Documents/select_data/') 

# Read in data ----

load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")

## Set coordinate system
epsg <- crs(stack_daily_unscaled)
wgs84 <- 4326

littoral_zone <- read_sf('VEG_OKE_1973/VEG_OKE_1973.shp') %>% 
  filter(is.na(VEG_CLASS_) == F) %>%
  filter(VEG_CLASS_ != 'Open Water') %>%
  st_transform(crs = epsg) %>% st_union()

lake_bound <-  read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200') %>%
  st_transform(epsg)

# Quadrants ----
lake_centroid <- st_centroid(st_transform(lake_bound,wgs84))
lake_extent <- extent(st_transform(lake_bound,wgs84))

xmn <- lake_extent[1]
xmx <- lake_extent[2]
ymn <- lake_extent[3]
ymx <- lake_extent[4]
xmd <- lake_centroid$geometry[[1]][1]
ymd <- lake_centroid$geometry[[1]][2]
p1 <- st_point(c(xmn,ymx)) 
p2 <- st_point(c(xmd,ymx)) 
p3 <- st_point(c(xmx,ymx)) 
p4 <- st_point(c(xmn,ymd))
p5 <- st_point(c(xmd,ymd)) 
p6 <- st_point(c(xmx,ymd)) 
p7 <- st_point(c(xmn,ymn)) 
p8 <- st_point(c(xmd,ymn)) 
p9 <- st_point(c(xmx,ymn)) 

nw_quad_84 <- st_polygon(x = list(rbind(p1,p2, p5,p4, p1))) %>% 
  st_sfc(crs = wgs84) %>% st_intersection(st_transform(lake_bound,wgs84))
ne_quad_84 <- st_polygon(x = list(rbind(p2,p3, p6,p5, p2))) %>% 
  st_sfc(crs = wgs84) %>% st_intersection(st_transform(lake_bound,wgs84))
sw_quad_84 <- st_polygon(x = list(rbind(p4,p5, p8,p7, p4))) %>% 
  st_sfc(crs = wgs84) %>% st_intersection(st_transform(lake_bound,wgs84))
se_quad_84 <- st_polygon(x = list(rbind(p5,p6, p9,p8, p5))) %>% 
  st_sfc(crs = wgs84) %>% st_intersection(st_transform(lake_bound,wgs84))

# Transform quadrants
nw_quad <- st_transform(nw_quad_84,epsg)
ne_quad <- st_transform(ne_quad_84,epsg)
sw_quad <- st_transform(sw_quad_84,epsg)
se_quad <- st_transform(se_quad_84,epsg)

## DBHYDRO ----
stations <- c('EASTSHORE','FEBIN','FEBOUT',
              'KBARSE','KISSR0.0',
              'L001','L004','L005','L006','L007','L008',
              'LZ2','LZ25A','LZ30','LZ40',
              'MBOXSOU','MH16000','MH24000','MH32000',
              'NCENTER','NES191','NES135', 'OISLAND',
              'PALMOUT','PALMOUT1','PALMOUT2','PALMOUT3','PELBAY3',
              'POLESOUT','POLESOUT1','POLESOUT2','POLESOUT3',
              'RITTAE2','S308C','TIN13700','TIN16100')

length(stations)

# I'll have to use ncol() to ensure the number of columns is length(stations) + 1
get_wq(station_id = stations,
       date_min = '2016-01-01',
       date_max = '2021-12-31',
       test_name = 'CHLOROPHYLL-A(LC)') %>% as_tibble() %>% 
  mutate(date = date(date)) %>%
  filter(date >= ymd('20160501') &
           date < ymd('20210501')) -> lake_chla

lake_chla %>%
  gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = chla) %>%
  mutate(obs = !is.na(chla)) %>%
  group_by(month(date),year(date),station_id) %>%
  summarize(n_obs = sum(obs,na.rm = T)) %>% 
  filter(n_obs > 1) %>% pull(`year(date)`) %>% unique() # There are 159 instances of stations being sampled twice in a month in 2020

# Spatial data ----
all_station_coords <- read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names()

all_station_coords %>% 
  filter(station %in% stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>% 
  dplyr::select(station, start_date, end_date, geometry) -> lake_station_coords

## Map ALL stations
ggplot() +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent') + 
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent') + 
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'steelblue4',size = 2) +
  # geom_sf(data = st_transform(sle_station_coords[2,],wgs84),color = '#b3bf75',size = 2) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Determine quadrants ----
tibble(
  site = lake_station_coords$station,
  nw_quad = as.vector(st_intersects(st_transform(lake_station_coords,wgs84),
                                    nw_quad_84,sparse = F)),
  ne_quad= as.vector(st_intersects(st_transform(lake_station_coords,wgs84),
                                   ne_quad_84,sparse = F)),
  sw_quad = as.vector(st_intersects(st_transform(lake_station_coords,wgs84),
                                    sw_quad_84,sparse = F)),
  se_quad = as.vector(st_intersects(st_transform(lake_station_coords,wgs84),
                                    se_quad_84,sparse = F)),
  quad = NA
) -> quad_assignments

quad_assignments$quad[which(quad_assignments$nw_quad == T)] <- 'Northwest'
quad_assignments$quad[which(quad_assignments$ne_quad == T)] <- 'Northeast'
quad_assignments$quad[which(quad_assignments$sw_quad == T)] <- 'Southwest'
quad_assignments$quad[which(quad_assignments$se_quad == T)] <- 'Southeast'

quad_assignments %>% rename(station_id = site) %>% dplyr::select(station_id,quad) %>%
  arrange(station_id) -> quad_assignments

## Plot ----
lake_chla %>%
  gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = chla) %>%
  mutate(station_id = factor(x = station_id, labels = stations))  %>%
  mutate(obs = !is.na(chla)) %>%
  group_by(station_id) %>%
  summarize(n_obs = sum(obs,na.rm = T)) %>%
  full_join(quad_assignments) %>% 
  ggplot() +
  geom_col(aes(station_id,n_obs,fill = quad)) +
  geom_hline(yintercept = 30,linetype = 'twodash',color = 'coral2',size = .9)+
  geom_hline(yintercept = 40,linetype = 'twodash',color = 'goldenrod2',size = .9) +
  xlab('Station ID') +
  ylab('Number of Observations \n(some bimonthly sampling in 2020)') +
  scale_fill_manual(name = 'Quadrant',
                    values =  as.vector(wes_palette(name = 'GrandBudapest2'))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))

# lake_chla %>%
#   gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = chla) %>%
#   mutate(station_id = factor(x = station_id, labels = stations))  %>%
#   mutate(obs = !is.na(chla)) %>%
#   group_by(station_id) %>%
#   summarize(n_obs = sum(obs,na.rm = T)) %>%
#   filter(n_obs > 30) %>% pull(station_id) %>% as.character() -> stations_30

lake_chla %>%
  gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = chla) %>%
  mutate(station_id = factor(x = station_id, labels = stations))  %>%
  mutate(obs = !is.na(chla)) %>%
  group_by(station_id) %>%
  summarize(n_obs = sum(obs,na.rm = T)) %>%
  filter(n_obs > 40) %>% pull(station_id) %>% as.character() -> stations_40

ggplot() +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent') + 
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent') + 
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = st_transform(
    filter(lake_station_coords,station %in% stations_40),
    wgs84), 
    color = 'coral2',size = 3.5) +
  geom_sf(data = st_transform(
    filter(lake_station_coords,station %in% stations_40),
    wgs84), 
    color = 'goldenrod2',size = 3.5) +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'steelblue4',size = 2) +
  # geom_sf(data = st_transform(sle_station_coords[2,],wgs84),color = '#b3bf75',size = 2) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Plot out data availability of stations_40 ----
stations_40

lake_chla %>%
  gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = chla) %>%
  mutate(station_id = factor(x = station_id, labels = stations))  %>%
  mutate(obs = !is.na(chla)) %>%
  filter(station_id %in% stations_40) %>% 
  group_by(year(date),month(date),station_id) %>%
  summarize(n_obs = sum(obs, na.rm = T)) %>%
  mutate(day = 1) %>% # nominal
  unite(`year(date)`,`month(date)`,day, col = 'date', sep = '-') %>%
  mutate(date = ymd(date)) %>%
  ggplot() +
  geom_tile(aes(x = date,y = station_id,fill = as.character(n_obs)),alpha = .8) +
  scale_fill_manual(values = rev(wes_palette('Cavalcanti1')[2:5])) +
  scale_x_date(date_breaks = '3 months',date_minor_breaks = '1 month',
               date_labels = '%b %Y') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

ggplot() +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent') + 
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent') + 
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = st_transform(
    filter(lake_station_coords,site %in% stations_40),
    wgs84), 
    color = 'plum2',size = 3.5) +
  geom_sf(data = st_transform(
    filter(lake_station_coords,site == 'PELBAY3'),
    wgs84), 
    color = 'coral2',size = 3.5) +
  geom_sf(data = st_transform(
    filter(lake_station_coords,site == 'RITTAE2'),
    wgs84), 
    color = 'goldenrod2',size = 3.5) +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'steelblue4',size = 2) +
  # geom_sf(data = st_transform(sle_station_coords[2,],wgs84),color = '#b3bf75',size = 2) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))
  
# Apply land qa mask ----

# Re-define stations as those that have no more than one month's worth of missing data
lake_stations <- c("KISSR0.0",
                   "L001",
                   "L004",
                   "L005",
                   "L006",
                   "L007",
                   "L008",
                   "LZ2",
                   "LZ25A",
                   "LZ30",
                   "LZ40",
                   "PALMOUT",
                   "POLESOUT",
                   "S308C")

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% lake_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> lake_station_coords

load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")

land_adjacency_flag <- raster('LandAdjacenyQA.tif') %>% crop(extent(stack_daily_scaled)) 

littoral_zone <- read_sf('VEG_OKE_1973/VEG_OKE_1973.shp') %>% 
  na.omit() %>%
  filter(VEG_CLASS_ != 'Open Water') %>%
st_transform(crs = epsg) %>% st_union()

littoral_zone %>% pull(VEG_CLASS) %>% unique()

raster::extract(land_adjacency_flag, lake_station_coords) %>% 
  as.tibble() %>% rename(qa_flag = value) %>%
  bind_cols(lake_station_coords) %>% 
  filter(qa_flag == 1) -> lake_station_coords_filtered

ggplot() +
  geom_sf(data = littoral_zone,color = 'transparent', fill = 'slategray') +
  layer_spatial(data = land_adjacency_flag) +
  geom_sf(data = lake_station_coords_filtered, 
          aes(geometry = geometry,color = station_id)) +
  scale_fill_gradientn(colors = c('slategray3','transparent'),
                       na.value = 'transparent') +
  theme_minimal()

raster::extract(focal(land_adjacency_flag,w = matrix(1/9, nc=3, nr=3)),lake_station_coords) %>% 
  as.tibble() %>% rename(qa_flag = value) %>%
  bind_cols(lake_station_coords)

# Blur land adjacency mask with a focal filter
