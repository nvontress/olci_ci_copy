# New script to figure out which objects in the NHD dataset I need :) 
# The goal: stop crashing R

# Date: 10/25/2021
# Author: Natalie Von Tress

# Set up work space ----

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
# library(lubridate)
library(raster)
library(sf)
library(ggspatial)
# library(ggnewscale)
# library(scales)
# library(beepr)
# library(colormap)
library(janitor)
# library(dbhydroR)
# library(ggpubr)
library(cowplot)

# Set working directory
setwd('/Users/natalie/Documents/select_data')

# Load data ----

# Loads raster stacks, processed and unprocessed, weekly and daily
load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")

# Spatial data ----
## Set coordinate system
epsg <- crs(stack_daily_unscaled)
wgs84 <- 4326
albers <- 5070

lake_bound  <- read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200') %>%
  st_transform(epsg)

florida_bound <- read_sf('cb_2019_us_state_500k/cb_2019_us_state_500k.shp') %>%
  filter(NAME == 'Florida') %>% 
  st_transform(crs = epsg)

not_conus <- c('United States Virgin Islands','Hawaii','Commonwealth of the Northern Mariana Islands','Puerto Rico','Guam','American Samoa','Alaska')

conus_bound <- read_sf('cb_2019_us_state_500k/cb_2019_us_state_500k.shp') %>% 
  filter(NAME != 'United States Virgin Islands' &
           NAME != 'Hawaii' &
           NAME != 'Commonwealth of the Northern Mariana Islands' &
           NAME != 'Puerto Rico' &
           NAME != 'Guam' &
           NAME != 'American Samoa' &
           NAME != 'Alaska' ) %>%
  st_union()%>% 
  st_transform(crs = epsg)

# read_sf('cb_2019_us_state_500k/cb_2019_us_state_500k.shp') %>% pull(NAME)

nhd_flowline <- read_sf('Shape_NHD/NHDFlowline.shp') %>% st_transform(crs = epsg)

canals_and_streams <- read_sf('AHED_Canals_and_Streams/AHED_Canals_and_Streams.shp') %>% st_transform(epsg)

canals <- read_sf('SFWMD_Canals/SFWMD_Canals.shp') %>% st_transform(epsg)

ggplot() + 
  geom_sf(data = florida_bound) +
  geom_sf(data = lake_bound) +
  geom_sf(data = canals) +
  theme_minimal()
  

watershed <- read_sf('Shape_NHD/WBDHU4.shp') %>% st_transform(crs = epsg)
huc_6 <- read_sf('Shape_NHD/WBDHU6.shp') %>% st_transform(crs = epsg)
huc_8 <- read_sf('Shape_NHD/WBDHU8.shp') %>% st_transform(crs = epsg)
huc_10 <- read_sf('Shape_NHD/WBDHU10.shp') %>% st_transform(crs = epsg)

ggplot() +
  geom_sf(data = florida_bound) +
  geom_sf(data = watershed) +
  geom_sf(data = filter(huc_10,name %in% c('Caloosahatchee River Headwaters','South Fork of the St. Lucie River-Indian River','St. Lucie Canal','North Fork of the St. Lucie River','Upper Caloosahatchee River')),aes(fill = name),alpha = 1) +
  geom_sf(data = huc_6,alpha = .3) +
  geom_sf(data = lake_bound)

littoral_zone <- read_sf('VEG_OKE_1973/VEG_OKE_1973.shp') %>% 
  filter(is.na(VEG_CLASS_) == F) %>%
  filter(VEG_CLASS_ != 'Open Water') %>%
  st_transform(crs = epsg) %>% st_union()

florida_lakes <- read_sf('Florida_Lakes/Florida_Lakes.shp') %>% st_transform(crs = epsg)

# State map ----
aliases <- c(#"BOLLES CANAL",
             "C-43","CALOOSAHATCHEE CANAL","CALOOSAHATCHEE RIVER",
             # "FEC",
             # "HARNEY POND CANAL",
             # "HILLSBORO CANAL",
             # "HILLSBORO RIVER",
             # "HOLLYWOOD CANAL",
             "INDIAN PRAIRIE CANAL",
             # "KISSIMMEE RIVER",
             # "L-2W CANAL",
             # "LOXAHATCHEE CANAL",
             # "MIAMI CANAL","MIAMI RIVER","NEW RIVER",
             # # "NINE MILE CANAL","NINE MILE STUB CANAL",
             # "NORTH FORK MIAMI RIVER","NORTH NEW RIVER CANAL",
             # "SOUTH NEW RIVER CANAL",
             "ST. LUCIE CANAL")#,"WEST PALM BEACH CANAL")

names <- c("NORTH FORK ST LUCIE RIVER","SOUTH FORK ST LUCIE RIVER",
           "ST LUCIE ESTUARY") #,
           # "TAYLOR CREEK","FISHEATING CREEK",
           # # "NEW RIVER",
           # "SOUTH FORK NEW RIVER",
           # "SOUTH NEW RIVER CANAL",
           # "C-39A","C-41A","C-41")

 # ggplot() +
 #   # State watershed, and lake boundaries
 #   geom_sf(data = florida_bound,alpha = .4) +
 #   geom_sf(data = st_intersection(florida_bound,watershed),alpha = .4) +
 #   geom_sf(data = lake_bound) +
 #  
 #   geom_sf(data = filter(canals_and_streams,ALIAS %in% aliases | NAME %in% names),
 #           aes(color = NAME)) +
 #   geom_sf(data = filter(nhd_flowline,gnis_name == "Kissimmee River")) +
 #   geom_sf(data = florida_lakes%>% filter(NAME == "Lake Kissimmee" |
 #                                            NAME == "Lake Istokpoga")) +
 #   theme_minimal()

ggplot() +
  # State, watershed, and lake bounds
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_intersection(florida_bound,watershed),
          color = 'transparent',
          fill = 'lightsteelblue1',
          alpha = .5) +
  geom_sf(data = st_transform(filter(huc_6,name == 'Kissimmee'),wgs84),
          alpha = .65,
          lwd = 0, fill = 'lightsteelblue') +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  # Canals
   geom_sf(data = st_transform(
   st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                            NAME %in% names),lake_bound)),
    wgs84),
   color = 'darkslategray4') +
   geom_sf(data = st_transform(
     st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
     wgs84),
     color = 'darkslategray4') +
  geom_sf(data = canals, color = 'darkslategray4') +
  # Misc lakes
  geom_sf(data = st_transform(st_intersection(florida_lakes,
                                              filter(huc_6,name == 'Kissimmee')),
                              wgs84),
          lwd = 0,fill = 'darkslategray4') +
  # Littoral zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # Add boundary outlines
  geom_sf(data = st_transform(florida_bound,wgs84), fill = 'transparent') +
  geom_sf(data = st_transform(lake_bound,wgs84), fill = 'transparent',
          color = 'gray30') +
  # geom_rect(aes(
  #   xmin = -81.17,
  #   xmax = -80.55,
  #   ymin = 26.63,
  #   ymax = 27.25
  # ), fill = 'transparent', color = 'lightcyan4', size = .5) +
  # coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() +
  theme(text = element_text(size = 12),        
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ggsave(
#   filename = 'state_map.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 4,
#   height = 4,
#   units = 'in'
# )

# Inset map ----
conus_map <-
  ggplot() + 
  geom_sf(data = st_transform(conus_bound,albers),fill = 'white', alpha = 0.4, size = .3) +
  geom_sf(data = st_transform(
    st_as_sfc(st_bbox(st_buffer(florida_bound,90000))),
    albers), fill = 'transparent',color = 'red3', size = .5) +
  geom_sf(data = st_transform(florida_bound,albers), alpha = .8, size = .3) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.1,pad_x = unit(-0.03,'in')) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(0.9, 'cm'),
                         pad_x = unit(-0.25,'in'), pad_y = unit(0.2,'in')) +
  theme_minimal() +
  theme(text = element_blank(),
        plot.margin=grid::unit(rep(2,4),"pt"),
        plot.background = element_rect(fill = 'white'))

state_map <- ggplot() +
  # State, watershed, and lake bounds
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8, size = 0.4) +
  geom_sf(data = st_intersection(florida_bound,watershed),
          color = 'transparent',
          fill = 'lightsteelblue1',
          alpha = .5, size = 0.4) +
  geom_sf(data = st_transform(filter(huc_6,name == 'Kissimmee'),wgs84),
          alpha = .65,
          lwd = 0, fill = 'lightsteelblue3', size = 0.4) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30', size = 0.4) +
  # Canals
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4', size = 0.4) +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4', size = 0.4) +
  geom_sf(data = canals, color = 'darkslategray4', size = 0.4) +
  # Misc lakes
  geom_sf(data = st_transform(st_intersection(florida_lakes,
                                              filter(huc_6,name == 'Kissimmee')),
                              wgs84),
          lwd = 0,fill = 'darkslategray4', size = 0.4) +
  # Littoral zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent', size = 0.4) +
  # Add boundary outlines
  geom_sf(data = st_transform(florida_bound,wgs84), fill = 'transparent', size = 0.5) +
  geom_sf(data = st_transform(lake_bound,wgs84), fill = 'transparent',
          color = 'gray30', size = 0.4) +
  # Formatting
  annotation_scale(location = 'tr', height = unit(0.2,'cm'), width_hint = 0.2) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(2,'cm'),
                         pad_x = unit(4.1,'in'), pad_y = unit(3.4,'in')) +
  theme_minimal() +
  theme(text = element_text(size = 12),        
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggdraw() +
  draw_plot(state_map) +
  draw_plot(conus_map, x = 0.2, y = 0.1, width = 0.3, height = 0.3)

ggsave(
  filename = 'inset_map_redo.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 5,
  units = 'in'
)

# Other figures ----

# sle_stations <- c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11')

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>% arrange(site) %>%
  filter(station == 'S308C',
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> s308c_coords

# sle_station_coords$station_id <- factor(sle_station_coords$station_id, 
#                                         levels =  c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
#                                         ordered = T)

lake_stations <- c("KISSR0.0","L001","L004","L005","L006","L007","L008",
                   # "LZ2",
                   "LZ25A","LZ30","LZ40"#,
                   # "PALMOUT","POLESOUT","S308C"
                   )

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% lake_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> lake_station_coords

near_qa <- c("LZ2","PALMOUT","POLESOUT","S308C")

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% near_qa,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> qa_stations

all_stations <- c('EASTSHORE','FEBIN','FEBOUT',
                  'KBARSE','KISSR0.0',
                  'L001','L004','L005','L006','L007','L008',
                  'LZ2','LZ25A','LZ30','LZ40',
                  'MBOXSOU','MH16000','MH24000','MH32000',
                  'NCENTER','NES191','NES135', 'OISLAND',
                  'PALMOUT','PALMOUT1','PALMOUT2','PALMOUT3','PELBAY3',
                  'POLESOUT','POLESOUT1','POLESOUT2','POLESOUT3',
                  'RITTAE2','S308C','TIN13700','TIN16100')

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% all_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> all_station_coords



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


ggplot() +
  # Boundaries
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  # Canals and rivers
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  # Littoral Zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # Quadrants
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  # Stations
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(color = 'palegreen4')) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = 'steelblue4')) +
  scale_color_identity(name = 'Sampling stations:',
                       breaks = c('palegreen4','steelblue4'),
                       labels = c('Lake stations\n(Chlorophyll data)',
                                  'S-308\n(Flow data)'),
                       #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
                       guide = 'legend') +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
  theme_minimal() +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_minimal,
                          pad_x = unit(0.08,'in'), pad_y = unit(0.3,'in')) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom'
        )

ggsave(
  filename = 'lake_stations_plot_with_quads.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5.5,
  height = 5,
  units = 'in'
)



ggplot() +
  # Boundaries
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  # Canals and rivers
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  # Littoral Zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # Stations
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(color = 'steelblue4')) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = 'transparent')) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = 'palegreen4')) +
  scale_color_identity(name = 'Sampling stations:',
                       breaks = c('steelblue4','transparent','palegreen4'),
                       labels = c('Lake stations\n(Chlorophyll data)',
                                  ' ',
                                  'S-308\n(Chlorophyll and flow data)'),
                       guide = 'legend') +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +

  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right'
  )

ggsave(
  filename = 'lake_stations_plot_without_quads.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5,
  height = 3,
  units = 'in'
)

# Site map for pptx ----
ggplot() +
  # Boundaries
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  # Canals and rivers
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  # Littoral Zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # Quadrants
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  # Stations
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(color = 'steelblue4')) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = 'palegreen4')) +
  scale_color_identity(name = 'Sampling stations:',
                       breaks = c('steelblue4','palegreen4'),
                       labels = c('Lake stations\n(Chlorophyll data)',
                                  'S-308\n(Chlorophyll and flow data)'),
                       #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
                       guide = 'legend') +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
  theme_minimal() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom'
  )

ggsave(
  filename = 'lake_stations_plot_with_quads_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 7,
  height = 6,
  units = 'in'
)

# Without quads
ggplot() +
  # Boundaries
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  # Canals and rivers
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  # Littoral Zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # Stations
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(color = 'steelblue4')) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = 'palegreen4')) +
  scale_color_identity(name = 'Sampling stations:',
                       breaks = c('steelblue4','palegreen4'),
                       labels = c('Lake stations\n(Chlorophyll data)',
                                  'S-308\n(Chlorophyll and flow data)'),
                       guide = 'legend') +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
  theme_minimal() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

ggsave(
  filename = 'lake_stations_plot_without_quads_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 7,
  height = 6,
  units = 'in'
)

# With all stations ----


ggplot() +
  # Boundaries
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  # Canals and rivers
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  geom_sf(data = st_difference(canals,lake_bound), color = 'darkslategray4', size = 0.4) +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4') +
  # Littoral Zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # Quadrants
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  # Stations
  # goldenrod,steelblue4,palegreen4,plum2,coral2
  geom_star(aes(x = unlist(pull(st_transform(s308c_coords,wgs84),geometry))[1],
                y = unlist(pull(st_transform(s308c_coords,wgs84),geometry))[2],
                # label="★",
                color = 'steelblue4'),
            fill = 'steelblue4',
     size=10
     ) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = 'steelblue4')) + 
  geom_sf(data = st_transform(filter(all_station_coords,
                                     !station_id %in% qa_stations$station_id &
                                       !station_id %in% lake_station_coords$station_id),wgs84),size = 2,
          aes(color = 'goldenrod'), shape = 25, fill = 'goldenrod') +  
  geom_sf(data = st_transform(qa_stations,wgs84),size = 2,
          aes(color = 'coral2'), shape = 24, fill = 'coral2') + 
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(color = 'palegreen4')) +
  scale_color_identity(name = 'Sampling stations:',
                       breaks = c('goldenrod','coral2','palegreen4','steelblue4'),
                       labels = c('Not enough samples\nduring study period',
                                  'Near mixed\nland/water pixels',
                                  'All requirements met',
                                  'S-308 (Flow data)'
                                  ),
                       #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
                       guide = 'legend') +
  coord_sf(xlim = c(-81.17,-80.55),ylim = c(26.68,27.2)) +
  # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.2) +
  annotation_north_arrow(style = north_arrow_minimal,
                         pad_x = unit(-0.2,'cm'), 
                         pad_y = unit(0.25,'in')) +
  guides(color = guide_legend(override.aes = list(starshape = c(23,11,15,1),
                                                  shape = c(25,24,19,42),
                                                  size = c(2,2,2,10),
                                                  fill = c('goldenrod','coral2',
                                                           'palegreen4','steelblue4')))) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        legend.text = element_text(
          margin = margin(b = 5,t = 5, unit = "pt")))

ggsave(
   filename = 'all_stations_plot_with_quads.png',
   path = '/Users/natalie/Documents/thesis_writing/r_plots',
   width = 6.5,
   height = 4.25,
   units = 'in'
 )

# blank outlines for conceptual figure
ggplot() +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
        color = 'gray30') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())

ggsave(
  filename = 'blank_map.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 3,
  height = 3,
  units = 'in'
)

ggplot() +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())

ggsave(
  filename = 'blank_quad.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 1.5,
  height = 1.5,
  units = 'in'
)

ggplot() +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = st_transform(filter(lake_station_coords,station_id %in% c("L004","L001")),wgs84),size = 1.5) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())

ggsave(
  filename = 'blank_quad_stations.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 1.5,
  height = 1.5,
  units = 'in'
)

# Figure for internship chapter
read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station == 'S79',
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(site, station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> s79_coords

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station == 'C44S80',
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(site, station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> s80_coords

ggplot() + 
  geom_sf(data = st_transform(florida_bound,wgs84),fill = 'gray100', alpha = .8, size = 0.4) +
  geom_sf(data = st_intersection(florida_bound,watershed),
          color = 'transparent',
          fill = 'lavenderblush2',
          alpha = .5, size = 0.4) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30', size = 0.4) +
  # Canals
  geom_sf(data = st_transform(
    st_union(st_difference(filter(canals_and_streams,ALIAS %in% aliases | 
                                    NAME %in% names),lake_bound)),
    wgs84),
    color = 'darkslategray4', size = 0.4) +
  geom_sf(data = st_transform(
    st_union(st_difference(filter(st_zm(nhd_flowline),gnis_name == "Kissimmee River"),lake_bound)),
    wgs84),
    color = 'darkslategray4', size = 0.4) +
  geom_sf(data = canals, color = 'darkslategray4', size = 0.4) +
  # Misc lakes
  geom_sf(data = st_transform(st_intersection(florida_lakes,
                                              filter(huc_6,name == 'Kissimmee')),
                              wgs84),
          lwd = 0,fill = 'darkslategray4', size = 0.4) +
  # Littoral zone
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent', size = 0.4) +
  # Add boundary outlines
  geom_sf(data = st_transform(st_intersection(florida_bound,watershed),wgs84), fill = 'transparent', size = 0.4, color = 'black') +
  geom_sf(data = st_transform(lake_bound,wgs84), fill = 'transparent',
          color = 'gray30', size = 0.4) +  
  # Stations
  geom_sf(data = s79_coords) + geom_sf(data = s80_coords) +
  # Formatting
  coord_sf(xlim = c(-83,-80),ylim = c(24.5,28.5)) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.1,pad_x = unit(-0.03,'in')) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(0.9, 'cm'),
                         pad_x = unit(-0.25,'in'), pad_y = unit(0.2,'in')) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        legend.text = element_text(
          margin = margin(b = 5,t = 5, unit = "pt")))

ggsave(
  filename = 'south_fl_map.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 4,
  height = 5,
  units = 'in'
)
