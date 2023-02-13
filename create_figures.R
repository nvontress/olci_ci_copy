# Create figures

# Set up work space ----

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(raster)
library(sf)
library(ggspatial)
library(ggnewscale)
library(scales)
library(beepr)
library(colormap)
library(janitor)
library(dbhydroR)
library(ggpubr)
library(cowplot)

# Set working directory
setwd('/Users/natalie/Documents/select_data')

# Read in data ----

# Loads raster stacks, processed and unprocessed, weekly and daily
load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")
# load("~/Documents/R projects/olci_ci_copy/frequency_analysis_ci_cyano_COPY_workspace.RData")

# Spatial data
## Set coordinate system
epsg <- crs(stack_daily_unscaled)
wgs84 <- 4326

lake_bound  <- read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200') %>%
  st_transform(epsg)

florida_bound <- read_sf('cb_2019_us_state_500k/cb_2019_us_state_500k.shp') %>% 
   filter(NAME == 'Florida') %>% 
   st_transform(crs = epsg)

waterways <- read_sf('Waterways_Florida-shp/Waterways_Florida.shp') 
# For the Caloosahatchee, I think you should be able to filter for OBJECTID == 79 (canal), OBJECTID == 399 (river), and OBJECTID == 349 (estuary). An interactive map of the data can be found here.

waterways %>% filter(
  NAME == "St. Lucie Canal  (C-44)" |
    NAME == "South Fork St. Lucie River" |
    NAME == "St. Lucie River" |
    OBJECTID == 762
) %>% st_transform(crs = epsg) %>% st_intersection(y = florida_bound)  -> sle

waterways %>% filter(
  OBJECTID == 79 |
    OBJECTID == 399 |
    OBJECTID == 349 
) %>% st_transform(crs = epsg) %>% st_intersection(y = florida_bound) -> caloosahatchee

land_adjacency_flag <- raster('LandAdjacenyQA.tif') %>% crop(extent(stack_daily_scaled)) 

littoral_zone <- read_sf('VEG_OKE_1973/VEG_OKE_1973.shp') %>% 
  filter(is.na(VEG_CLASS_) == F) %>%
  filter(VEG_CLASS_ != 'Open Water') %>%
  st_transform(crs = epsg) %>% st_union()

# # sle_station_coords <- st_transform(sle_station_coords, epsg)
 sle_stations <- c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11')

 read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>% arrange(site) %>%
   filter(station %in% sle_stations,
          activity_t == 'Chemistry',
          activity_s == 'Surface Water Grab') %>%
   dplyr::select(station, geometry) %>%
   rename(station_id = station) %>%
   st_transform(crs = epsg) -> sle_station_coords
 
 sle_station_coords$station_id <- factor(sle_station_coords$station_id, 
                                         levels =  c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                                         ordered = T)

# Chosen lake stations
lake_stations <- c("KISSR0.0",
                   "L001",
                   "L004",
                   "L005",
                   "L006",
                   "L007",
                   "L008",
                   # "LZ2",
                   "LZ25A",
                   "LZ30",
                   "LZ40"
                   # "PALMOUT",
                   # "POLESOUT",
                   # "S308C"
                   )

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% lake_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> lake_station_coords

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station == 'S308C',
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> s308c_coords

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

# Transform quadrants]
nw_quad <- st_transform(nw_quad_84,epsg)
ne_quad <- st_transform(ne_quad_84,epsg)
sw_quad <- st_transform(sw_quad_84,epsg)
se_quad <- st_transform(se_quad_84,epsg)

# Chlorophyll data

get_wq(station_id = sle_stations,
       date_min = '2016-01-01',
       date_max = '2020-12-31',
       test_name = 'CHLOROPHYLL-A(LC)') %>% as_tibble() -> sle_chla

sle_chla %>%
  gather(colnames(sle_chla)[2:ncol(sle_chla)],key = station_id,value = value) %>%
  mutate(station_id = str_extract(station_id,'.+(?=_CHLOROPHYLL)'),
         date = date(date)) %>% 
  mutate(station_id = factor(station_id, 
                             levels =  c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                             ordered = T)) %>%
  rename(collection_date = date) %>%
  filter(!is.na(value)) -> sle_chla

get_wq(station_id = lake_stations,
       date_min = '2016-01-01',
       date_max = '2020-12-31',
       test_name = 'CHLOROPHYLL-A(LC)') %>% as_tibble() -> lake_chla

lake_chla %>%
  gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = value) %>%
  mutate(station_id = str_extract(station_id,'.+(?=_CHLOROPHYLL)'),
         date = date(date)) %>% 
  rename(collection_date = date) %>%
  filter(!is.na(value)) -> lake_chla

s308_usgs <- read_csv('NWISData/20160101-20201231_s308.csv')

# Pull dates from file names ----

## Daily
dates_daily <- rep(NA, length(filenames_daily))
for(i in 1:length(dates_daily)){
  # i <- 1
  
  yr <- str_split(filenames_daily,pattern = '\\.')[[i]] %>% first() %>% 
    str_extract('[:digit:]{4}')
  
  yr_prev <- as.character(as.numeric(yr)-1)
  
  origin_date <- str_c(yr_prev,'-12-31') %>% ymd()
  
  str_split(filenames_daily,pattern = '\\.')[[i]] %>% first() %>% 
    str_extract('[:digit:]{3}$') %>% as.numeric() %>%
    as.Date(origin = origin_date) -> dates_daily[i]
}

dates_daily %>% as_date() -> dates_daily

## Weekly
dates_weekly_start <- rep(NA, length(filenames_weekly))
dates_weekly_end <- rep(NA, length(filenames_weekly))

for(i in 1:length(dates_weekly_end)){
  # i <- 1
  
  yr <- str_split(filenames_weekly,pattern = '\\.')[[i]] %>% first() %>% 
    str_extract('[:digit:]{4}')
  
  yr_prev <- as.character(as.numeric(yr)-1)
  
  origin_date <- str_c(yr_prev,'-12-31') %>% ymd()
  
  str_split(filenames_weekly,pattern = '\\.')[[i]] %>% first() %>% 
    str_extract('[:digit:]{3}$') %>% as.numeric() %>%
    as.Date(origin = origin_date) -> dates_weekly_end[i]
  
  dates_weekly_end[i] - 6 -> dates_weekly_start[i]
}

dates_weekly_start <- as_date(dates_weekly_start)
dates_weekly_end <- as_date(dates_weekly_end)
beep(2)

# Create Plots ----

# Time Series ----
plot_dbhydro <- function(dbhydro_report = sle_chla) {
  dbhydro_report %>% 
    ggplot() + 
    geom_point(aes(x = collection_date, 
                   y = value, 
                   color = station_id),
               size = 2,
               alpha = 0.85)
}

plot_chla_and_flow <- function(year_bounds = year_bounds, 
                               scale_flow = 1/20,
                               chla = sle_chla, flow = s308_usgs){
  chla$station_id <- factor(chla$station_id, 
                            levels =  c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                            ordered = T)
  
  chla %>% filter(year(collection_date) >= year_bounds[1],
                  year(collection_date) <= year_bounds[2]) %>%
    plot_dbhydro() + 
    geom_hline(yintercept = c(11,40), color = 'seagreen',linetype = 'longdash') + 
    geom_line(data = s308_usgs %>%
                  group_by(date(dateTime)) %>% 
                  summarize(daily_flow = mean(flow_cfs)) %>%
                  mutate(avg_7d = zoo::rollmean(daily_flow, k = 7, fill = NA)),
                mapping = aes(x = `date(dateTime)`,y = avg_7d*scale_flow),
                color = 'black', size = .75) + 
    xlab('Date') +
    scale_y_continuous(name = 'Clorophyll-a concentation (µg/L)',
                       sec.axis = sec_axis(trans = ~. / scale_flow, 
                                           name = 'Discharge to SLC (cfs)')) }

plot_dbhydro_shape <- function(dbhydro_report = sle_chla) {
  dbhydro_report %>% 
    ggplot() + 
    geom_point(aes(x = collection_date, 
                   y = value, 
                   shape = station_id),
               size = 2,
               alpha = 0.7)
}

plot_chla_and_flow_shape <- function(year_bounds = year_bounds, 
                                     scale_flow = 1/20,
                                     chla = sle_chla, flow = s308_usgs){
  chla$station_id <- factor(chla$station_id, 
                            levels =  c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                            ordered = T)
  
  chla %>% filter(year(collection_date) >= year_bounds[1],
                  year(collection_date) <= year_bounds[2]) %>%
    plot_dbhydro_shape() + 
    geom_hline(yintercept = c(11,40), color = 'seagreen',linetype = 'longdash') + 
    stat_smooth(data = (flow %>% filter(year(dateTime) >= year_bounds[1],
                                        year(dateTime) <= year_bounds[2])),
                method = 'loess',
                se = F,
                span = .05,
                mapping = aes(x = date(dateTime),y = flow_cfs*scale_flow),
                color = 'black', size = .75) + 
    xlab('Date') +
    scale_y_continuous(name = 'Clorophyll-a concentation (µg/L)',
                       sec.axis = sec_axis(trans = ~. / scale_flow, 
                                           name = 'Discharge to SLC (cfs)')) }

# Mess with colors ----
scales::show_col(
 colormap(
   colormap = colormaps$velocity_blue,
   nshades = 7
 )
)
scales::show_col(
  colormap(
    colormap = colormaps$salinity,
    nshades = 12
  )
)

scales::show_col(c(colors_6,
                   '#10014a'))

scales::show_col(
  colorspace::rainbow_hcl(
    n = 6,
    start = 95,
    end = 330,
    l = 75
  )
)

scales::show_col(
  colorspace::rainbow_hcl(
    n = 5,
    start = 95,
    end = 330,
    l = 85
  )
)

scales::show_col(
  wes_palette(
    name = 'Royal2',
    n = 12,
    type = c('continuous')
  ) #%>% as.vector()
)

# Plot ----
colors <- c(
  'palegreen3',
  'steelblue3',
  'plum',
  'hotpink',
  'chocolate1',
  'tomato3'
)

# New colors
plot_chla_and_flow(year_bounds = c(2016,2020),
                   scale_flow = 1/20) + 
  scale_color_manual(name = 'Sampling \nLocations:',
                     labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                     values = colors) +
  theme_minimal() + 
  coord_cartesian(
    xlim = c(date('2016-01-01'),date('2021-01-02')),
    # ylim = c(-10,170),
    expand = F) +
  scale_x_date(date_labels = '%Y') +
  theme(text = element_text(size = 14), 
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y.right = element_text(margin = margin(0,10,0,5)),
        legend.position = 'bottom',legend.direction = 'horizontal')

# Save
ggsave(
  filename = 'time_series_intro.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 11.5,
  height = 5.5,
  units = 'in'
)

# Using shapes
shapes <- c(8, 2, 15:18)
plot_chla_and_flow_shape(year_bounds = c(2016,2020),
                   scale_flow = 1/20) + 
  scale_shape_manual(name = 'Sampling \nLocations:',
                      labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                      values = shapes) +
  theme_minimal() + 
  coord_cartesian(
    xlim = c(date('2016-01-01'),date('2021-01-02')),
    # ylim = c(-10,170),
    expand = F) +
  scale_x_date(date_labels = '%Y') +
  theme(text = element_text(size = 14), 
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y.right = element_text(margin = margin(0,10,0,5)),
        legend.position = 'bottom',legend.direction = 'horizontal')

# Original colors
plot_chla_and_flow(year_bounds = c(2016,2020),
                   scale_flow = 1/20) + 
  scale_color_manual(name = 'Sampling Locations:',
                     labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                     values = rev(colormap(colormap = colormaps$YIGnBu,nshades = 12)[4:9])) +
  theme_minimal() + 
  coord_cartesian(
    xlim = c(date('2016-01-01'),date('2021-01-02')),
    # ylim = c(-10,170),
    expand = F) +
  scale_x_date(date_labels = '%Y') +
  theme(text = element_text(size = 14), 
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y.right = element_text(margin = margin(0,10,0,5)),
        legend.position = 'bottom',legend.direction = 'horizontal')


plot_chla_and_flow(year_bounds = c(2018,2018),
                   scale_flow = 1/25) +
  scale_color_manual(name = 'Sampling Locations:',
                     labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                     values = rev(colormap(colormap = colormaps$YIGnBu,nshades = 12)[4:9])) +
  theme_minimal() + 
  coord_cartesian(
    xlim = c(date('2018-01-01'),date('2018-12-31')),
    ylim = c(-10,155),
    expand = F) +
  scale_x_date(date_labels = '%B', date_breaks = '1 months') +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = 'none',legend.direction = 'horizontal')

chla <- sle_chla
chla$station_id <- factor(chla$station_id, 
                          levels =  c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                          ordered = T)
year_bounds <- c(2018,2018)
scale_flow <- 1/25
chla %>% filter(year(collection_date) >= year_bounds[1],
                year(collection_date) <= year_bounds[2],
                month(collection_date) >= 5,
                month(collection_date) < 11) %>%
  plot_dbhydro() + geom_hline(yintercept = c(10,50), color = 'seagreen') + 
  stat_smooth(data = (s308_usgs %>% filter(year(dateTime) >= year_bounds[1],
                                           year(dateTime) <= year_bounds[2],
                                           month(dateTime) >= 5,
                                           month(dateTime) < 11) %>%
                        mutate(date = date(dateTime))),
              se = F,
              mapping = aes(x = date,y = s308_usgs*scale_flow),
              color = 'black', size = .75) + 
  xlab('Date') +
  scale_y_continuous(name = 'Clorophyll-a concentation (µg/l)',
                     sec.axis = sec_axis(trans = ~. / scale_flow, 
                                         name = 'Discharge to SLC (cfs)'))+
  scale_color_manual(name = 'Sampling Locations:',
                     labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                     values = rev(colormap(colormap = colormaps$YIGnBu,nshades = 12)[4:9])) +
  theme_minimal() + 
  coord_cartesian(xlim = c(date('2018-05-01'),date('2018-10-31')),ylim = c(-20,95),expand = F) +
  scale_x_date(date_labels = '%B', date_breaks = '1 month') +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = 'bottom',legend.direction = 'horizontal')


# Plots for time-series ----

tibble(week_date_start = dates_weekly_start) %>% 
  mutate(index = row_number()) %>%
  filter(
    week_date_start == date('2018-05-27') |
      week_date_start == date('2018-06-17') |
      week_date_start == date('2018-08-05') |
      week_date_start == date('2018-09-16') |
      week_date_start == date('2018-10-14') 
  ) %>% pull(index) %>% as.vector() -> image_subset

subset_chla_stack <- stack(
  stack_weekly_chla_tom[[image_subset[1]]],
  stack_weekly_chla_tom[[image_subset[2]]],
  stack_weekly_chla_tom[[image_subset[3]]],
  stack_weekly_chla_tom[[image_subset[4]]],
  stack_weekly_chla_tom[[image_subset[5]]]
)

plot_subset <- function(index, raster_stack = subset_chla_stack){
  ggplot() + 
    geom_sf(data = st_transform(lake_bound,wgs84), color = 'grey30', alpha = 0) +
    layer_spatial(data = projectRaster(raster_stack[[index]],crs = wgs84)) + 
    # scale_fill_viridis_c(na.value = 'transparent', limits = c(20,round(max(maxValue(raster_stack)))),
    #                      name = expression(paste('Chlorophyll-',italic(a),' (μg/L)'))) +
     scale_fill_gradientn(limit = c(0, round(max(maxValue(raster_stack)))), na.value = 'transparent', name = expression(paste('Chlorophyll-',italic(a),' (μg/L)')), 
                          colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                          values = rescale(c(0,20,30,round(max(maxValue(raster_stack))))),
                          breaks = c(10,20,30,40)) +
    geom_sf(data = st_transform(sle_station_coords[2,],wgs84),aes(color = station_id), alpha = 1, size = 3) +
    scale_color_manual(values = c("#bbe4b5ff"),
                       name = "Water Quality\nGrab Station ID") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_colorbar(barheight = 10)) +
    # labs(x = 'Longitude', y = 'Latitude') +
    # annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
    # annotation_north_arrow(style = north_arrow_minimal,
    #                        pad_x = unit(0.08,'in'), pad_y = unit(0.3,'in')) +
    theme_minimal() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
}



plot_subset(1)
plot_subset(2)
plot_subset(3)
plot_subset(4)
plot_subset(5)

# Lake O sampling locations with imagery ----
'2018-08-12'

which(dates_weekly_start == '2018-08-12')

stack_weekly_scaled 
 
ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_scaled[[which(dates_weekly_start == '2018-08-12')]],
    crs = wgs84)) + 
  scale_fill_continuous( na.value = 'transparent', name = 'CI multi') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +  
  labs(title = 'CI multi composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 


ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(10,20,40,90)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +
  labs(title = 'Tomlinson chlorophyll composite, \n2016-07-10 to 2016-07-16') +
  theme_minimal() +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(10,20,40,90)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +
  labs(title = 'Tomlinson chlorophyll composite, \n2016-07-10 to 2016-07-16') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_seeg[[which(dates_weekly_start == '2018-08-12')]],
    crs = wgs84)) + 
  scale_fill_gradientn(limit = c(-4, 130), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)',
                       # name = expression(paste('Chlorophyll-',italic(a),' (μg/L)')), 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,50,130)),
                       breaks = c(0,20,50,100)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +
  labs(title = 'Seegers chlorophyll composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10))

# Conceptual model WITH GGSAVE for clicking through in presentations ---- 
which(dates_weekly_start == '2018-08-05')

radius <- 10000
radius <- st_buffer(s308c_coords,dist = radius) 

ggplot() +
  geom_sf(data = st_transform(florida_bound, crs = wgs84), fill = 'grey90', color = 'grey30', alpha = .5) +
  geom_sf(data = st_transform(lake_bound, crs = wgs84),fill = 'grey98', color = 'grey30', alpha = .9) +
  geom_sf(data = st_transform(st_intersection(caloosahatchee,florida_bound),crs = wgs84), color = 'grey30') +
  geom_sf(data = st_transform(st_intersection(sle,florida_bound),crs = wgs84), color = 'grey30') +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(20,40,90)) +
  geom_sf(data = st_transform(s308c_coords,crs = wgs84),aes(color = station_id), size = 2) +
  scale_color_manual(values = 'violetred3',
                     labels = 'S308',
                     name = "Water Quality\nGrab Site ID") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = 'roi_conceptual_1.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 7,
  height = 3,
  units = 'in'
)

ggplot() +
  geom_sf(data = st_transform(florida_bound, crs = wgs84), fill = 'grey90', color = 'grey30', alpha = .5) +
  geom_sf(data = st_transform(lake_bound, crs = wgs84),fill = 'grey98', color = 'grey30', alpha = .9) +
  geom_sf(data = st_transform(st_intersection(caloosahatchee,florida_bound),crs = wgs84), color = 'grey30') +
  geom_sf(data = st_transform(st_intersection(sle,florida_bound),crs = wgs84), color = 'grey30') +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(20,40,90)) +
  geom_sf(data = st_transform(st_intersection(radius,lake_bound), crs = wgs84), alpha = .25, 
          color = 'violetred4',fill = 'red') +
  geom_sf(data = st_transform(lake_bound, crs = wgs84),fill = 'grey98', color = 'grey30', alpha = 0) +
  geom_sf(data = st_transform(s308c_coords,crs = wgs84),aes(color = station_id), size = 2) +
  scale_color_manual(values = 'violetred3',
                     labels = 'S308',
                     name = "Water Quality\nGrab Site ID") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = 'roi_conceptual_2.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 7,
  height = 3,
  units = 'in'
)

# Zoomed in ROI and average val determination - should redo this with classified image ----
stack_weekly_chla_tom[[which(dates_weekly_start == '2018-09-23')]]

rad_1435 <- st_buffer(s308c_coords,dist = 1435) #% >% st_transform(crs = wgs84)

ggplot() + 
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = expression(paste('Chlorophyll-',italic(a),' (μg/L)')), 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(20,40,90)) +
  geom_sf(data = st_transform(st_intersection(rad_1435, lake_bound),crs = wgs84), color = 'red', alpha = 0.1) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  guides(fill = guide_colorbar(barheight = 10)) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() + 
  theme(text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

ggsave(
  filename = 'roi_zoomed_out.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5,
  height = 4,
  units = 'in'
)

ggplot() + 
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (µg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(20,40,90)) +
  geom_sf(data = st_transform(st_intersection(rad_1435, lake_bound),crs = wgs84), color = 'red', alpha = 0.1) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  guides(fill = guide_colorbar(barheight = 10)) +
  coord_sf(xlim = c(-80.64,-80.6),ylim = c(26.96,27),expand = F) +
  annotation_scale(location = 'tr', height = unit(0.2,'cm'), width_hint = 0.3) +
  theme_minimal() + theme(text = element_text(size = 12), 
                          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = 'roi_zoomed_in.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5,
  height = 4,
  units = 'in'
)

ggplot() + 
  geom_sf(data = st_transform(st_intersection(rad_1435, lake_bound),crs = wgs84), 
          fill = '#407480', alpha = 1, color = 'red') +
  scale_fill_gradientn(limit = c(0, 45), 
                       na.value = 'transparent', 
                       name = expression(paste('Chlorophyll-',italic(a),' (μg/L)')), 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,30,45)),
                       breaks = c(10,20,30,40)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  guides(fill = guide_colorbar(barheight = 10)) +
  coord_sf(xlim = c(-80.64,-80.6),ylim = c(26.96,27),expand = F) +
  annotation_scale(location = 'tr', height = unit(0.2,'cm'), width_hint = 0.3) +
  theme_minimal() + theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))

# Chl-a raster chla to quintiles ----
ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_scaled[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) + 
  scale_fill_continuous(na.value = 'transparent', name = 'CI multi',
                        type = 'viridis') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  # geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +  
  # labs(title = 'CI multi composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 


ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84)) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(10,20,40,90)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  # geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +
  # labs(title = 'Tomlinson chlorophyll composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  guides(fill = guide_colorbar(barheight = 10)) +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

# Classify and plot
stack_weekly_scaled[[which(year(dates_weekly_start) != 2021)]] -> stack_weekly_scaled_sub

quantile(x = as.vector(values(stack_weekly_scaled_sub)),
         probs = seq(0, 1, .2),
         na.rm = T) -> quintiles_weekly

cut(stack_weekly_scaled[[which(dates_weekly_start == '2016-07-10')]],
    breaks = as.vector(quintiles_weekly)) %>% as.factor() 


colors_5 <- c(colormap(colormaps$salinity,12)[11:10],
              colormap(colormaps$salinity,12)[6:4])

ggplot() +
  layer_spatial(projectRaster(
    cut(stack_weekly_scaled[[which(dates_weekly_start == '2016-07-10')]],
        breaks = as.vector(quintiles_weekly)),    
    crs = wgs84)) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Quintile', 
                       colors = rev(colors_5), 
                       values = rescale(c(1,2,3,4,5)),
                       breaks = c(1,2,3,4,5)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  # geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +
  # labs(title = 'Tomlinson chlorophyll composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  guides(fill = guide_colorbar(barheight = 10)) +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

# Raw CI images ----
par(mfrow = c(1,1))
plot(stack_weekly_unscaled[[which(dates_weekly_start == '2016-07-03')]])

# Scaled image plot
plot(crop(stack_weekly_unscaled[[which(dates_weekly_start == '2016-07-03')]],lake_buffer_extent))
plot(stack_weekly_scaled[[which(dates_weekly_start == mdy('08-12-2018'))]])

plot(crop(stack_weekly_unscaled[[which(dates_weekly_start == '2016-07-10')]],lake_buffer_extent))

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_unscaled_cropped[[which(dates_weekly_start == mdy('08-12-2018'))]],
    crs = wgs84
  )) +
  scale_fill_continuous(na.value = 'transparent', name = 'Digital Number') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_scaled[[which(dates_weekly_start == mdy('08-12-2018'))]],
    crs = wgs84
  )) +
  scale_fill_continuous(na.value = 'transparent', 
                        name = 'Scaled CI value',
                        limits = c(0,round(maxValue(stack_weekly_scaled[[which(dates_weekly_start == mdy('08-12-2018'))]]),digits = 3))) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == mdy('08-12-2018'))]],
    crs = wgs84
  )) +
  scale_fill_continuous(na.value = 'transparent', 
                        name = 'Chl-a (μg/L)',
                        limits = c(-2,127),
                        type = 'viridis') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  theme_minimal() +
  labs(title = 'Tomlinson et al. estimation \nof chlorophyll concentration') +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_seeg[[which(dates_weekly_start == mdy('08-12-2018'))]],
    crs = wgs84
  )) +
  scale_fill_continuous(na.value = 'transparent', 
                        name = 'Chl-a (μg/L)',
                        type = 'viridis') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  theme_minimal() +
  labs(title = 'Seegers estimation of \nchlorophyll concentration') +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

# Maps ----
## SLE Sampling locations - NEW COLORS
ggplot() +
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),fill = 'azure2',alpha = .5) +
  geom_sf(data = st_transform(sle,wgs84)) +
  geom_sf(data = st_transform(sle_station_coords,wgs84), aes(color = station_id), size = 2.5) +
  # coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) + # for whole lake
  coord_sf(xlim = c(-80.7,-80),ylim = c(26.9,27.25)) +
  scale_color_manual(name = 'Sampling \nLocations:',
                    labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                    values = colors) +
  theme_minimal() +
  annotation_scale(location = 'tl', height = unit(0.2,'cm'), width_hint = 0.3) +
  theme(text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.background = element_rect(fill = 'white',color = 'white'),
        panel.grid = element_blank())

ggsave(
   filename = 'sle_map.png',
   path = '/Users/natalie/Documents/thesis_writing/r_plots',
   width = 4.5,
   height = 2,
   units = 'in'
 )

## SLE Sampling locations - ORIGINAL COLORS
ggplot() +
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),fill = 'azure',alpha = .5) +
  geom_sf(data = st_transform(sle,wgs84)) +
  geom_sf(data = st_transform(sle_station_coords,wgs84), aes(color = station_id), size = 2.5) +
  # coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) + # for whole lake
  coord_sf(xlim = c(-80.7,-80),ylim = c(26.9,27.25)) +
  scale_color_manual(name = 'Sampling \nLocations:',
                     labels = c('S308C','C44S80','SE 08B','SE 02','SE 01','SE 11'),
                     values = rev(colormap(colormap = colormaps$YIGnBu,nshades = 12)[4:9])) +
  theme_minimal() +
  annotation_scale(location = 'tl', height = unit(0.2,'cm'), width_hint = 0.3) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))
  
## Lake Okeechobee sampling locations - NEW COLORS
ggplot() + 
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'grey40',size = 2) +
  geom_sf(data = st_transform(s308c_coords,wgs84),color = 'blue',size = 2) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

### Full map -- color
ggplot() +
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  geom_sf(data = st_transform(sle,wgs84),
          color = 'gray30') +
  geom_sf(data = st_transform(caloosahatchee,wgs84),
          color = 'gray30') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(color = 'steelblue4')) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(color = colors[1])) +
  geom_sf(data = st_transform(sle_station_coords[2,],wgs84),size = 2,
          aes(color = colors[2])) +
  geom_sf(data = st_transform(sle_station_coords[3,],wgs84),size = 2,
          aes(color = colors[3])) +
  geom_sf(data = st_transform(sle_station_coords[4,],wgs84),size = 2,
          aes(color = colors[4])) +
  geom_sf(data = st_transform(sle_station_coords[5,],wgs84),size = 2,
          aes(color = colors[5])) +
  geom_sf(data = st_transform(sle_station_coords[6,],wgs84),size = 2,
          aes(color = colors[6])) +
  scale_color_identity(name = 'Sampling station',
                       breaks = c('steelblue4',colors),
                       labels = c('Lake stations',levels(sle_station_coords$station_id)),
                       guide = 'legend') +
  coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) + # for whole lake
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

## Full map - grayscale -- NEED TO FIGURE OUT HOW TO PAIR SCALE_IDENTITY WITH FACTORS SO MY STATIONS ARE IN THE CORRECT ORDER
ggplot() +
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30') +
  geom_sf(data = st_transform(sle,wgs84),
          color = 'gray30') +
  geom_sf(data = st_transform(caloosahatchee,wgs84),
          color = 'gray30') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'gray30') +
  geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
          aes(shape = 1)) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
          aes(shape = shapes[1])) +
  geom_sf(data = st_transform(sle_station_coords[2,],wgs84),size = 2,
          aes(shape = shapes[2])) +
  geom_sf(data = st_transform(sle_station_coords[3,],wgs84),size = 2,
          aes(shape = shapes[3])) +
  geom_sf(data = st_transform(sle_station_coords[4,],wgs84),size = 2,
          aes(shape = shapes[4])) +
  geom_sf(data = st_transform(sle_station_coords[5,],wgs84),size = 2,
          aes(shape = shapes[5])) +
  geom_sf(data = st_transform(sle_station_coords[6,],wgs84),size = 2,
          aes(shape = shapes[6])) +
  scale_shape_identity(name = 'Sampling station',
                       breaks = c(1,shapes),
                       labels = c('Lake stations',
                                  levels(sle_station_coords$station_id)),
                       guide = 'legend') +
  coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) + # for whole lake
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

## Lake Okeechobee sampling locations - ORIGINAL COLORS
ggplot() +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure') +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'grey40',size = 2) +
  geom_sf(data = st_transform(s308c_coords,wgs84),color = '#bbe4b5ff',size = 2) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# PLOTS FOR MANUSCRIPT AND PRESENTATIONS - with ggsave ----
# # Lake O stations without quadrants  ----
# ggplot() +
#   geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
#   geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
#           color = 'gray30') +
#   geom_sf(data = st_transform(sle,wgs84),
#           color = 'gray30') +
#   geom_sf(data = st_transform(caloosahatchee,wgs84),
#           color = 'gray30') +
#   geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
#   geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
#           aes(color = 'steelblue4')) +
#   geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
#           aes(color = 'palegreen')) +
#   scale_color_identity(name = 'Sampling station:',
#                        breaks = c('steelblue4','palegreen'),
#                        labels = c('Lake stations\n(Chlorophyll data)',
#                                   'S-308\n(Chlorophyll and flow data)'),
#                        #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
#                        guide = 'legend') +
#   coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
#   # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
#   theme_minimal() +
#   theme(text = element_text(size = 12),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = 'bottom')
# 
# ggsave(
#   filename = 'lake_stations_plot.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 5.5,
#   height = 5,
#   units = 'in'
# )
# 
# ### Lake Okeechobee with quadrants ----
# ggplot() +
#   geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
#   geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
#           color = 'gray30') +
#   geom_sf(data = st_transform(sle,wgs84),
#           color = 'gray30') +
#   geom_sf(data = st_transform(caloosahatchee,wgs84),
#           color = 'gray30') +
#   geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
#   geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
#           color = 'gray30') +
#   geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
#           color = 'gray30') +
#   geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
#           color = 'gray30') +
#   geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
#           color = 'gray30') +
#   geom_sf(data = st_transform(lake_station_coords,wgs84),size = 2,
#           aes(color = 'steelblue4')) +
#   geom_sf(data = st_transform(s308c_coords,wgs84),size = 2,
#           aes(color = colors[1])) +
#   scale_color_identity(name = 'Sampling station:',
#                        breaks = c('steelblue4',colors[1]),
#                        labels = c('Lake stations\n(Chlorophyll data)',
#                                   'S-308\n(Chlorophyll and flow data)'),
#                        #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
#                        guide = 'legend') +
#   coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
#   # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
#   theme_minimal() +
#   theme(text = element_text(size = 12),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = 'bottom')
# 
# ROI plots for methods section ----

## Create geom for final rad 
rad_b <- st_buffer(s308c_coords,dist = 7000) #% >% st_transform(crs = wgs84)

## Calculate associated roi
st_intersection(rad_b, lake_bound) %>% st_area() %>% as.numeric() -> aoi

## Create geom for initial rad
rad_a <- st_buffer(s308c_coords,dist = sqrt(aoi/pi))

ggplot() +
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30',size = 0.3) +
  geom_sf(data = st_transform(sle,wgs84),
          color = 'gray30') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # ROI HERE
  # geom_sf(data = st_transform(rad_a,crs = wgs84), color = 'red3', fill = 'red3', alpha = 0.5) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2.5,
          color = 'steelblue') +
  coord_sf(xlim = c(-81.12,-80.5),ylim = c(26.68,27.2)) +
  # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), 
                   width_hint = 0.2) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1,'cm'), width = unit(1,'cm'),
                         pad_x = unit(-0.2,'cm'), 
                         pad_y = unit(0.25,'in')
                         ) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave(
  filename = 'roi_methods_1.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 2.25,
  height = 2.25,
  units = 'in'
)

# plot a - uncropped area

# plot_a <- 
  ggplot() +
  geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
          color = 'gray30',size = 0.3) +
  geom_sf(data = st_transform(sle,wgs84),
          color = 'gray30') +
  geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
  # ROI HERE
  geom_sf(data = st_transform(rad_a,crs = wgs84), color = 'red3', fill = 'red3', alpha = 0.5) +
  # geom_sf(data = st_transform(lake_bound,wgs84),
  #         fill = 'transparent',
  #         color = 'gray30') +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 2.5,
          color = 'steelblue') +
  coord_sf(xlim = c(-81.12,-80.5),ylim = c(26.68,27.2)) +
  # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
    annotation_scale(location = 'bl', height = unit(0.2,'cm'), 
                     width_hint = 0.2) +
    annotation_north_arrow(style = north_arrow_minimal,
                           height = unit(1,'cm'), width = unit(1,'cm'),
                           pad_x = unit(-0.2,'cm'), 
                           pad_y = unit(0.25,'in')
    ) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  ggsave(
    filename = 'roi_methods_2.png',
    path = '/Users/natalie/Documents/thesis_writing/r_plots',
    width = 2.5,
    height = 2.5,
    units = 'in'
  )
  
# plot b - scaled roi
# plot_b <- 
  ggplot() +
    geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
    geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
            color = 'gray30',size = 0.3) +
    geom_sf(data = st_transform(sle,wgs84),
            color = 'gray30') +
    geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
    # ROI HERE
    geom_sf(data = st_transform(rad_a,crs = wgs84), color = 'gray40', fill = 'gray40', alpha = 0.5,linetype = 'dashed') +
    geom_sf(data = st_transform(st_intersection(rad_b,lake_bound),crs = wgs84), color = 'red3', fill = 'red3', alpha = 0.5) +
    # geom_sf(data = st_transform(lake_bound,wgs84),
    #         fill = 'transparent',
    #         color = 'gray30') +
    geom_sf(data = st_transform(s308c_coords,wgs84),size = 2.5,
            color = 'steelblue') +
    coord_sf(xlim = c(-81.12,-80.5),ylim = c(26.68,27.2)) +
    # guides(color = guide_legend(override.aes = list(size = c(2,3.5)))) +
    annotation_scale(location = 'bl', height = unit(0.2,'cm'), 
                     width_hint = 0.2) +
    annotation_north_arrow(style = north_arrow_minimal,
                           height = unit(1,'cm'), width = unit(1,'cm'),
                           pad_x = unit(-0.2,'cm'), 
                           pad_y = unit(0.25,'in')
    ) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())

  ggsave(
    filename = 'roi_methods_3.png',
    path = '/Users/natalie/Documents/thesis_writing/r_plots',
    width = 2.5,
    height = 2.5,
    units = 'in'
  )
  
# ggarrange(plot_a,plot_b,
#           labels = c('A','B'),
#           ncol = 2,nrow = 1)

# ggsave(
#   filename = 'roi_methods.png',
#     path = '/Users/natalie/Documents/thesis_writing/r_plots',
#      width = 6.5,
#      height = 3.75,
#      units = 'in'
#    )

# ggarrange( 
#   plot_a + theme(text = element_text(size = 16),
#                  axis.text.x = element_text(angle = 45, hjust = 1)),
#   plot_b + theme(text = element_text(size = 16),
#                  axis.text.x = element_text(angle = 45, hjust = 1)),
#   labels = c('AOI prior to adjusting radius:',
#              'AOI with adjusted ROI:'),
#   ncol = 2,nrow = 1,
#   hjust = -.05,
#   font.label = list(size = 16, color = "black", face = "plain", family = NULL)
# )
# 
# ggsave(
#   filename = 'roi_methods_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 7,
#   height = 4,
#   units = 'in'
# )

# Moved this code to figuring_out_spatial_data script
# ggsave(
#   filename = 'lake_stations_plot_with_quads.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 5.5,
#   height = 5,
#   units = 'in'
# )

# # State map with Lake O, Caloosahatchee, and SLE ----
# ggplot() +
#   geom_sf(data = st_transform(florida_bound,wgs84), alpha = .8) +
#   geom_sf(data = st_transform(lake_bound,wgs84),alpha = .8, fill = 'azure2',
#           color = 'gray30') +
#   geom_sf(data = st_transform(sle,wgs84),
#           color = 'gray30') +
#   geom_sf(data = st_transform(caloosahatchee,wgs84),
#           color = 'gray30') +
#   geom_sf(data = st_transform(st_intersection(littoral_zone,lake_bound),wgs84),alpha = .8, fill = 'azure3',color = 'transparent') +
#    # geom_rect(aes(
#    #   xmin = -81.17,
#    #   xmax = -80.55,
#    #   ymin = 26.63,
#    #   ymax = 27.25
#    # ), fill = 'transparent', color = 'lightcyan4', size = .5) +
#   # coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
#   theme_minimal() +
#   theme(text = element_blank())
# 
# ggsave(
#   filename = 'state_map_with_boundary_box.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 3,
#   height = 3,
#   units = 'in'
# )
# 
# In situ stations over tomlinson image ----
ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84,
    method = 'ngb')) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(20,40,90)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2.5) +
  labs(title = 'Tomlinson chlorophyll composite, \n2016-07-10 to 2016-07-16') +
  theme_minimal() +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 

ggsave(
  filename = 'in_situ_vs_tom.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6,
  height = 6,
  units = 'in'
)

# CI multi example plot ----
ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_scaled[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84,
    method = 'ngb')) + 
  scale_fill_continuous(na.value = 'transparent', name = 'CI_cyano',
                        type = 'viridis') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  # geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +  
  # labs(title = 'CI multi composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 

ggsave(
  filename = 'ci_example_plot.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 4,
  height = 3,
  units = 'in'
)

ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_scaled[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84,
    method = 'ngb')) + 
  scale_fill_continuous(na.value = 'transparent', name = 'CI_cyano',
                        type = 'viridis') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  # geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +  
  # labs(title = 'CI multi composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  guides(fill = guide_colorbar(barheight = 10)) 

ggsave(
  filename = 'ci_example_plot_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6,
  height = 5,
  units = 'in'
)


# Tomlinson example plot ----
ggplot() +
  layer_spatial(projectRaster(
    stack_weekly_chla_tom[[which(dates_weekly_start == '2016-07-10')]],
    crs = wgs84,
    method = 'ngb')) +
  scale_fill_gradientn(limit = c(0, 90), 
                       na.value = 'transparent', 
                       name = 'Chl-a (μg/L)', 
                       colors = c('#440154ff','#31668dff','#37b578ff','#fde725ff'), 
                       values = rescale(c(0,20,40,90)),
                       breaks = c(20,40,90)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  # geom_sf(data = st_transform(lake_station_coords, crs = wgs84),color = 'violetred1', alpha = 0.85,size = 2) +
  # labs(title = 'Tomlinson chlorophyll composite, \n2018-08-12 to 2018-08-18') +
  theme_minimal() +
  guides(fill = guide_colorbar(barheight = 10)) +
  theme(text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

ggsave(
  filename = 'tomlinson_example_plot.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 4,
  height = 3,
  units = 'in'
)

# Classified example plot ----
stack_weekly_scaled[[which(year(dates_weekly_start) != 2021)]] -> stack_weekly_scaled_sub

breaks_weekly <- c(min(minValue(stack_weekly_scaled)),
                   (((40*0.7)-20)/4050),
                   (((40*1.3)-20)/4050),
                   max(maxValue(stack_weekly_scaled)))

mask(stack_classified,
     land_adjacency_flag,
     maskvalue = 0) -> stack_classified

scale_ex_img <- function(raster_image) {
  
  # Create template raster
  x <- raster_image
  scaled_image <- raster(x = extent(x), res = res(x), crs = crs(x), vals = rep(NA, ncell(x)))
  
  # Scale image to CI
  scaled_image <- 10^((3/250)*x - 4.2)
  
  # Remove NA values
  na_sub <- which(x[] == 254 | x[] == 255)
  zero_sub <- which(raster_image[] == 0)
  scaled_image[na_sub] <- NA
  scaled_image[zero_sub] <- 0
  
  return(scaled_image)
}

cut(scale_ex_img(stack_weekly_unscaled[[which(dates_weekly_start == '2016-07-10')]]),
    breaks = as.vector(breaks_weekly)) %>%
  crop(lake_buffer_extent) %>%
  mask(land_adjacency_flag,
       maskvalue = 0) -> class_ex 

colors_3 <- c(colormap(colormaps$salinity,12)[11],
              colormap(colormaps$salinity,12)[9],
              colormap(colormaps$salinity,12)[5])

ggplot() +
  layer_spatial(data = projectRaster(class_ex,    
                                     crs = wgs84,
                                     method = 'ngb')) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Classification', 
                       colors = rev(colors_3), 
                       values = rescale(1:3),
                       breaks = 1:3) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  theme_minimal() +
  guides(fill = guide_colorsteps()) +
  theme(text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

ggsave(
  filename = 'quintile_example_plot.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 4,
  height = 3,
  units = 'in'
)

ggplot() +
  layer_spatial(data = projectRaster(class_ex,    
                                     crs = wgs84,
                                     method = 'ngb')) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Classification', 
                       colors = rev(colors_3), 
                       values = rescale(1:3),
                       breaks = 1:3) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  theme_minimal() +
  guides(fill = guide_colorsteps()) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) 

ggsave(
  filename = 'quintile_example_plot_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6,
  height = 5,
  units = 'in'
)

# Zoomed in ROI plot with quintiles for presentations ----
# Code for legend is in frequency_analysis_ci_cyano; just search "legend" and you should find it 

rad_1502 <- st_buffer(s308c_coords,dist = 1502) #% >% st_transform(crs = wgs84)

ggplot() + 
  layer_spatial(data = projectRaster(class_ex,    
                                     crs = wgs84,
                                     method = 'ngb')) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Quintile', 
                       colors = rev(colors_3), 
                       values = rescale(1:3),
                       breaks = 1:3) +
  guides(fill = guide_colorsteps()) +
  geom_sf(data = st_transform(st_intersection(rad_1502, lake_bound),crs = wgs84), color = 'red', alpha = 0.1) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  guides(fill = guide_colorbar(barheight = 10)) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() + 
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') 

ggsave(
  filename = 'roi_zoomed_out_quints_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6,
  height = 5,
  units = 'in'
)

ggplot() + 
  layer_spatial(data = projectRaster(class_ex,    
                                     crs = wgs84,
                                     method = 'ngb')) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Quintile', 
                       colors = rev(colors_3), 
                       values = rescale(1:3),
                       breaks = 1:3) +
  geom_sf(data = st_transform(st_intersection(rad_1502, lake_bound),crs = wgs84), color = 'red', alpha = 0.1) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  guides(fill = guide_colorsteps()) +
  coord_sf(xlim = c(-80.64,-80.6),ylim = c(26.96,27),expand = F,crs = wgs84) +
  annotation_scale(location = 'tr', height = unit(0.2,'cm'), width_hint = 0.3) +
  theme_minimal() + theme(text = element_text(size = 16), 
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(-80.64,-80.6,.01)) +
  scale_y_continuous(breaks = seq(26.96,27,.01)) 

ggsave(
  filename = 'roi_zoomed_in_quints_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6,
  height = 5,
  units = 'in'
)

ggplot() + 
  geom_sf(data = st_transform(st_intersection(rad_1435, lake_bound),crs = wgs84), 
          fill = colors_5[1], alpha = 1, color = 'red') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  guides(fill = guide_colorbar(barheight = 10)) +
  coord_sf(xlim = c(-80.64,-80.6),ylim = c(26.96,27),expand = F) +
  annotation_scale(location = 'tr', height = unit(0.2,'cm'), width_hint = 0.3) +
  theme_minimal() + theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(-80.64,-80.6,.01)) +
  scale_y_continuous(breaks = seq(26.96,27,.01)) 

ggsave(
  filename = 'roi_zoomed_in_val_quints_pptx.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6,
  height = 5,
  units = 'in'
)

# ROI zoomed plot for manuscript----
zoomed_out <- 
  ggplot() + 
  geom_sf(data = st_transform(st_sym_difference(florida_bound,lake_bound),wgs84), alpha = .8) +
  layer_spatial(data = projectRaster(class_ex,    
                                     crs = wgs84,
                                     method = 'ngb')) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Classification', 
                       colors = rev(colors_3), 
                       values = rescale(c(1,2,3)),
                       breaks = 1:3) +
  guides(fill = guide_colorsteps()) +
  geom_sf(data = st_transform(sle,wgs84)) +
  geom_sf(data = st_transform(st_intersection(rad_1502, lake_bound),crs = wgs84), color = 'red', alpha = 0.1) +
  # coord_sf(xlim = c(-80.64,-80.6),ylim = c(26.96,27),expand = F) + # from zoomed in plot

  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  annotate('rect',
           xmin = -80.65,
           xmax = -80.6,
           ymin = 26.96,
           ymax = 27.01,
           fill = 'transparent',
           color = 'black') +
  guides(fill = guide_colorbar(barheight = 10)) +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_minimal,
                         pad_x = unit(-0.2,'cm'), 
                         pad_y = unit(0.25,'in')) +
  coord_sf(xlim = c(-81.15,-80),ylim = c(26.65,27.25)) + # for whole lake
  # coord_sf(xlim = c(-80.7,-80),ylim = c(26.9,27.25)) +
  # coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() + 
  theme(text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') 

zoomed_in <- 
ggplot() + 
  geom_sf(data = st_transform(st_sym_difference(florida_bound,lake_bound),wgs84), alpha = .8) +
  layer_spatial(data = projectRaster(class_ex,    
                                     crs = wgs84,
                                     method = 'ngb')) +
  scale_fill_gradientn(na.value = 'transparent', 
                       name = 'Quintile', 
                       colors = rev(colors_3), 
                       values = rescale(c(1,2,3)),
                       breaks = 1:3) +
  geom_sf(data = st_transform(sle,wgs84)) +
  geom_sf(data = st_transform(st_intersection(rad_1502, lake_bound),crs = wgs84), color = 'red', alpha = 0.1) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84), color = 'grey30', alpha = 0) +
  geom_sf(data = st_transform(s308c_coords,wgs84),size = 3) +
  # guides(fill = guide_colorsteps()) +
  coord_sf(xlim = c(-80.65,-80.6),ylim = c(26.96,27.01)) + # for whole lake 
  annotation_scale(location = 'tr', height = unit(0.2,'cm'), width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1,'cm'),
                         pad_x = unit(1.32,'in'), pad_y = unit(1.4,'in')) +
  theme_minimal() + 
  theme(text = element_text(size = 12), 
        axis.text = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.background = element_rect(fill = 'white'),
        legend.position = 'none',
        plot.margin=grid::unit(rep(2,4),"pt")) +
  scale_x_continuous(breaks = seq(-80.64,-80.6,.01)) +
  scale_y_continuous(breaks = seq(26.96,27,.01)) 

ggdraw() +
  draw_plot(zoomed_out) +
  draw_plot(zoomed_in, x = 0.5, y = 0.3, width = .6, height = .6)

# ggarrange(zoomed_out,zoomed_in,
#           labels = c('A','B'),
#           ncol = 2,nrow = 1)

ggsave(
  filename = 'roi_zoomed_plot_inset.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5.25,
  height = 3.5,
  units = 'in'
)
