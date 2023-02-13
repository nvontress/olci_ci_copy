# Redoing roi analysis with all images from Blake

# Set up work space ----
# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(sf)
library(janitor)
library(raster)
library(exactextractr)
library(lubridate)
library(beepr)
library(caret)

# Set Working Directory
# setwd('/Volumes/Hard Drive/Data/') # for hard drive
# -- or --
setwd('/Users/natalie/Documents/select_data/') # for data on laptop

# Load work space - last updated 7/8/21 -- goes to modeling (round 2)
# load("~/Documents/R projects/olci_ci/roi_analysis_ci_cyano_workspace.RData")

# Set roi
# roi <- 1000

# Read in data ----

# Loads raster stacks, processed and unprocessed, weekly and daily
load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")

# Set coordinate system
epsg <- crs(stack_daily_unscaled)
wgs84 <- 4326

lake_bound  <- read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200') %>%
  st_transform(epsg)

lake_stations <- c(#'S308C',
  "KISSR0.0","L001","L004","L005","L006","L007","L008",
  #"LZ2",
  "LZ25A","LZ30","LZ40"
  # ,"PALMOUT","POLESOUT"
)

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% lake_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> lake_station_coords

s308c_coords <- read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station == 'S308C',
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg)


#   filter(date(collection_date) >= '2016-04-01')

# Hydrologic data
s308_usgs <- read_csv('NWISdata/20160501-20210430_s308.csv')

# Calculate roi ----

# Calculate total water volume discharged each day
calc_tot_vol <- function(day,flow_tib = s308_usgs) {
  # day <- "2017-08-03"
  # flow_tib <- s308
  
  flow_tib %>%
    filter(date(dateTime) == day,
           !is.na(flow_cfs)) %>%
    mutate(cumulative_volume_discharged = NA) -> day_i_gauge
  
  # Zero out negative flow data
  # flow_tib$flow_cfs[which(flow_tib$flow_cfs < 0)] <- 0
  
  # Calculate the total volume discharged
  day_i_gauge$cumulative_volume_discharged[1] <- 0
  
  if (nrow(day_i_gauge) == 0) {
    return(NA)
  } else {
    for(j in 2:nrow(day_i_gauge)) {
      day_i_gauge$cumulative_volume_discharged[j] <- day_i_gauge$flow_cfs[j]*int_length(int_diff(c(day_i_gauge$dateTime[j-1],day_i_gauge$dateTime[j]))) + day_i_gauge$cumulative_volume_discharged[j-1] 
    }
    
    return(last(day_i_gauge$cumulative_volume_discharged))
  }
}

tibble(
  date = unique(date(s308_usgs$dateTime)),
  tot_vol_ft3 = NA
) -> roi_tibble

pb = txtProgressBar(min = 0, 
                    max = nrow(roi_tibble), 
                    initial = 0)
for(i in 1:nrow(roi_tibble)) {
  roi_tibble$tot_vol_ft3[i] <- calc_tot_vol(roi_tibble$date[i])
  
  setTxtProgressBar(pb,i)
}
beep(2)

# Go ahead and remove NA values
roi_tibble %>%
  filter(!is.na(tot_vol_ft3)) -> roi_tibble

## Take a moment to visualize with and without negative flow
# View(roi_tibble)

# ggplot(data = s308) +
#   geom_point(aes(x = dateTime, y = flow_cfs), alpha = 0.01) + 
#   geom_hline(yintercept = 0,color = 'red') +
#   theme_minimal()
# 
# ggplot(data = roi_tibble) +
#   geom_point(aes(x = date, y = tot_vol_ft3), alpha = .5) +
#   geom_hline(yintercept = 0,color = 'red') +
#   theme_minimal()

# Calculate average height of water over sill
calc_h <- function(day,flow_tib = s308_usgs) {
  flow_tib %>% filter(date(dateTime) == day) %>% pull(gauge_ft) %>% mean(na.rm = T) -> avg_gauge
  return(avg_gauge - 9.1)
}

roi_tibble %>%
  mutate(h_ft = NA) -> roi_tibble

for(i in 1:nrow(roi_tibble)) {
  roi_tibble$h_ft[i] <- calc_h(day = roi_tibble$date[i])
}
beep(2)

# Visualize
# ggplot(data = roi_tibble) +
#   geom_point(aes(x = tot_vol_ft3, y = h_ft), alpha = .6) +
#   theme_minimal()

# No days where the lake height was lower than the sill
# s308 %>% filter(gauge_ft < 9.1)

## How can there be negative flow, yet the gauge height never be lower than the sill? I suppose the downstream gauge is still higher than the sill height

# Calculate area of influence; filter out negative flow days (obviously no export on these days)
roi_tibble %>%
  mutate(aoi_ft2 = tot_vol_ft3/h_ft) %>%
  filter(aoi_ft2 >= 0) -> roi_tibble

## Visualize
# ggplot(data = roi_tibble) +
#   geom_point(aes(x = aoi_ft2, y = h_ft), alpha = .6) +
#   theme_minimal()

# Calculate radius of influence
calc_roi_m <- function(aoi_ft2,
                       lake_geom = lake_bound,
                       sle_outlet = s308c_coords) {
  # aoi_ft2 <- roi_tibble$aoi_ft2[1]
  # lake_geom <- lake_bound
  # sle_outlet <- sle_station_coords[2,]
  
  
  # Convert aoi to m2
  aoi_m2 <- aoi_ft2/(3.28^2)
  
  # Calculate initial radius
  rad_low <-  sqrt(aoi_m2 / pi)
  # Create un upper bound
  rad_high <- rad_low + 10000
  # Calculate median rad
  rad_med <- mean(c(rad_low,rad_high))
  
  iterate_radius <- function(radius) {
    
    st_buffer(sle_outlet,dist = radius) -> circle_m
    st_intersection(circle_m,lake_geom) %>% pull(geometry) %>% st_area() %>% as.numeric() -> est_area_m2
    
    error <- abs(aoi_m2 - est_area_m2)
    
    return(error)
  }
  
  iterate_radius(rad_med)
  
  optim(
    par = rad_med,
    fn = iterate_radius,
    method = 'Brent',
    lower = rad_low,
    upper = rad_high
  ) -> opt
  
  return(opt$par)
  
}

roi_tibble %>% mutate(roi_m = NA) -> roi_tibble
pb = txtProgressBar(min = 0, 
                    max = nrow(roi_tibble), 
                    initial = 0)
for(i in 1:nrow(roi_tibble)) {
  roi_tibble$roi_m[i] <- calc_roi_m(aoi_ft2 = roi_tibble$aoi_ft2[i])
  setTxtProgressBar(pb,i)
}
beep(2)

# Max roi was 1502 m

# Extract values ---- 

# Pull file dates from image file paths

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
beep(2)

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

## Write functions to pull raster images
pull_image_daily<- function(day,
                            image_stack = stack_daily_scaled,
                            image_dates = dates_daily) {
  # day <- roi_tibble$date[60]
  # image_stack <- chla_individual_stack
  tibble(image_dates = ymd(image_dates)) %>%
    mutate(index = row_number()) %>%
    filter(image_dates == day) %>% pull(index) -> index
  
  day_adj <- day
  
  if(length(index) == 0) {
    repeat{
      day_adj <- day_adj - 1
      
      tibble(image_dates = ymd(image_dates)) %>%
        mutate(index = row_number()) %>%
        filter(image_dates == day_adj) %>% pull(index) -> index
      
      image <- image_stack[[index]]
      
      if(length(index) > 0) {break}
    }
  } else {
    image <- image_stack[[index]]
  } 
    
  return(as.list(image))
}
  
pull_image_weekly <- function(day,
                              image_stack = stack_weekly_scaled,
                              image_dates_start = dates_weekly_start,
                              image_dates_end = dates_weekly_end) {

  tibble(image_dates_start = ymd(image_dates_start),
         image_dates_end = ymd(image_dates_end)) %>%
    mutate(index = row_number()) %>%
    filter(image_dates_start <= day,
           image_dates_end >= day) %>% pull(index) -> index
  
  day_adj <- day
  
  if(length(index) == 0) {
    repeat{
      day_adj <- day_adj - 1
      
      tibble(image_dates_start = ymd(image_dates_start),
             image_dates_end = ymd(image_dates_end)) %>%
        mutate(index = row_number()) %>%
        filter(image_dates_start <= day_adj &
               image_dates_end >= day_adj) %>% pull(index) -> index
      
      image <- image_stack[[index]]
      
      if(length(index) > 0) {break}
    }
  } else {
    image <- image_stack[[index]]
  } 
  
  return(as.list(image))
}

pull_image_index_daily <-  function(day,
                                    image_dates = dates_daily) {
  # day <- roi_tibble$date[60]
  # image_stack <- chla_individual_stack
  tibble(image_dates = ymd(image_dates)) %>%
    mutate(index = row_number()) %>%
    filter(image_dates == day) %>% pull(index) -> index
  
  day_adj <- day
  
  if(length(index) == 0) {
    repeat{
      day_adj <- day_adj - 1
      
      tibble(image_dates = ymd(image_dates)) %>%
        mutate(index = row_number()) %>%
        filter(image_dates == day_adj) %>% pull(index) -> index
      
      if(length(index) > 0) {break}
    }
  } else {
    index <- index
  } 
  
  return(index)
}

pull_image_index_weekly <- function(day,
                              image_dates_start = dates_weekly_start,
                              image_dates_end = dates_weekly_end) {
  
  tibble(image_dates_start = ymd(image_dates_start),
         image_dates_end = ymd(image_dates_end)) %>%
    mutate(index = row_number()) %>%
    filter(image_dates_start <= day,
           image_dates_end >= day) %>% pull(index) -> index
  
  day_adj <- day
  
  if(length(index) == 0) {
    repeat{
      day_adj <- day_adj - 1
      
      tibble(image_dates_start = ymd(image_dates_start),
             image_dates_end = ymd(image_dates_end)) %>%
        mutate(index = row_number()) %>%
        filter(image_dates_start <= day_adj &
                 image_dates_end >= day_adj) %>% pull(index) -> index
      
      if(length(index) > 0) {break}
    }
  } else {
    image <- index
  } 
  
  return(index)
}

# Test functions
# pull_image_daily(ymd('2019-06-20'))
# pull_image_weekly(ymd('2019-06-20'))

# roi_tibble %>%
#   filter(date >= dates_daily[1]) %>%
#   mutate(raster_daily = NA,
#          raster_weekly = NA) -> roi_tibble
# 
# # Pull images
# pb = txtProgressBar(min = 0, 
#                     max = nrow(roi_tibble), 
#                     initial = 0)
# for(i in 1:nrow(roi_tibble)) {
#   roi_tibble$raster_daily[i] <- pull_image_daily(day = roi_tibble$date[i])
#   roi_tibble$raster_weekly[i] <- pull_image_weekly(day = roi_tibble$date[i])
#   setTxtProgressBar(pb,i)
# }
roi_tibble %>%
 filter(date >= dates_daily[1]) %>%
 mutate(raster_index_daily = NA,
        raster_index_weekly = NA) -> roi_tibble

pb = txtProgressBar(min = 0, 
                    max = nrow(roi_tibble), 
                    initial = 0)
for(i in 1:nrow(roi_tibble)) {
  roi_tibble$raster_index_daily[i] <- pull_image_index_daily(day = roi_tibble$date[i])
  roi_tibble$raster_index_weekly[i] <- pull_image_index_weekly(day = roi_tibble$date[i])
  
  setTxtProgressBar(pb,i)
}
beep(2)

# Extract ---- SKIP SECTION ----
roi <- round(max(roi_tibble$roi_m)) # NOTE: max roi did not change when I considered negative flow values

# roi_tibble %>% 
  # mutate(extracted_daily = NA,
         # extracted_weekly = NA) -> roi_tibble

pb = txtProgressBar(min = 0, 
                    max = nrow(roi_tibble), 
                    initial = 0)
for(i in 1:nrow(roi_tibble)){
  roi_tibble$extracted_daily[i] <- raster::extract(
    x = roi_tibble$raster_daily[[i]],
    y = sle_station_coords[2,],
    method = 'simple',
    buffer = roi
  )
  
  roi_tibble$extracted_weekly[i] <- raster::extract(
    x = roi_tibble$raster_weekly[[i]],
    y = sle_station_coords[2,],
    method = 'simple',
    buffer = roi
  )
  
  setTxtProgressBar(pb,i)
}
beep(2)

# Calc average ci
roi_tibble %>%
  mutate(average_ci_daily = NA,
         average_ci_weekly = NA) -> roi_tibble

for(i in 1:nrow(roi_tibble)) {
  roi_tibble$average_ci_daily[i] <- mean(roi_tibble$extracted_daily[[i]], na.rm = T)
  roi_tibble$average_ci_weekly[i] <- mean(roi_tibble$extracted_weekly[[i]], na.rm = T)
}

# Extract with constant roi ----

## Set up columns 
exact_extract(
  x = stack_daily_chla_tom[[1]],
  y = st_buffer(sle_station_coords[2,],dist = 1435),
  include_cell = T
) %>%  first() %>% as.tibble() %>% 
  pull(cell) -> cell_list

# Determining the number of cells within the lake that would be extracted
st_buffer(sle_station_coords[2,],dist = 1435) %>% st_intersection(lake_bound)

exact_extract(
  x = stack_daily_chla_tom[[1]],
  y = st_intersection(st_buffer(sle_station_coords[2,],dist = 1435),lake_bound),
  include_cell = T
) %>%  first() %>% as.tibble() %>% 
  pull(cell) %>% length()

exact_extract(
  x = stack_daily_chla_tom[[1]],
  y = st_buffer(sle_station_coords[2,],dist = 1435),
  include_cell = T
) %>%  first() %>% as.tibble() %>% 
  mutate(cell = str_c('cell_',as.character(cell))) %>%
  select(-coverage_fraction) %>%
  spread(cell,value) %>% tibble(.rows = nrow(roi_tibble)) -> cells

exact_extract(
  x = stack_daily_chla_tom[[1]],
  y = st_buffer(sle_station_coords[2,],dist = 1435),
  include_cell = T
) %>%  first() %>% as.tibble() %>% 
  mutate(cell = str_c('coverage_',as.character(cell))) %>%
  select(-coverage_fraction) %>%
  spread(cell,value) %>% tibble(.rows = nrow(roi_tibble)) -> coverage

## Daily ----
roi_tibble %>% 
  mutate(extracted_mean = NA,
         extracted_sum = NA,
         extracted_med = NA) %>% 
  add_column(cells,coverage) -> roi_tibble_daily

## Populate roi_tibble with extracted values
pb = txtProgressBar(min = 0, 
                    max = nrow(roi_tibble_daily), 
                    initial = 0)
for(i in 1:nrow(roi_tibble_daily)){
  # i <- 1
  img_i <- stack_daily_scaled[[roi_tibble_daily$raster_index_daily[i]]]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    fun = 'mean'
  ) -> roi_tibble_daily$extracted_mean[i]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    fun = 'sum'
  ) -> roi_tibble_daily$extracted_sum[i]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    fun = 'median'
  ) -> roi_tibble_daily$extracted_med[i]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    coverage_area = T,
    include_cell = T
  ) %>% first() %>% as_tibble() -> extracted_tib
  
  which(cell_list %in% pull(extracted_tib,cell)) -> cell_subset
  extracted_tib %>% pull(value) -> vals
  extracted_tib %>% pull(coverage_area) -> area
  
  for(j in 1:length(cell_subset)){
    roi_tibble_daily[i,cell_subset[j] + 10] <- vals[j]
    roi_tibble_daily[i,cell_subset[j] + 102] <- area[j]
  }
  
  setTxtProgressBar(pb,i)
}
beep(2)

## Weekly ----
roi_tibble %>% 
  mutate(extracted_mean = NA,
         extracted_sum = NA,
         extracted_med = NA) %>% 
  add_column(cells,coverage) -> roi_tibble_weekly

## Populate roi_tibble with extracted values
pb = txtProgressBar(min = 0, 
                    max = nrow(roi_tibble_weekly), 
                    initial = 0)
for(i in 1:nrow(roi_tibble_weekly)){
  # i <- 1
  img_i <- stack_weekly_scaled[[roi_tibble_weekly$raster_index_weekly[i]]]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    fun = 'mean'
  ) -> roi_tibble_weekly$extracted_mean[i]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    fun = 'sum'
  ) -> roi_tibble_weekly$extracted_sum[i]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    fun = 'median'
  ) -> roi_tibble_weekly$extracted_med[i]
  
  exact_extract(
    x = img_i,
    y = st_buffer(x = sle_station_coords[2,],dist = 1435), #roi_tibble_daily$roi_m[i]),
    coverage_area = T,
    include_cell = T
  ) %>% first() %>% as_tibble() -> extracted_tib
  
  which(cell_list %in% pull(extracted_tib,cell)) -> cell_subset
  extracted_tib %>% pull(value) -> vals
  extracted_tib %>% pull(coverage_area) -> area
  
  for(j in 1:length(cell_subset)){
    roi_tibble_weekly[i,cell_subset[j] + 10] <- vals[j]
    roi_tibble_weekly[i,cell_subset[j] + 102] <- area[j]
  }
  
  setTxtProgressBar(pb,i)
}
beep(2)

## SAVE WORKSPACE
save.image("~/Documents/R projects/olci_ci/roi_analysis_ci_cyano_workspace.RData")

# Heat map ----

# Whole lake heat map

sum(is.na(stack_daily_chla_tom)) %>% plot()

ggplot() +
  layer_spatial(projectRaster(
    sum(is.na(stack_daily_scaled)),
    crs = wgs84)) + 
  scale_fill_viridis_c(option = 'A', na.value = 'transparent', name = 'NA tally') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  labs(title = 'Daily NA heat map, whole lake',
       subtitle = 'Number of layers in stack: 1793') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_colorbar(barheight = 10)) 

ggplot() +
  layer_spatial(projectRaster(
    sum(is.na(stack_weekly_scaled)),
    crs = wgs84)) + 
  scale_fill_viridis_c(option = 'A', na.value = 'transparent', name = 'NA tally') +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  labs(title = 'Weekly NA heat map, whole lake',
       subtitle = 'Number of layers in stack: 263') +  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_colorbar(barheight = 10)) 

# Extracted values heat map
exact_extract(
  x = stack_daily_chla_tom[[1]],
  y = st_buffer(sle_station_coords[2,],dist = 1435),
  include_xy = T
) %>% first() %>% as.data.frame() %>% select(x,y) %>% mutate(na_count = as.vector(colSums(is.na(roi_tibble_daily[,11:102])))) %>% rasterFromXYZ(res = c(300,300),crs = epsg) -> raster_na_daily

exact_extract(
  x = stack_weekly_chla_tom[[1]],
  y = st_buffer(sle_station_coords[2,],dist = 1435),
  include_xy = T,
  include_cell = T
) %>% first() %>% as.data.frame() %>% select(x,y) %>% mutate(na_count = as.vector(colSums(is.na(roi_tibble_weekly[,11:102])))) %>% rasterFromXYZ(res = c(300,300),crs = epsg) -> raster_na_weekly

ggplot() +
  layer_spatial(projectRaster(
  raster_na_daily,
    crs = wgs84)) + 
  scale_fill_viridis_c(option = 'A', na.value = 'transparent', name = 'NA tally', limits = c(500,1000)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(sle_station_coords[2,], crs = wgs84),color = 'violetred1', alpha = 0.85) +  
  geom_sf(data = st_transform(st_buffer(x = sle_station_coords[2,],
                           dist = 1435),crs = wgs84), alpha = 0, color = 'red') +
  labs(title = 'NA heat map, extracted daily images\n(995 total images)') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-80.63902 ,-80.60098),ylim = c(26.96274,26.99725)) +
  guides(fill = guide_colorbar(barheight = 10)) 
  
ggplot() +
  layer_spatial(projectRaster(
    raster_na_weekly,
    crs = wgs84)) + 
  scale_fill_viridis_c(option = 'A', na.value = 'transparent', name = 'NA tally', limits = c(500,1000)) +
  geom_sf(data = st_transform(lake_bound,crs = wgs84),alpha = 0) +
  geom_sf(data = st_transform(sle_station_coords[2,], crs = wgs84),color = 'violetred1', alpha = 0.85) +  
  geom_sf(data = st_transform(st_buffer(x = sle_station_coords[2,],
                                        dist = 1435),crs = wgs84), alpha = 0, color = 'red') +
  labs(title = 'NA heat map, extracted weekly images\n(995 total images)') +
  theme_minimal() +
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_sf(xlim = c(-80.63902 ,-80.60098),ylim = c(26.96274,26.99725)) +
  guides(fill = guide_colorbar(barheight = 10)) 

exact_extract(
  x = stack_weekly_chla_tom[[1]],
  y = st_buffer(sle_station_coords[2,],dist = 1435),
  include_xy = T,
  include_cell = T
) %>% first() %>% as.data.frame() %>% 
  mutate(na_count = as.vector(colSums(is.na(roi_tibble_weekly[,11:102])))) %>% 
  dplyr::select(cell,na_count) %>% as_tibble() -> na_count_tibble

sle_chla %>%
  filter(station_id == 'S308C') %>%
  rename(date = collection_date) %>%
  inner_join(roi_tibble_weekly) %>%
  dplyr::select(date,sample_type_new,value,colnames(roi_tibble_weekly)) -> joined_data_weekly

# which(na_count_tibble$na_count < 700)
# colnames(joined_data_weekly)

colSums(is.na(joined_data_weekly)) %>% as.vector()

tibble(
  cell = as.character(na_count_tibble$cell),
  na_count = as.vector(colSums(is.na(joined_data_weekly[,13:104])))
) %>%
  filter(na_count < 86) %>%
  ggplot() + 
  geom_col(aes(x = cell, y = na_count),alpha = .8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

tibble(
  date = joined_data_weekly$date
) %>%
  add_column(rowSums(is.na(subset(joined_data_weekly[,13:104],select = -which(as.vector(colSums(is.na(joined_data_weekly[,13:104]))) == 86))))) %>% 
  rename(na_count = `rowSums(...)`) %>%
  ggplot() +
  geom_col(aes(as.character(date),na_count), alpha = .8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

(colnames(joined_data_weekly[,13:104]))[which(as.vector(colSums(is.na(joined_data_weekly[,13:104]))) != 86)] %>% str_extract('[:digit:]{5}') -> ci_cells

tibble(
  date = joined_data_weekly$date
) %>%
  add_column(as.tibble(is.na(joined_data_weekly[,13:104]))) %>%
  gather(colnames(joined_data_weekly[,13:104]), 
         key = cell, value = missing_data) %>%
  mutate(cell = str_extract(cell, '[:digit:]{5}'),
         missing_data = factor(x = missing_data,
                               levels = c(F,T),
                               labels = c('Not Missing','Missing'))) %>%
  filter(cell %in% ci_cells) %>%
  ggplot() +
  geom_tile(aes(x = as.character(cell),y = as.character(date),fill = missing_data))  +
  scale_fill_discrete(
    name = 'Fill',
    type = c('lightblue2', 'lightskyblue4')) +
  theme_minimal() +
  xlab('Cell Number') +
  ylab('Date') +
  theme(axis.text.x = element_text(angle = 90))


  
  
is.na(joined_data_weekly[,13:104]) %>% as.tibble()


  
  
  filter()
  ggplot() +
  geom_bar(aes(x = as.character(date), y = `rowSums(...)`))



# Model fitting (Round 2) ----
sle_chla %>%
  filter(station_id == 'S308C') %>%
  rename(date = collection_date) %>%
  inner_join(roi_tibble_daily) %>%
  dplyr::select(date,sample_type_new,value,colnames(roi_tibble_daily)) -> joined_data_daily

sle_chla %>%
  filter(station_id == 'S308C') %>%
  rename(date = collection_date) %>%
  inner_join(roi_tibble_weekly) %>%
  dplyr::select(date,sample_type_new,value,colnames(roi_tibble_weekly)) -> joined_data_weekly

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

# Looks horrendous, but sum is by far the best predictor for something this simple
mod_lm_daily <- lm(value ~ extracted_sum, data = filter(joined_data_daily, 
                                                         !is.na(extracted_mean)))
summary(mod_lm_daily)
par(mfrow = c(2,2))
plot(mod_lm_daily)
ggplot() + 
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point(aes(x = predict(mod_lm_daily),y = dplyr::pull(
    filter(joined_data_daily,!is.na(extracted_mean)),value))) +
  xlim(0,55) +
  ylim(0,55) +
  theme_minimal()

mod_lm_weekly <- lm(value ~ extracted_mean, data = filter(joined_data_weekly, 
                                                        !is.na(extracted_mean)))
summary(mod_lm_weekly)
par(mfrow = c(2,2))
plot(mod_lm_weekly)
ggplot() + 
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point(aes(x = predict(mod_lm_weekly),y = dplyr::pull(
    filter(joined_data_weekly,!is.na(extracted_mean)),value))) +
  xlim(0,55) +
  ylim(0,55) +
  theme_minimal()

# Random forest??
set.seed(99)
rf_daily <- train(
  value ~ .,
  data = joined_data_daily,
  method = 'rf',
  importance = T,
  trControl = trainControl(method = 'none'),
  na.action = na.omit
  )

randomForest::rfImpute(x = dplyr::select(joined_data_daily, -value),
         y = joined_data_daily$value,
         inter = 5,
         ntree = 300)

# Try fitting a model (First attempts) ----
sle_chla %>%
  filter(station_id == 'S308C') %>% 
  rename(date = collection_date) %>%
  inner_join(roi_tibble) %>%
  dplyr::select(date,sample_type_new,value,colnames(roi_tibble)) -> joined_data

joined_data %>% filter(!is.na(average_ci_daily)) %>% filter(!is.na(value)) -> joined_data_daily

mod_daily <- lm(value ~ average_ci_daily, data = joined_data_daily)
summary(mod_daily)
par(mfrow = c(2,2))
plot(mod_daily)
ggplot() + 
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point(aes(x = predict(mod_daily),y = joined_data_daily$value)) +
  theme_minimal()

joined_data %>% filter(!is.na(average_ci_weekly)) %>% filter(!is.na(value)) -> joined_data_weekly

mod_weekly <- lm(sqrt(value) ~ average_ci_weekly, data = joined_data_weekly)
summary(mod_weekly)
par(mfrow = c(2,2))
plot(mod_weekly)
ggplot() + 
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point(aes(x = (predict(mod_weekly))^2,y = joined_data_weekly$value)) +
  xlim(0,155) +
  ylim(0,155) +
  theme_minimal()

ggplot() + 
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_point(aes(x =(joined_data_weekly$average_ci_weekly*coef_seeg + const_seeg),y = joined_data_weekly$value)) +
  xlim(0,155) +
  ylim(0,155) +
  theme_minimal()

### Previous analysis ----
roi_tibble <- tibble(
  image_date = file_dates,
  extracted_vals_unscaled = NA,
  extracted_vals_ci = NA,
  extracted_vals_chla = NA,
  mean_chla = NA,
  daily_discharge = NA
)

for(i in 1:nrow(roi_tibble)){
  roi_tibble$extracted_vals_unscaled[i] <- raster::extract(
    x = stack_unscaled_cropped[[i]],
    y = sle_station_coords[2,],
    method = 'simple',
    buffer = roi
  )
}
beep(2)

# roi_tibble %>% mutate(extracted_vals_ci = NA) -> roi_tibble

for(i in 1:nrow(roi_tibble)){
  unscaled_vals_list_i <- roi_tibble$extracted_vals_unscaled[i] %>% first()
  scaled_vals_list_i <- rep(NA,length(unscaled_vals_list_i))
  
  for(j in 1:length(scaled_vals_list_i)){
    
    if(unscaled_vals_list_i[j] == 0 |
       unscaled_vals_list_i[j] == 254 |
       unscaled_vals_list_i[j] == 255 
    ) {scaled_vals_list_i[j] <- NA} else {
      scaled_vals_list_i[j] <- 10^((3/350)*unscaled_vals_list_i[j] - 4.2)
    }
  }
  
  roi_tibble$extracted_vals_ci[i] <- list(scaled_vals_list_i)
}

# roi_tibble %>% mutate(extracted_vals_chla = NA) -> roi_tibble

for(i in 1:nrow(roi_tibble)){
  # 6620*CI - 3.07(+/-5)
  ci_vals_i <- roi_tibble$extracted_vals_ci[[i]]
  
  chla_vals_i <- 6620*ci_vals_i - 3.07 + 5

  roi_tibble$extracted_vals_chla[i] <- list(chla_vals_i)
}

# Calculate the average chla within the radius
# roi_tibble %>% mutate(mean_chla = NA) -> roi_tibble

for(i in 1:nrow(roi_tibble)){
  roi_tibble$mean_chla[i] <- mean(roi_tibble$extracted_vals_chla[[i]], na.rm = T)
}

# roi_tibble %>% mutate(daily_discharge = NA) -> roi_tibble
for(i in 1:nrow(roi_tibble)){
  s308 %>%
    filter(date(dateTime) == roi_tibble$image_date[i],
           !is.na(flow_cfs)) %>%
    mutate(cumulative_volume_discharged = NA) -> day_i_gauge
  
  if(nrow(day_i_gauge) == 0) {
    roi_tibble$daily_discharge[i] <- NA
  } else {
    # Zero out negative flow values
    for(j in 1:nrow(day_i_gauge)) {
      flow <- day_i_gauge$flow_cfs[j]
      if(flow < 0) {
        day_i_gauge$flow_cfs[j] <- 0
      } else {day_i_gauge$flow_cfs[j] <- flow}
    }
    
    # Calculate the total volume discharged
    day_i_gauge$cumulative_volume_discharged[1] <- 0
    for(j in 2:nrow(day_i_gauge)) {
      day_i_gauge$cumulative_volume_discharged[j] <- day_i_gauge$flow_cfs[j]*int_length(int_diff(c(day_i_gauge$dateTime[j-1],day_i_gauge$dateTime[j]))) + day_i_gauge$cumulative_volume_discharged[j-1] 
    }
    
    last(day_i_gauge$cumulative_volume_discharged) -> roi_tibble$daily_discharge[i]
  }
  
}
beep(2)

# Visualize ----

# Plot remotely sensed values over time
roi_tibble %>%
  dplyr::select(image_date,mean_chla) %>%
  na.omit() %>%
  ggplot() +
  geom_line(aes(x = image_date, y = mean_chla), color = 'blue') +
  geom_point(aes(x = image_date, y = mean_chla)) +
  xlab('Date (April 2016 - Present)') +
  ylab('Estimated chlorophyll concentration within ROI') + 
  theme_minimal()

# Plot remotely sensed and in situ values over time
sle_chla %>%
  filter(station_id == 'S308C') %>%
  rename(image_date = collection_date) %>%
  full_join(roi_tibble) %>% 
  dplyr::select(image_date,value,mean_chla) %>% 
  rename(grab_sample = value,
         avg_satellite_roi = mean_chla,
         date = image_date) %>%
  gather(grab_sample, avg_satellite_roi,
         key = 'sample_type',value = 'value') %>% na.omit() %>%
  ggplot() +
  geom_point(aes(x = date,y = value, color = sample_type), alpha = .5) +
  xlab('Date (April 2016 - Present)') +
  ylab('Chlorophyll value') +
  theme_minimal()


# Plot remotely sensed and in situ values over time
sle_chla %>%
  filter(station_id == 'S308C') %>%
  rename(image_date = collection_date) %>%
  inner_join(na.omit(roi_tibble)) %>% 
  dplyr::select(image_date,value,mean_chla) %>% 
  rename(grab_sample = value,
         avg_satellite_roi = mean_chla,
         date = image_date) %>%
  gather(grab_sample, avg_satellite_roi,
         key = 'sample_type',value = 'value') %>% na.omit() %>%
  ggplot() +
  geom_point(aes(x = date,y = value, color = sample_type)) +
  xlab('Date (April 2016 - Present)') +
  ylab('Chlorophyll value') +
  theme_minimal()

# 1:1 plot of remotely sensed and grab values
sle_chla %>%
  filter(station_id == 'S308C') %>%
  rename(image_date = collection_date) %>%
  left_join(roi_tibble) %>% 
  dplyr::select(image_date,value,mean_chla) %>% na.omit() %>%
  rename(grab_sample = value,
         avg_satellite_roi = mean_chla,
         date = image_date) %>%
  # gather(grab_sample, avg_satellite_roi,
  #        key = 'sample_type',value = 'value') %>%
  ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(grab_sample,avg_satellite_roi)) + 
  xlim(0,35) + 
  ylim(0,35) +
  labs(title = '1:1 plot for concentrations from individual images') +
  xlab('Observed chlorophyll concentration') +
  ylab('Simulated chlorophyll concentration') +
  theme_minimal()

# Flow on x and sensed chla on y
roi_tibble %>%
  full_join(dplyr::select((rename(filter(sle_chla,station_id == 'S308C'),image_date = collection_date)),value,image_date)) %>%
  rename(simulated = mean_chla,
         observed = value) %>%
  gather(simulated,observed,
         key = 'value_type',value = 'chla_value') %>%
  ggplot() + 
  geom_point(aes(x = daily_discharge, y = chla_value, color = value_type), alpha = .75) + 
  scale_color_manual(values = c('paleturquoise3','gray40'), name = 'Value \nclassification') +
  xlab('Daily discharge (ft^3/d)') +
  ylab('Chlorophyll value') +
  theme_minimal()
  
  


# Save work space
# save.image("~/Documents/R projects/olci_ci/roi_analysis_ci_cyano_workspace.RData")

