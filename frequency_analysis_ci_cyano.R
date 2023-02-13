# Classify raster values into quintiles and examine frequency across different quadrants of the lake
# Frequencies were calculated for each classification level by year and month
# Frequency = number of classified pixels across the time range and boundary / (the total number of lake pixels in the boundary * the number of images in the time span) * 100 

# Set up work space ----
# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(sf)
library(janitor)
library(raster)
library(exactextractr)
library(colormap)
library(lubridate)
library(beepr)
library(AnalystHelper)
library(scales)
library(effectsize)

# Set Working Directory
# setwd('/Volumes/Hard Drive/Data/') # for hard drive
# -- or --
setwd('/Users/natalie/Documents/select_data/') # for data on laptop

# Load workspace -- Last saved 08/09/22 to check all calculations
load("~/Documents/R projects/olci_ci_copy/frequency_analysis_ci_cyano_COPY_workspace.RData")

# Read in data ----

# Load raster stacks, processed and unprocessed, weekly and daily
load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")

# Re-crop images
## Daily -- can't ditch the daily analysis because of S3A and S3B drama
stack_daily_unscaled_cropped <- crop(stack_daily_unscaled, lake_buffer_extent)
beep(2)

# Read in land qa flags and apply
land_adjacency_flag <- raster('LandAdjacenyQA.tif') %>% crop(extent(stack_daily_scaled)) 

stack_daily_unscaled_cropped_masked_2 <- crop(stack_daily_unscaled[[2]], lake_buffer_extent) %>%
  mask(land_adjacency_flag,maskvalue = 0) # This image contains no data

## Weekly
stack_weekly_unscaled_cropped <- crop(stack_weekly_unscaled, lake_buffer_extent)
beep(2)

# Load S3A and S3B flyover dates
load("~/Documents/R projects/sentinel_records_workspace.RData")

# Set coordinate system
epsg <- crs(stack_daily_unscaled)
wgs84 <- 4326

lake_bound  <- read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200') %>%
  st_transform(epsg)

## Daily
# As a quick reference -- "scaled" referes to CIcyano values and "unscaled" to Digital Numbers
mask(stack_daily_unscaled_cropped,
     land_adjacency_flag,
     maskvalue = 0) -> stack_daily_unscaled_cropped_masked

mask(stack_daily_scaled,
  land_adjacency_flag,
  maskvalue = 0) -> stack_daily_scaled_masked

## Weekly
mask(stack_weekly_unscaled_cropped,
     land_adjacency_flag,
     maskvalue = 0) -> stack_weekly_unscaled_cropped_masked

mask(stack_weekly_scaled,
     land_adjacency_flag,
     maskvalue = 0) -> stack_weekly_scaled_masked
beep(2)

# Read in lake station coords
lake_stations <- c("KISSR0.0","L001","L004","L005","L006","L007","L008","LZ25A","LZ30","LZ40")

read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% lake_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> lake_station_coords

# ... and S308 coords
s308c_coords <- read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station == 'S308C',
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg)

# Create ROI
roi <- st_buffer(x = s308c_coords,dist = 1502)

# Divide lake into quadrants ----
lake_centroid <- st_centroid(st_transform(lake_bound,wgs84))
lake_extent <- extent(st_transform(lake_bound,wgs84))

xmn <- lake_extent[1]; xmx <- lake_extent[2]; xmd <- lake_centroid$geometry[[1]][1]
ymn <- lake_extent[3]; ymx <- lake_extent[4]; ymd <- lake_centroid$geometry[[1]][2]
p1 <- st_point(c(xmn,ymx)); p2 <- st_point(c(xmd,ymx)); p3 <- st_point(c(xmx,ymx)) 
p4 <- st_point(c(xmn,ymd)); p5 <- st_point(c(xmd,ymd)); p6 <- st_point(c(xmx,ymd)) 
p7 <- st_point(c(xmn,ymn)); p8 <- st_point(c(xmd,ymn)); p9 <- st_point(c(xmx,ymn)) 

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

# Pull image dates from image file paths ----
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
    str_sub(9,12)
  
  yr_prev <- as.character(as.numeric(yr)-1)
  
  origin_date <- str_c(yr_prev,'-12-31') %>% ymd()
  
  str_split(filenames_weekly,pattern = '\\.')[[i]] %>% first() %>% 
    str_extract('[:digit:]{3}$') %>% as.numeric() %>%
    as.Date(origin = origin_date) -> dates_weekly_end[i]
  
  dates_weekly_end[i] - 6 -> dates_weekly_start[i]
}

dates_weekly_start <- as_date(dates_weekly_start)
dates_weekly_end <- as_date(dates_weekly_end)

# Remove images from WY 2016 & WY 2022 ----
# The WY year number is determined by the calendar year in which the WY ends
dates_daily[which(WY(dates_daily) != 2016 & WY(dates_daily) != 2022)] -> dates_daily
stack_daily_scaled_masked[[which(WY(dates_daily) != 2016 & WY(dates_daily) != 2022)]] -> stack_daily_scaled_sub
stack_daily_unscaled_cropped_masked[[which(WY(dates_daily) != 2016 & WY(dates_daily) != 2022)]] -> stack_daily_unscaled_sub

stack_daily_scaled_masked[[which(dates_daily %in% s3a_dates)]] -> stack_daily_scaled_sub_s3a
stack_daily_unscaled_cropped_masked[[which(dates_daily %in% s3a_dates)]] -> stack_daily_unscaled_sub_s3a

# This is different than stack_daily_scaled_sub because it ONLY includes images where ALL of Lake Okeechobee was captured by the flyover event. Because the images were pieced together in the processing by CyAN,  there will be more in stack_daily_scaled_sub because there will be images for days when only a portion of the lake was captured by the satellite.
stack_daily_scaled_masked[[which(dates_daily %in% s3a_s3b_dates)]] -> stack_daily_scaled_sub_s3a_s3b
stack_daily_unscaled_cropped_masked[[which(dates_daily %in% s3a_s3b_dates)]] -> stack_daily_unscaled_sub_s3a_s3b

stack_daily_unscaled_cropped_masked[[which(dates_daily %in% s3b_dates)]] 

dates_weekly_start[which(WY(dates_weekly_start) != 2016 & WY(dates_weekly_start) != 2022)] -> dates_weekly_start
dates_weekly_end[which(WY(dates_weekly_start) != 2016 & WY(dates_weekly_start) != 2022)] -> dates_weekly_end
stack_weekly_scaled_masked[[which(WY(dates_weekly_start) != 2016 & WY(dates_weekly_start) != 2022)]] -> stack_weekly_scaled_sub
stack_weekly_unscaled_cropped_masked[[which(WY(dates_weekly_start) != 2016 & WY(dates_weekly_start) != 2022)]] -> stack_weekly_unscaled_sub

# Determine breaks ----

## Classify daily raster images
# Define breaks
# Original
breaks_daily <- c(min(minValue(stack_daily_scaled_masked),na.rm = T),
                   (((40*0.7)-20)/4050),
                   (((40*1.3)-20)/4050),
                   max(maxValue(stack_daily_scaled_masked),na.rm = T))

# After Naomi's comments -- Tomlinson split, corresponds to "redo_3" tag on figures
# breaks_daily <- c(min(minValue(stack_daily_scaled_masked),na.rm = T),
#                 (40-20-3)/(4050+271), # lower bound of 95% CI
#                 (40-20+3)/(4050-271), # upper bound of 95% CI
#                 max(maxValue(stack_daily_scaled_masked),na.rm = T))

# After Naomi's comments -- Seegers split, corresponds to "redo_4" tag on figures
# breaks_daily <- c(min(minValue(stack_daily_scaled_masked),na.rm = T),
#                 (40 + 3.1 - 5.2)/(6620 + 646), # lower bound of 95% CI
#                 (40 + 3.1 + 5.2)/(6620 - 646), # upper bound of 95% CI
#                 max(maxValue(stack_daily_scaled_masked),na.rm = T))

## Classify weekly raster images
cut(stack_daily_scaled_sub,breaks = breaks_daily) -> stack_daily_classified
names(stack_daily_classified) <- names(stack_daily_unscaled_sub)

cut(stack_daily_scaled_sub_s3a,breaks = breaks_daily) -> stack_daily_classified_s3a
names(stack_daily_classified_s3a) <- names(stack_daily_scaled_sub_s3a)

cut(stack_daily_scaled_sub_s3a_s3b,breaks = breaks_daily) -> stack_daily_classified_s3a_s3b
names(stack_daily_classified_s3a_s3b) <- names(stack_daily_scaled_sub_s3a_s3b)
beep(2)

# Weekly
# Define breaks
# Original
breaks_weekly <- c(min(minValue(stack_weekly_scaled_masked)),
                   (((40*0.7)-20)/4050),
                   (((40*1.3)-20)/4050),
                   max(maxValue(stack_weekly_scaled_masked)))

# After Naomi's comments -- Tomlinson split, corresponds to "redo_3" tag on figures
# breaks_weekly <- c(min(minValue(stack_weekly_scaled_masked)),
#                    (40-20-3)/(4050+271), # lower bound of 95% CI
#                    (40-20+3)/(4050-271), # upper bound of 95% CI
#                    max(maxValue(stack_weekly_scaled_masked)))

# After Naomi's comments -- Seegers split, corresponds to "redo_4" tag on figures
# breaks_weekly <- c(min(minValue(stack_weekly_scaled_masked),na.rm = T),
                  # (40 + 3.1 - 5.2)/(6620 + 646), # lower bound of 95% CI
                  # (40 + 3.1 + 5.2)/(6620 - 646), # upper bound of 95% CI
                  # max(maxValue(stack_weekly_scaled_masked),na.rm = T))

## Classify weekly raster images
cut(stack_weekly_scaled_sub,breaks = breaks_weekly) -> stack_weekly_classified
names(stack_weekly_classified) <- names(stack_weekly_unscaled_sub)
stack_weekly_classified[[1]] # check your work

beep(2)

# Extract Lake Okeechobee values and summarize into a tibble -----

generate_vals_tibble <- function(daily_or_weekly, bound) {
  # daily_or_weekly <- 'weekly'
  # bound <- roi
  
  # Set values based on daily or weekly images
  if(daily_or_weekly == 'daily') {
    stack_classified <- stack_daily_classified
    stack_unscaled_sub <- stack_daily_unscaled_sub
    img_dates <- dates_daily
  } else if(daily_or_weekly == 'daily_s3a') {
    stack_classified <- stack_daily_classified_s3a
    stack_unscaled_sub <- stack_daily_unscaled_sub_s3a
    img_dates <- dates_daily[which(dates_daily %in% s3a_dates)]
  } else if(daily_or_weekly == 'daily_s3a_s3b') {
    stack_classified <- stack_daily_classified_s3a_s3b
    stack_unscaled_sub <- stack_daily_unscaled_sub_s3a_s3b
    img_dates <- dates_daily[which(dates_daily %in% s3a_s3b_dates)]
  } else {
    stack_classified <- stack_weekly_classified
    stack_unscaled_sub <- stack_weekly_unscaled_sub
    img_dates <- dates_weekly_start
  }
  
  # Extract classified values
  exact_extract(stack_classified,bound) %>% 
    first() %>% as_tibble() %>% dplyr::select(-coverage_fraction) %>%
    gather(names(stack_classified),key = 'name',value = 'class') %>%
    group_by(name,class) %>% count() %>% spread(class,n) %>% 
    add_column(tibble(date = img_dates)) %>%
    # Here, NAs are values that weren't classified by the cut function (i.e. DNs of 255 and 254)
    dplyr::select(-`<NA>`) %>% relocate(date,.before = name) -> tib
  
  exact_extract(stack_unscaled_sub,bound) %>%
    first() %>% as_tibble() %>% dplyr::select(-coverage_fraction) %>%
    gather(names(stack_unscaled_sub),key = 'name',value = 'class') %>%
    group_by(name,class) %>% count() %>% spread(class,n) %>% 
    add_column(tibble(date = img_dates)) %>%
    # 0 are a subset of 1 (both indicate no bloom); 255 are land pixels (i.e., no obs possible)
    dplyr::select(date,name,`0`,`255`) %>% relocate(date,.before = name) %>% 
    inner_join(tib,by = 'date') %>% 
    rename(`NA` = `255`) %>% relocate(`NA`,.after = `3`) %>% return()
  
}

## Daily
generate_vals_tibble('daily',lake_bound) -> lake_vals_daily
# beep(2)

## Daily S3A
generate_vals_tibble('daily_s3a',lake_bound) -> lake_vals_daily_s3a
# beep(2)

## Daily S3A and S3B
generate_vals_tibble('daily_s3a_s3b',lake_bound) -> lake_vals_daily_s3a_s3b 
# beep(2)

## Weekly
generate_vals_tibble('weekly',lake_bound) -> lake_vals_weekly
beep(2)

# Visualize classification for entire lake ----

count_pixels <- function(bound) {
  exact_extract(
    stack_daily_unscaled_cropped_masked_2, # a non-flyover date, so all detected water vals are 255
    bound
  ) %>% first() %>% as_tibble() %>% 
    filter(value != 254) %>% nrow() %>% ## Remove land pixels if any not in land mask
    return()
}

# ## Daily
# lake_vals_daily %>%
#   gather(`0`,`1`,`2`,`3`,`4`,`5`,key = class, value = count) %>%
#   ggplot() +
#   geom_tile(aes(x = as.character(date),y = as.character(class),fill = count/count_pixels(lake_bound)*100)) +
#   scale_fill_gradientn(na.value = 'transparent',
#                        colors = colormap(colormap=colormaps$picnic, nshades=8),
#                        name = 'Percent') +
#   theme_minimal() +
#   xlab('Date') + ylab('Classification') +
#   theme(axis.text.x = element_text(angle = 90,size = 4,vjust = 0.5,hjust = 1))
# 
# ## Weekly
# lake_vals_weekly %>%
#   gather(`0`,`1`,`2`,`3`,key = class, value = count) %>%
#   ggplot() +
#   geom_tile(aes(x = as.character(date),y = as.character(class),fill = count/count_pixels(lake_bound)*100)) +
#   scale_fill_gradientn(na.value = 'transparent',
#                        colors = colormap(colormap=colormaps$picnic, nshades=8),
#                        name = 'Percent') +
#   theme_minimal() +
#   xlab('Date') + ylab('Classification') +
#   theme(axis.text.x = element_text(angle = 90,size = 4,vjust = 0.5,hjust = 1))

# Extract quadrant values ----

## Daily - don't need to do S3A and S3B
generate_vals_tibble('daily',nw_quad) %>% 
  mutate(quad = 'Northwest Quadrant') -> nw_vals_daily
generate_vals_tibble('daily',ne_quad) %>% 
  mutate(quad = 'Northeast Quadrant') -> ne_vals_daily
generate_vals_tibble('daily',sw_quad) %>% 
  mutate(quad = 'Southwest Quadrant') -> sw_vals_daily
generate_vals_tibble('daily',se_quad) %>% 
  mutate(quad = 'Southeast Quadrant') -> se_vals_daily

## Weekly
generate_vals_tibble('weekly',nw_quad) %>% 
  mutate(quad = 'Northwest Quadrant') -> nw_vals_weekly
generate_vals_tibble('weekly',ne_quad) %>% 
  mutate(quad = 'Northeast Quadrant') -> ne_vals_weekly
generate_vals_tibble('weekly',sw_quad) %>% 
  mutate(quad = 'Southwest Quadrant') -> sw_vals_weekly
generate_vals_tibble('weekly',se_quad) %>% 
  mutate(quad = 'Southeast Quadrant') -> se_vals_weekly

# Combine and visualize quadrant values----

## WILL HAVE TO REVISIT TO FIGURE OUT HOW TO DO PERCENTAGES
## Daily
# bind_rows(nw_vals_daily,ne_vals_daily,sw_vals_daily,se_vals_daily) %>%
#   gather(`0`,`1`,`2`,`3`,key = class, value = count) %>%
#   ggplot() +
#   geom_tile(aes(x = as.character(date),y = as.character(class),fill = count)) +
#   scale_fill_gradientn(na.value = 'transparent',
#                        colors = colormap(colormap=colormaps$picnic, nshades=8),
#                        name = 'Percent') +
#   facet_wrap(~quad) +
#   theme_minimal() +
#   xlab('Date') + ylab('Classification') +
#   theme(axis.text.x = element_text(angle = 90,size = 4,vjust = 0.5,hjust = 1))

## Weekly
# bind_rows(nw_vals_weekly,ne_vals_weekly,sw_vals_weekly,se_vals_weekly) %>%
#   gather(`0`,`1`,`2`,`3`,key = class, value = count) %>%
#   ggplot() +
#   geom_tile(aes(x = as.character(date),y = as.character(class),fill = count)) +
#   scale_fill_gradientn(na.value = 'transparent',
#                        colors = colormap(colormap=colormaps$picnic, nshades=8),
#                        name = 'Pixel count') +
#   facet_wrap(~quad) +
#   theme_minimal() +
#   xlab('Date') + ylab('Classification') +
#   theme(axis.text.x = element_text(angle = 90,size = 4,vjust = 0.5,hjust = 1))

# Extract ROI values ----

## Daily
generate_vals_tibble('daily',roi) -> roi_vals_daily 

## Weekly
generate_vals_tibble('weekly',roi) -> roi_vals_weekly
beep(2)

# Visualize ROI values ----

## Daily
# roi_vals_daily %>%
#   gather(`0`,`1`,`2`,`3`,key = class, value = count) %>%
#   ggplot() +
#   geom_tile(aes(x = as.character(date),y = as.character(class),fill = count/count_pixels(roi)*100)) +
#   scale_fill_gradientn(na.value = 'transparent',
#                        colors = colormap(colormap=colormaps$picnic, nshades=8),
#                        name = 'Pixel count') +
#   theme_minimal() +
#   xlab('Date') + ylab('Classification') +
#   theme(axis.text.x = element_text(angle = 90,size = 4,vjust = 0.5,hjust = 1))

## Weekly
# roi_vals_weekly %>%
#   gather(`0`,`1`,`2`,`3`,key = class, value = count) %>%
#   ggplot() +
#   geom_tile(aes(x = as.character(date),y = as.character(class),fill = count)) +
#   scale_fill_gradientn(na.value = 'transparent',
#                        colors = colormap(colormap=colormaps$picnic, nshades=8),
#                        name = 'Pixel count') +
#   theme_minimal() +
#   xlab('Date') + ylab('Classification') +
#   theme(axis.text.x = element_text(angle = 90,size = 4,vjust = 0.5,hjust = 1))

# Calculate daily frequencies ----

## Annual ----
annual_freq <- function(tib,bound){
  # tib <- lake_vals_daily
  # bound <- lake_bound
  
  tib %>% group_by(WY(date)) %>% 
    summarize(
      total_1 = sum(`0`,na.rm = T) + sum(`1`,na.rm = T), # Combined "below limit" and "no bloom"
      total_2 = sum(`2`,na.rm = T),
      total_3 = sum(`3`,na.rm = T),
      total_NA = sum(`NA`,na.rm = T),
      n = sum(!is.na(name.x)) # Number of images
    ) %>% mutate(
      freq_1 = total_1/((n*count_pixels(bound))-total_NA)*100,
      freq_2 = total_2/((n*count_pixels(bound))-total_NA)*100,
      freq_3 = total_3/((n*count_pixels(bound))-total_NA)*100,
    ) %>% return()
}

annual_freq_lake_daily <- annual_freq(tib = lake_vals_daily,bound = lake_bound)
annual_freq_nw_quad_daily <- annual_freq(tib = nw_vals_daily,bound = nw_quad)
annual_freq_ne_quad_daily <- annual_freq(tib = ne_vals_daily,bound = ne_quad)
annual_freq_sw_quad_daily <- annual_freq(tib = sw_vals_daily,bound = sw_quad)
annual_freq_se_quad_daily <- annual_freq(tib = se_vals_daily,bound = se_quad)
annual_freq_roi_daily <- annual_freq(tib = roi_vals_daily,bound = roi)

annual_freq_lake_daily_s3a <- annual_freq(tib = lake_vals_daily_s3a, bound = lake_bound)
annual_freq_lake_daily_s3a_s3b <- annual_freq(tib = lake_vals_daily_s3a_s3b, bound = lake_bound)

# S3A/3B comparison ----
cor.test(x = (annual_freq_lake_daily_s3a$freq_2 + annual_freq_lake_daily_s3a$freq_3),
         y = (annual_freq_lake_daily_s3a_s3b$freq_2 + annual_freq_lake_daily_s3a_s3b$freq_3),
         method = 'kendall',
         alternative = 'two.sided')

wilcox.test(x = (annual_freq_lake_daily_s3a$freq_2 + annual_freq_lake_daily_s3a$freq_3),
            y = (annual_freq_lake_daily_s3a_s3b$freq_2 + annual_freq_lake_daily_s3a_s3b$freq_3),
            alternative = 'two.sided',
            paired = T)

cohens_d(x = (annual_freq_lake_daily_s3a$freq_2 + annual_freq_lake_daily_s3a$freq_3),
              y = (annual_freq_lake_daily_s3a_s3b$freq_2 + annual_freq_lake_daily_s3a_s3b$freq_3),
              alternative = 'two.sided',
              paired = T)

# Table 2 in manuscript
tibble(WY = seq(2017,2021,1),
       S3A = (annual_freq_lake_daily_s3a$freq_2 + annual_freq_lake_daily_s3a$freq_3),
       S3AB = (annual_freq_lake_daily_s3a_s3b$freq_2 + annual_freq_lake_daily_s3a_s3b$freq_3),
       diff = S3AB - S3A
) 

beep(2)
# 
# wilcox_effsize(data = df,
#                formula = x ~ y,
#                paired = T) # I don't understand


### Plot ----
# scales::show_col(
#   colormap(
#     colormap = colormaps$salinity,
#     nshades = 12
#   )
# )

# scales::show_col(colormap(colormaps$salinity,12))
# 
colors_3 <- c(colormap(colormaps$salinity,12)[11],
              colormap(colormaps$salinity,12)[9],
              colormap(colormaps$salinity,12)[5])
# 
# colors_4 <- c(colormap(colormaps$salinity,12)[11],
#               colormap(colormaps$salinity,12)[9],
#               colormap(colormaps$salinity,12)[5],
#               colormap(colormaps$salinity,12)[1])
# 
# colors_5 <- c(colormap(colormaps$salinity,12)[11:10],
#               colormap(colormaps$salinity,12)[6:4])
# 
# colors_6 <- c(colormap(colormaps$salinity,12)[11:10],
#               colormap(colormaps$salinity,12)[6:4],
#               colormap(colormaps$salinity,12)[1])
# 
plot_annual_freq <- function(annual_freq_tib){
  annual_freq_tib %>%
    dplyr::select(-total_1,-total_2,-total_3,-n) %>%
    rename(
      Year = `WY(date)`,
      `1` = freq_1,
      `2` = freq_2,
      `3` = freq_3,
      # `4` = freq_4,
      # `5` = freq_5
    ) %>%
    mutate(Year = factor(Year,levels = unique(Year),ordered = T)) %>%
    gather(`1`,`2`,`3`,key = 'Class',value = 'Frequency') %>%
    mutate(Class = factor(Class,levels = c(3,2,1),labels = c('Bloom','Bloom Possible','No Bloom'), ordered = T)) %>%
    ggplot(aes(x = Year,y = Frequency,fill = Class)) +
    geom_bar(position = 'stack',stat = 'identity') +
    scale_fill_manual(values = colors_3,name = 'Remote\nsensing\nclassifications') +
    xlab('Water Year') +
    scale_y_continuous(breaks = seq(0,100,10)) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(margin = margin(-10,0,5,0)),
          axis.text.y = element_text(margin = margin(0,0,0,5)),
          panel.grid.major.x = element_blank()) %>% return()
}
# 
# plot_annual_freq_pptx <- function(annual_freq_tib){
#   annual_freq_tib %>% 
#     dplyr::select(-total_1,-total_2,-total_3,-n) %>%
#     rename(
#       Year = `year(date)`,
#       `1` = freq_1,
#       `2` = freq_2,
#       `3` = freq_3,
#       `4` = freq_4,
#       `5` = freq_5
#     ) %>%
#     mutate(Year = factor(Year,levels = unique(Year),ordered = T)) %>%
#     gather(`1`,`2`,`3`,`4`,`5`,key = 'Class',value = 'Frequency') %>%
#     mutate(Class = factor(Class,levels = c(3,2,1),labels = c('Bloom','Bloom Possible','No Bloom'), ordered = T)) %>%
#     ggplot(aes(x = Year,y = Frequency,fill = Class)) +
#     geom_bar(position = 'stack',stat = 'identity') +
#     scale_fill_manual(values = colors_3,name = 'Remote\nsensing\nclassifications') +
#     scale_y_continuous(breaks = seq(0,100,10)) +
#     theme_minimal() +
#     theme(text = element_text(size = 16),
#           axis.text.x = element_text(margin = margin(-10,0,5,0)),
#           axis.text.y = element_text(margin = margin(0,0,0,5)),
#           panel.grid.major.x = element_blank()) %>% return()
# }
# 
# plot_annual_freq_with_0 <- function(annual_freq_tib){
#   annual_freq_tib %>% 
#     dplyr::select(-total_1,-total_2,-total_3,-n) %>%
#     rename(
#       Year = `WY(date)`,
#       `0` = freq_0,
#       `1` = freq_1,
#       `2` = freq_2,
#       `3` = freq_3
#     ) %>%
#     mutate(Year = factor(Year,levels = unique(Year),ordered = T)) %>%
#     gather(`0`,`1`,`2`,`3`, key = 'Class',value = 'Frequency') %>%
#     mutate(Class = factor(Class,levels = c('3','2','1','0'), 
#                           labels = c('Bloom','Bloom Possible','No Bloom','Below \nCI_cyano'),ordered = T)) %>%
#     ggplot(aes(x = Year,y = Frequency,fill = Class)) +
#     geom_bar(position = 'stack',stat = 'identity') +
#     scale_fill_manual(values = colors_4,name = 'Quintile') +
#     scale_y_continuous(breaks = seq(0,100,10)) +
#     theme_minimal()+
#     theme(text = element_text(size = 12),
#           axis.text.x = element_text(margin = margin(-10,0,5,0)),
#           axis.text.y = element_text(margin = margin(0,0,0,5)),
#           panel.grid.major.x = element_blank()) %>% return()
# }

library(ggpubr)
ggarrange(plot_annual_freq(annual_freq_lake_daily_s3a),
          plot_annual_freq(annual_freq_lake_daily_s3a_s3b),
          labels = c('A. Annual frequencies from \nall daily CI images',
                     'B. Annual frequencies from \ndays when S3A flew over'),
          hjust = 0, vjust = 1,
          common.legend = T,
          legend = 'right')

# plot_annual_freq(annual_freq_lake_daily)

# ggsave(
#   filename = 'use_for_quint_legend_325in_height_word_95.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 6.5,
#   height = 3.25,
#   units = 'in'
# )

# ggsave(
#   filename = 'use_for_quint_legend_5x9_word_height.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 9,
#   height = 5,
#   units = 'in'
# )

# plot_annual_freq_pptx(annual_freq_lake_daily)

# ggsave(
#   filename = 'use_for_quint_legend_5x6in_hxw_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 6,
#   height = 5,
#   units = 'in'
# )

# ggsave(
# filename = 'use_for_quint_legend_5x9_hxw_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#    width = 9,
#    height = 5,
#    units = 'in'
#  )

# ggsave(
#   filename = 'use_for_quint_legend_55x10_hxw_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 10,
#   height = 5.5,
#   units = 'in'
# )

# plot_annual_freq_with_0(annual_freq_lake_daily)
# 
# plot_annual_freq(annual_freq_roi_daily)
# plot_annual_freq_with_0(annual_freq_roi_daily)
# 
# bind_rows(mutate(annual_freq_nw_quad_daily,quad = 'Northwest Quadrant'),
#           mutate(annual_freq_ne_quad_daily,quad = 'Northeast Quadrant'),
#           mutate(annual_freq_sw_quad_daily,quad = 'Southwest Quadrant'),
#           mutate(annual_freq_se_quad_daily,quad = 'Southeast Quadrant')) %>%
#   plot_annual_freq() + facet_wrap(~quad)
# 
# bind_rows(mutate(annual_freq_nw_quad_daily,quad = 'Northwest Quadrant'),
#           mutate(annual_freq_ne_quad_daily,quad = 'Northeast Quadrant'),
#           mutate(annual_freq_sw_quad_daily,quad = 'Southwest Quadrant'),
#           mutate(annual_freq_se_quad_daily,quad = 'Southeast Quadrant')) %>%
#   plot_annual_freq_with_0() + facet_wrap(~quad)

## Monthly ----
monthly_freq <- function(tib,bound){
  tib %>% group_by(month(date)) %>% 
    summarize(
      total_1 = sum(`0`,na.rm = T) + sum(`1`,na.rm = T), # Combined "below limit" and "no bloom"
      total_2 = sum(`2`,na.rm = T),
      total_3 = sum(`3`,na.rm = T),
      total_NA = sum(`NA`,na.rm = T),
      n = sum(!is.na(name.x))
    ) %>% mutate(
      freq_1 = total_1/((n*count_pixels(bound))-total_NA)*100,
      freq_2 = total_2/((n*count_pixels(bound))-total_NA)*100,
      freq_3 = total_3/((n*count_pixels(bound))-total_NA)*100
    ) %>% return()
}

monthly_freq_lake_daily <- monthly_freq(tib = lake_vals_daily,lake_bound)
monthly_freq_nw_quad_daily <- monthly_freq(tib = nw_vals_daily,nw_quad)
monthly_freq_ne_quad_daily <- monthly_freq(tib = ne_vals_daily,ne_quad)
monthly_freq_sw_quad_daily <- monthly_freq(tib = sw_vals_daily,sw_quad)
monthly_freq_se_quad_daily <- monthly_freq(tib = se_vals_daily,se_quad)
monthly_freq_roi_daily <- monthly_freq(tib = roi_vals_daily,roi)

### Plot ----
# plot_monthly_freq <- function(monthly_freq_tib){
#   monthly_freq_tib %>% 
#     dplyr::select(-total_1,-total_2,-total_3,-n) %>%
#     rename(
#       Month = `month(date)`,
#       `1` = freq_1,
#       `2` = freq_2,
#       `3` = freq_3
#     ) %>%
#     mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
#     gather(`1`,`2`,`3`,key = 'Class',value = 'Frequency') %>%
#     mutate(Class = factor(Class,levels = c(3,2,1), ordered = T)) %>%
#     ggplot(aes(x = Month,y = Frequency,fill = Class)) +
#     geom_bar(position = 'stack',stat = 'identity') +
#     scale_fill_manual(values = colors_3,name = 'Quintile') +
#     scale_y_continuous(breaks = seq(0,100,10)) +
#     theme_minimal() + 
#     theme(text = element_text(size = 14),
#           axis.text.x = element_text(angle = 45,
#                                      margin = margin(-5,0,0,0)),
#           axis.text.y = element_text(margin = margin(0,0,0,5)),
#           panel.grid.major.x = element_blank()) %>% return()
# }
# 
# plot_monthly_freq_with_0 <- function(monthly_freq_tib){
#   monthly_freq_tib %>% 
#     dplyr::select(-total_1,-total_2,-total_3,-n) %>%
#     rename(
#       Month = `month(date)`,
#       `0` = freq_0,
#       `1` = freq_1,
#       `2` = freq_2,
#       `3` = freq_3
#     ) %>%
#     mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
#     gather(`0`,`1`,`2`,`3`,key = 'Class',value = 'Frequency') %>%
#     mutate(Class = factor(Class,levels = c(3,2,1,0), 
#                           labels = c('3','2','1','Below \nCI_cyano'), ordered = T)) %>%
#     ggplot(aes(x = Month,y = Frequency,fill = Class)) +
#     geom_bar(position = 'stack',stat = 'identity') +
#     scale_fill_manual(values = colors_4,name = 'Quintile') +
#     scale_y_continuous(breaks = seq(0,100,10)) +
#     theme_minimal() + 
#     theme(text = element_text(size = 14),
#           axis.text.x = element_text(angle = 45,
#                                      margin = margin(-5,0,0,0)),
#           axis.text.y = element_text(margin = margin(0,0,0,5)),
#           panel.grid.major.x = element_blank()) %>% return()
# }
# 
# plot_monthly_freq(monthly_freq_lake_daily)
# plot_monthly_freq_with_0(monthly_freq_lake_daily)
# 
# plot_monthly_freq(monthly_freq_roi_daily)
# plot_monthly_freq_with_0(monthly_freq_roi_daily)
# 
# bind_rows(mutate(monthly_freq_nw_quad_daily,quad = 'Northwest Quadrant'),
#           mutate(monthly_freq_ne_quad_daily,quad = 'Northeast Quadrant'),
#           mutate(monthly_freq_sw_quad_daily,quad = 'Southwest Quadrant'),
#           mutate(monthly_freq_se_quad_daily,quad = 'Southeast Quadrant')) %>%
#   plot_monthly_freq() + facet_wrap(~quad)
# 
# bind_rows(mutate(monthly_freq_nw_quad_daily,quad = 'Northwest Quadrant'),
#           mutate(monthly_freq_ne_quad_daily,quad = 'Northeast Quadrant'),
#           mutate(monthly_freq_sw_quad_daily,quad = 'Southwest Quadrant'),
#           mutate(monthly_freq_se_quad_daily,quad = 'Southeast Quadrant')) %>%
#   plot_monthly_freq_with_0() + facet_wrap(~quad)

# Calculate weekly frequencies ----

## Annual ----

annual_freq_lake_weekly <- annual_freq(tib = lake_vals_weekly,lake_bound)
annual_freq_nw_quad_weekly <- annual_freq(tib = nw_vals_weekly,nw_quad)
annual_freq_ne_quad_weekly <- annual_freq(tib = ne_vals_weekly,ne_quad)
annual_freq_sw_quad_weekly <- annual_freq(tib = sw_vals_weekly,sw_quad)
annual_freq_se_quad_weekly <- annual_freq(tib = se_vals_weekly,se_quad)
annual_freq_roi_weekly <- annual_freq(tib = roi_vals_weekly,roi)

### Plot ----

# plot_annual_freq(annual_freq_lake_weekly) 
# plot_annual_freq_with_0(annual_freq_lake_weekly)
# 
# plot_annual_freq(annual_freq_roi_weekly)
# plot_annual_freq_with_0(annual_freq_roi_weekly)

# bind_rows(mutate(annual_freq_nw_quad_weekly,quad = 'Northwest Quadrant'),
#           mutate(annual_freq_ne_quad_weekly,quad = 'Northeast Quadrant'),
#           mutate(annual_freq_sw_quad_weekly,quad = 'Southwest Quadrant'),
#           mutate(annual_freq_se_quad_weekly,quad = 'Southeast Quadrant')) %>%
#   plot_annual_freq() + facet_wrap(~quad) +
#   theme(axis.text.x = element_text(margin = margin(0,0,0,0)))

# bind_rows(mutate(annual_freq_nw_quad_weekly,quad = 'Northwest Quadrant'),
#           mutate(annual_freq_ne_quad_weekly,quad = 'Northeast Quadrant'),
#           mutate(annual_freq_sw_quad_weekly,quad = 'Southwest Quadrant'),
#           mutate(annual_freq_se_quad_weekly,quad = 'Southeast Quadrant')) %>%
#   plot_annual_freq_with_0() + facet_wrap(~quad) +
#   theme(axis.text.x = element_text(margin = margin(0,0,0,0)))

## Monthly ----

monthly_freq_lake_weekly <- monthly_freq(tib = lake_vals_weekly,bound = lake_bound)
monthly_freq_nw_quad_weekly <- monthly_freq(tib = nw_vals_weekly,bound = nw_quad)
monthly_freq_ne_quad_weekly <- monthly_freq(tib = ne_vals_weekly,bound = ne_quad)
monthly_freq_sw_quad_weekly <- monthly_freq(tib = sw_vals_weekly,bound = sw_quad)
monthly_freq_se_quad_weekly <- monthly_freq(tib = se_vals_weekly,bound = se_quad)
monthly_freq_roi_weekly <- monthly_freq(tib = roi_vals_weekly,bound = roi)

### Plot ----

# plot_monthly_freq(monthly_freq_lake_weekly) 
# plot_monthly_freq_with_0(monthly_freq_lake_weekly)

# plot_monthly_freq(monthly_freq_roi_weekly)
# plot_monthly_freq_with_0(monthly_freq_roi_weekly)

# bind_rows(mutate(monthly_freq_nw_quad_weekly,quad = 'Northwest Quadrant'),
#           mutate(monthly_freq_ne_quad_weekly,quad = 'Northeast Quadrant'),
#           mutate(monthly_freq_sw_quad_weekly,quad = 'Southwest Quadrant'),
#           mutate(monthly_freq_se_quad_weekly,quad = 'Southeast Quadrant')) %>%
#   plot_monthly_freq() + facet_wrap(~quad) +
#   theme(axis.text.x = element_text(margin = margin(0,0,0,0)))

# bind_rows(mutate(monthly_freq_nw_quad_weekly,quad = 'Northwest Quadrant'),
#           mutate(monthly_freq_ne_quad_weekly,quad = 'Northeast Quadrant'),
#           mutate(monthly_freq_sw_quad_weekly,quad = 'Southwest Quadrant'),
#           mutate(monthly_freq_se_quad_weekly,quad = 'Southeast Quadrant')) %>%
#   plot_monthly_freq_with_0() + facet_wrap(~quad) +
#   theme(axis.text.x = element_text(margin = margin(0,0,0,0)))

# Tile tibble ----

# Used for hydrograp in frequency_analysis_in_situ.R
generate_tile_tibble <- function(tib,bound){
  # tib <- lake_vals_weekly
  # bound <- lake_bound
  
  tib %>% group_by(month(date),WY(date)) %>% 
    rename(month = `month(date)`,
           year = `WY(date)`) %>%
    arrange(year,month) %>%
    summarise(total_1 = sum(`1`,na.rm = T) + sum(`0`,na.rm = T), # combine no bloom and below algorithm limit
              total_2 = sum(`2`,na.rm = T),
              total_3 = sum(`3`,na.rm = T),
              total_NA = sum(`NA`,na.rm = T),
              n = sum(!is.na(name.x))) %>% # number of images in that month
    mutate(
      total_pot_bloom_freq = (total_2 + total_3)/((n*count_pixels(bound))-total_NA)*100
      ) %>% return()

}

# lake_tile_vals <- generate_tile_tibble(lake_vals_weekly,lake_bound)
roi_tile_vals <- generate_tile_tibble(roi_vals_weekly,roi)
# nw_tile_vals <- generate_tile_tibble(nw_vals_weekly,nw_quad)
# ne_tile_vals <- generate_tile_tibble(ne_vals_weekly,ne_quad)
# sw_tile_vals <- generate_tile_tibble(sw_vals_weekly,sw_quad)
# se_tile_vals <- generate_tile_tibble(se_vals_weekly,se_quad)
# beep(2)

# roi_tile_vals %>% 
#   mutate(total_bloom = total_2 + total_3,
#          bloom_freq = total_bloom/(count_pixels(roi)*n)*100) %>% 
#   rename(WY = `WY(date)`,
#          month = `month(date)`) %>%
#   arrange(WY, month) %>%
#   ggplot() +
#   geom_tile(aes(x = month,y = WY,fill = bloom_freq))

# Heat maps for WY 2019 and July ----

## Prep data ##

stack_classified_WY2019 <- stack_weekly_classified[[which(WY(dates_weekly_start) == 2019)]] 
stack_count_WY2019 <- stack_classified_WY2019
for(i in 1:nlayers(stack_classified_WY2019)) {
  stack_count_WY2019[[i]][stack_count_WY2019[[i]] < 2] <- 0
  stack_count_WY2019[[i]][stack_count_WY2019[[i]] >=  2] <- 1
}
calc(x = stack_count_WY2019,fun = sum, na.rm = T)/nlayers(stack_count_WY2019)*100 -> heat_map_WY2019
heat_map_WY2019[heat_map_WY2019 == 0] <- NA

stack_classified_july <- stack_weekly_classified[[which(month(dates_weekly_start) == 7)]] 
stack_count_july <- stack_classified_july
for(i in 1:nlayers(stack_classified_july)) {
  stack_count_july[[i]][stack_count_july[[i]] < 2] <- 0
  stack_count_july[[i]][stack_count_july[[i]] >=  2] <- 1
}
calc(x = stack_count_july,fun = sum, na.rm = T)/nlayers(stack_count_july)*100 -> heat_map_july
heat_map_july[heat_map_july == 0] <- NA

## Plot ##
florida_bound <- read_sf('cb_2019_us_state_500k/cb_2019_us_state_500k.shp') %>%
  filter(NAME == 'Florida') %>% 
  st_transform(crs = epsg)

WY2019_plot <- ggplot() +
  layer_spatial(data = crop(projectRaster(heat_map_WY2019,    
                                    crs = wgs84,
                                    method = 'ngb'),st_transform(lake_bound,wgs84))) +
  scale_fill_viridis_c(na.value = 'transparent', 
                       name = 'Potential bloom\nfrequencies:\n ', 
                       option = 'D') +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'white') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'white') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'white') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'white')  +
  geom_sf(data = st_transform(st_sym_difference(florida_bound,lake_bound),wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),color = 'black', fill = 'transparent') +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'white', shape = 1) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1,'cm'),
                         pad_x = unit(-0.15,'in'),
                         pad_y = unit(0.3,'in')) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

july_plot <- ggplot() +
  # geom_sf(data = st_transform(lake_bound,wgs84)) +
  layer_spatial(data = crop(projectRaster(heat_map_july,    
                                     crs = wgs84,
                                     method = 'bilinear'),st_transform(lake_bound,wgs84))) +
  scale_fill_viridis_c(na.value = 'transparent', 
                       name = 'Potential bloom\nfrequencies:\n ', 
                       option = 'D',
                       limits = c(0,100)) +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
          color = 'white') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
          color = 'white') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
          color = 'white') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
          color = 'white')  +
  geom_sf(data = st_transform(st_sym_difference(florida_bound,lake_bound),wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),color = 'black', fill = 'transparent') +
  geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'white', shape = 1) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1,'cm'),
                         pad_x = unit(-0.15,'in'),
                         pad_y = unit(0.3,'in')) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

july_plot

ggarrange(WY2019_plot,july_plot,
          common.legend = T,legend = 'bottom'
          )

ggsave(
  filename = 'heat_maps.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 4,
  units = 'in'
)

make_heat_maps <- function(annual_or_monthly,year_or_month){
  if(annual_or_monthly == 'annual') {
    stack_classified <- stack_weekly_classified[[which(WY(dates_weekly_start) == year_or_month)]] 
    stack_count <- stack_classified
    for(i in 1:nlayers(stack_classified)) {
      stack_count[[i]][stack_count[[i]] < 2] <- 0
      stack_count[[i]][stack_count[[i]] >=  2] <- 1
    }
    calc(x = stack_count,fun = sum, na.rm = T)/nlayers(stack_count)*100 -> heat_map
    heat_map[heat_map == 0] <- NA
  } else {
    stack_classified <- stack_weekly_classified[[which(month(dates_weekly_start) == year_or_month)]] 
    stack_count <- stack_classified
    for(i in 1:nlayers(stack_classified)) {
      stack_count[[i]][stack_count[[i]] < 2] <- 0
      stack_count[[i]][stack_count[[i]] >=  2] <- 1
    }
    calc(x = stack_count,fun = sum, na.rm = T)/nlayers(stack_count)*100 -> heat_map
    heat_map[heat_map == 0] <- NA
  }
  
  ggplot() +
    layer_spatial(data = crop(projectRaster(heat_map,    
                                            crs = wgs84,
                                            method = 'bilinear'),st_transform(lake_bound,wgs84))) +
    scale_fill_viridis_c(na.value = 'transparent', 
                         name = 'Total potential bloom\nfrequency for each pixel\n ', 
                         option = 'D',
                         limits = c(0,100)) +
    geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent',
            color = 'white') +
    geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent',
            color = 'white') +
    geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent',
            color = 'white') +
    geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent',
            color = 'white')  +
    geom_sf(data = st_transform(st_sym_difference(florida_bound,lake_bound),wgs84), alpha = .8) +
    geom_sf(data = st_transform(lake_bound,wgs84),color = 'black', fill = 'transparent') +
    geom_sf(data = st_transform(lake_station_coords,wgs84), color = 'white', shape = 1) +
    coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
    theme_minimal() +
    annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
    annotation_north_arrow(style = north_arrow_minimal,
                           height = unit(1,'cm'),
                           pad_x = unit(-0.15,'in'),
                           pad_y = unit(0.3,'in')) +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)
    ) 
}

make_heat_maps('monthly', 7)

ggsave(
  filename = 'july_heat_map.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5.5,
  height = 4,
  units = 'in'
)

heat_map_july_proj <- projectRaster(heat_map_july, crs = wgs84,method = 'bilinear')

lake_station_coords %>%
  mutate(freq = raster::extract(heat_map_july_proj,st_transform(lake_station_coords,wgs84),method = 'bilinear')
)
raster::extract(heat_map_july_proj,st_transform(lake_station_coords,wgs84),method = 'bilinear')

lake_station_coords %>%
mutate(bloom_freq = raster::extract(heat_map_july_proj,
                                    st_transform(lake_station_coords,wgs84),
                                    method = 'bilinear')
)

ggplot() +
  # layer_spatial(data = crop(projectRaster(heat_map,    
  #                                         crs = wgs84,
  #                                         method = 'ngb'),st_transform(lake_bound,wgs84))) +
  # scale_fill_viridis_c(na.value = 'transparent', 
  #                      name = 'Total potential bloom\nfrequency for each pixel\n ', 
  #                      option = 'D',
  #                      limits = c(0,100)) +
  geom_sf(data = nw_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = ne_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = sw_quad_84, alpha = 0, fill = 'transparent') +
  geom_sf(data = se_quad_84, alpha = 0, fill = 'transparent')  +
  geom_sf(data = st_transform(st_sym_difference(florida_bound,lake_bound),wgs84), alpha = .8) +
  geom_sf(data = st_transform(lake_bound,wgs84),color = 'black', fill = 'transparent') +
  geom_sf(data = st_transform(lake_station_coords,wgs84), aes(color = station_id)) +
  coord_sf(xlim = c(-81.12,-80.6),ylim = c(26.68,27.2)) +
  theme_minimal() +
  annotation_scale(location = 'bl', height = unit(0.2,'cm'), width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_minimal,
                         height = unit(1,'cm'),
                         pad_x = unit(-0.15,'in'),
                         pad_y = unit(0.3,'in')) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))


# NOT used in analysis -- very early stages ----
# plot_dominant_class <- function(tile_vals){
#   tile_vals %>%
#     rename(month = `month(date)`,year = `year(date)`) %>%
#     mutate(
#       month = factor(month,labels = month.abb[unique(month)]),
#       year = factor(year, levels = seq(2016,2021,1),ordered = T),
#       dominant_class = factor(dominant_class,levels = c(5,4,3,2,1),ordered = T)) %>%
#     ggplot() +
#     geom_tile(aes(x = month,y = year,fill = dominant_class)) +
#     xlab('Month') + ylab('Year') +
#     scale_fill_manual(values = colors_5,name = 'Quintile') +
#     theme_minimal() + 
#     theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) %>% return()
# }
# 
# plot_percent_bloom <- function(tile_vals) {
#   tile_vals %>%
#     rename(month = `month(date)`,year = `year(date)`) %>%
#     mutate(
#       month = factor(month,labels = month.abb[unique(month)]),
#       year = factor(year, levels = seq(2016,2021,1),ordered = T)) %>%
#     ggplot() +
#     geom_tile(aes(x = month,y = year,fill = percent_bloom)) +
#     scale_fill_viridis_c(name = 'Percent') +
#     xlab('Month') + ylab('Year') +
#     theme_minimal() + 
#     theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) %>% return()
# }
# 
# plot_dominant_class(lake_tile_vals)
# plot_dominant_class(roi_tile_vals)
# 
# bind_rows(mutate(nw_tile_vals,quad = 'Northwest Quadrant'),
#           mutate(ne_tile_vals,quad = 'Northeast Quadrant'),
#           mutate(sw_tile_vals,quad = 'Southwest Quadrant'),
#           mutate(se_tile_vals,quad = 'Southeast Quadrant')) %>%
#   mutate(quad = factor(quad,
#                        levels = c('Northwest Quadrant','Northeast Quadrant',
#                                   'Southwest Quadrant','Southeast Quadrant')
#                        , ordered = T)
#   ) %>%
#   plot_dominant_class() + facet_wrap(~quad) +
#   theme(
#     axis.text.x = element_text(margin = margin(0,0,0,0)),
#     strip.background = element_rect(fill = 'transparent'),
#     panel.spacing = unit(1,'line')
#   )
# 
# 
# plot_percent_bloom(lake_tile_vals)
# plot_percent_bloom(roi_tile_vals)
# 
# bind_rows(mutate(nw_tile_vals,quad = 'Northwest Quadrant'),
#           mutate(ne_tile_vals,quad = 'Northeast Quadrant'),
#           mutate(sw_tile_vals,quad = 'Southwest Quadrant'),
#           mutate(se_tile_vals,quad = 'Southeast Quadrant')) %>%
#   plot_percent_bloom() + facet_wrap(~quad)

# Save work space - last saved 08/10/22 to check all calculations ---
save.image("~/Documents/R projects/olci_ci_copy/frequency_analysis_ci_cyano_COPY_workspace.RData")

