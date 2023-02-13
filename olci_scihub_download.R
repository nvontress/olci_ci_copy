# Author: Natalie Von Tress
# Last edited: 06/30/2020
# Description: This code serves to retrieve OLCI LFR/EFR images from SciHub and unzip files to a separate directory

# A few notes: WFR images cannot be accessed via SciHub, but rather have to be downloaded via Creodias by Eumetsat
# That being said, this code CAN be used to create a query of image names, that can be slightly edited to quickly call WFR images during processing

# Prepare workspace ----

# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(getSpatialData)
library(sf)
library(beepr)

# Set working directory
setwd('/Users/natalie/Documents/select_data')

# Load THIS workspace
load("~/Documents/R projects/olci_ci_copy/sentinel_records_workspace.RData")

# Retrieve data ----

# Log in to SciHub
login_CopHub(username = 'natnatvt') # Enter password to prompt

# Set path to download
set_archive('OLCI_LFR_zipped', create = T)

# Set AOI to Lake Okeechobee polygon
## Load Lake Okeechobee shape file
# sbtools::item_get(sb_id = '5a1632b9e4b09fc93dd17214') # This item was downloaded from USGS ScienceBase to serve as Lake Okeechobee boundary
lake_bound  <- read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200')
set_aoi(lake_bound$geometry)
# view_aoi()

# Retrieve products
records <- getSentinel_records(time_range = c('2016-05-01','2021-05-01'),
                             products = 'sentinel-3',
                             aoi = lake_bound$geometry) 
beep(2)

records %>% as_tibble() %>% colnames()

records %>% as_tibble() %>% filter(sensor_id == 'OLCI') %>%
  pull(record_id) %>% unique() 

# sensor_id == 'OLCI'

records %>% filter(sensor_id == 'OLCI',
                   product_type == 'OL_2_LFR___',
                   str_detect(record_id,'S3A') == T,
                   str_detect(record_id,'2520') == T) %>%
  arrange(record_id) %>%
  pull(start_time) %>% as.Date() %>% unique() %>% length()

records %>% filter(sensor_id == 'OLCI',
                   product_type == 'OL_2_LFR___',
                   str_detect(record_id,'S3A') == T,
                   str_detect(record_id,'2520') == T)  %>%
  pull(start_time) %>% as.Date() %>% unique() -> s3a_dates_verified

records %>% filter(sensor_id == 'OLCI',
                   product_type == 'OL_2_LFR___',
                   str_detect(record_id,'S3B') == T,
                   str_detect(record_id,'2520') == T) %>% 
  pull(start_time) %>% as.Date %>% unique() -> s3b_dates_verified

records %>% filter(sensor_id == 'OLCI',
                   product_type == 'OL_2_LFR___',
                   str_detect(record_id,'S3A') == T | 
                     str_detect(record_id,'S3B') == T,
                   str_detect(record_id,'2520') == T) %>% 
  pull(start_time) %>% as.Date() %>% unique() -> s3a_s3b_dates_verified

# Save workspace -- Last saved 1/20 to pull S3A AND S3B flyover dates
# save.image("~/Documents/R projects/sentinel_records_workspace_verify.RData")

  # filter(str_detect(filename,pattern = 'LFR') == TRUE) %>% 
  # arrange(filename)
# write_csv(x = records, path = '/Users/natalie/Documents/Data/OLCI_LFR_Records.csv') 

# Download images ----

# WARNING: This code may use up all your computer storage!
# ONLY RUN IF YOU ARE CERTAIN YOU HAVE THE MEMORY SPACE AND NEED THESE IMAGES

# for(i in 1:nrow(records)) {
#   skip_to_next_index <- F
#   tryCatch(expr = getSentinel_data(records = records[i,]), 
#            error = function(e) {skip_to_next_index <<- T})
#   if(skip_to_next_index) {next}
# } # Downloads all available images, skipping those which are archived

# Unzip files ----
# 
# for(i in 1:nrow(records)) {
#   zipfile_name <- str_split(records$filename[i], pattern = "\\.")[[1]][1] %>% str_c('zip', sep = '.')
#   zipfile <- str_c('OLCI_LFR/get_data/Sentinel-3/', zipfile_name, sep = '')
#   
#   skip_to_next_index <- F
#   tryCatch(expr = unzip(zipfile = zipfile, exdir = "OLCI_LFR_Unzipped"), 
#            error = function(e) {skip_to_next_index <<- T})
#   if(skip_to_next_index) {next}
# }
