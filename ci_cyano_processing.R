# Processing CI cyano images from Blake

# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(raster)
library(sf)
library(beepr)
library(progress)

# Set Working Directory
setwd('/Users/natalie/Documents/select_data')

# Read in data ----

## Rasters
# Daily
folder_daily <- 'ci_cyano_images/Daily/'
filenames_daily <- list.files(folder_daily)
filepaths_daily <- str_c(folder_daily,filenames_daily)
stack_daily_unscaled <- raster::stack(filepaths_daily)
beep(2)

# Weekly
folder_weekly <- 'ci_cyano_images/7D/'
filenames_weekly <- list.files(folder_weekly)
filepaths_weekly <- str_c(folder_weekly,filenames_weekly)
stack_weekly_unscaled <- raster::stack(filepaths_weekly)
beep(2)

## Shape files
lake_bound <-  read_sf('Shape_WBD/WBDHU12.shp') %>% filter(HUC12 == '030902010200') %>%
  st_transform(crs(stack_daily_unscaled))

# Crop images ----

lake_buffer_extent <- extent(st_buffer(lake_bound, 5000)) 

## Daily
stack_daily_unscaled_cropped <- crop(stack_daily_unscaled, lake_buffer_extent)
beep(2)

## Weekly
stack_weekly_unscaled_cropped <- crop(stack_weekly_unscaled, lake_buffer_extent)
beep(2)

# Scale images ----

scale_images <- function(raster_image) {
  
  # Create template raster
  x <- raster_image
  scaled_image <- raster(x = extent(x), res = res(x), crs = crs(x), vals = rep(NA, ncell(x)))
  
  # Scale image to CI
  scaled_image <- 10^((3/250)*x - 4.2)
  
  # Remove NA values
  na_sub <- which(x[] == 0 | x[] == 254 | x[] == 255)
  scaled_image[na_sub] <- NA
  
  return(scaled_image)
}

## Daily
# Create an empty raster stack and progress bar
stack_daily_scaled <- raster::stack()
pb = txtProgressBar(min = 0, 
                    max = nlayers(stack_daily_unscaled_cropped), 
                    initial = 0)

# Populate -- ONLY RUN ONCE
for(i in 1:nlayers(stack_daily_unscaled_cropped)) {
  image_i <- stack_daily_unscaled_cropped[[i]]
  
  image_i_scaled <- scale_images(image_i)
  
  stack_daily_scaled <- raster::stack(
    stack_daily_scaled,
    image_i_scaled
  )
  
  setTxtProgressBar(pb,i)
}
beep(2)

## Weekly
# Create an empty raster stack and progress bar
stack_weekly_scaled <- raster::stack()
pb = txtProgressBar(min = 0, 
                    max = nlayers(stack_weekly_unscaled_cropped), 
                    initial = 0)

# Populate -- ONLY RUN ONCE
for(i in 1:nlayers(stack_weekly_unscaled_cropped)) {
  image_i <- stack_weekly_unscaled_cropped[[i]]
  
  image_i_scaled <- scale_images(image_i)
  
  stack_weekly_scaled <- raster::stack(
    stack_weekly_scaled,
    image_i_scaled
  )
  
  setTxtProgressBar(pb,i)
}
beep(2)

# Process images - Tomlinson and Seegers ---- 

ci_to_chla <- function(raster_image,coef,const) {
  # x <- stack_daily_scaled[[1]]
  # coef <- 4050
  # const <- 20
  
  # Create template raster
  x <- raster_image
  chla_image <- raster(x = extent(x), res = res(x), crs = crs(x), vals = rep(NA, ncell(x)))
  
  # Grab subset of values that are NOT NA values
  val_sub <- which(!is.na(x[]))
  
  # Process image
  chla_image[val_sub] <- coef*x[val_sub] + const
  
  return(chla_image)
}

coef_tom <- 4050
const_tom <- 20

coef_seeg <- 6620
const_seeg <- -3.07

## Daily
# Create an empty raster stack and progress bar
stack_daily_chla_tom <- raster::stack()
stack_daily_chla_seeg <- raster::stack()
pb = txtProgressBar(min = 0, 
                    max = nlayers(stack_daily_scaled), 
                    initial = 0)

# Populate -- ONLY RUN ONCE
for(i in 1:nlayers(stack_daily_scaled)) {
  image_i <- stack_daily_scaled[[i]]
  
  image_i_tom <- ci_to_chla(image_i,coef = coef_tom,const = const_tom)
  image_i_seeg <- ci_to_chla(image_i,coef = coef_seeg,const = const_seeg)
  
  stack_daily_chla_tom <- raster::stack(
    stack_daily_chla_tom,
    image_i_tom
  )
  
  stack_daily_chla_seeg <- raster::stack(
    stack_daily_chla_seeg,
    image_i_seeg
  )
  
  setTxtProgressBar(pb,i)
}
beep(2)

## Weekly
# Create an empty raster stack and progress bar
stack_weekly_chla_tom <- raster::stack()
stack_weekly_chla_seeg <- raster::stack()
pb = txtProgressBar(min = 0, 
                    max = nlayers(stack_weekly_scaled), 
                    initial = 0)

# Populate -- ONLY RUN ONCE
for(i in 1:nlayers(stack_weekly_scaled)) {
  image_i <- stack_weekly_scaled[[i]]
  
  image_i_tom <- ci_to_chla(image_i,coef = coef_tom,const = const_tom)
  image_i_seeg <- ci_to_chla(image_i,coef = coef_seeg,const = const_seeg)
  
  stack_weekly_chla_tom <- raster::stack(
    stack_weekly_chla_tom,
    image_i_tom
  )
  
  stack_weekly_chla_seeg <- raster::stack(
    stack_weekly_chla_seeg,
    image_i_seeg
  )
  
  setTxtProgressBar(pb,i)
}
beep()

## Rename layers ----
names(stack_daily_chla_tom) <- names(stack_daily_unscaled_cropped)
names(stack_daily_chla_seeg) <- names(stack_daily_unscaled_cropped)

names(stack_weekly_chla_tom) <- names(stack_weekly_unscaled_cropped)
names(stack_weekly_chla_seeg) <- names(stack_weekly_unscaled_cropped)

# Save environment ----
# rm(image_i,image_i_scaled,image_i_seeg,image_i_tom,chla_image,scaled_image,x,i,na_sub,val_sub,pb)
# save.image("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")