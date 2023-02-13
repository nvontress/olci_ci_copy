# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(sf)
library(janitor)
library(raster)
library(lubridate)
library(beepr)
library(ggpubr)
library(modelr)
library(dbhydroR)


# Set Working Directory
setwd('/Users/natalie/Documents/select_data')

# Read in data ----

# Load ci_cyano_processing environment
load("~/Documents/R projects/olci_ci/ci_cyano_processing_workspace.RData")

# Set coordinate system
epsg <- crs(stack_daily_unscaled)
wgs84 <- 4326

# Seegers' 2021 national data
ntl_matchups <- read_csv('CyANChlBS_Matchups_2021June23.csv')
ntl_validation <- read_csv('CIcyanotestSet.csv')

ntl_validation %>% 
  mutate(ci_chla = coef_tom*`CIcyano Test Set` + const_tom) %>%
  ggplot() + 
  geom_point(aes(x = log10(`In situ Chl (ug/ml)`),y = log10(ci_chla))) +
  xlim(0,3) + ylim(0,3)
  
# Okay so it turns out that Bridget didn't actually send me all the info I needed, but I can definitely work backwards from this.  
# in her excel sheet, each X1 should be unique. I also have her resulting MAE and bias. SO I can try different combos of unique X1 values until I get the resulting MAE and bias. Is this a waste of time?? I think I can write up a script pretty quickly to accomplish this, honestly.
  
rename(ntl_validation, MERIS_ci_cyano = `CI Test Set`) %>% 
  left_join(dplyr::select(ntl_matchups, MERIS_ci_cyano,chl)) %>% distinct() %>%
  group_by(X1) %>% count() %>% filter(n>1)

## DBHYDRO
lake_stations_all <- c('S308C',"KISSR0.0","L001","L004","L005","L006","L007","L008","LZ2","LZ25A","LZ30","LZ40","PALMOUT","POLESOUT")

lake_stations <- c(#'S308C',
  "KISSR0.0","L001","L004","L005","L006","L007","L008",
  # "LZ2",
  "LZ25A","LZ30","LZ40"
  # "PALMOUT","POLESOUT"
)

get_wq(station_id = lake_stations,
       date_min = '2016-01-01',
       date_max = '2020-12-31',
       test_name = 'CHLOROPHYLL-A(LC)') %>% as_tibble() %>%
  mutate(date = date(date)) -> lake_chla

## Lake station coords
read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
  filter(station %in% lake_stations,
         activity_t == 'Chemistry',
         activity_s == 'Surface Water Grab') %>%
  dplyr::select(station, geometry) %>%
  rename(station_id = station) %>%
  st_transform(crs = epsg) -> lake_station_coords

# Pull dates ----

## Raster dates
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

## Sampling dates
lake_chla %>% 
  pull(date) -> sampling_dates

# Tibble setup ----

## Shared dates
## Create a tibble of dates shared by image capture and sampling dates
semi_join(
  tibble(date = dates_daily),
  tibble(date = sampling_dates)
) -> shared_dates

## Select corresponding raster images
tibble(date = dates_daily) %>%
  mutate(index = row_number()) %>% 
  semi_join(shared_dates) %>% 
  pull(index) -> raster_index

stack_daily_scaled[[raster_index]] -> stack_daily_scaled_filtered

## Select corresponding chla samples
lake_chla %>% 
  filter(date %in% shared_dates$date) -> lake_chla_filtered

## Set up tibble
df <- lake_chla_filtered %>%
  gather(colnames(lake_chla_filtered)[2:ncol(lake_chla_filtered)],key = station_id,value = value) %>%
  mutate(station_id = str_extract(station_id,'.+(?=_CHLOROPHYLL)')) %>% 
  rename(field_val = value)

ci_extracted_vals <- tibble(
  station_id = lake_stations
)

for(i in 1:nrow(shared_dates)){
  bind_cols(ci_extracted_vals,
            as_tibble_col(x = as.numeric(NA), 
                          column_name = as.character(pull(shared_dates,date)[[i]]))) -> ci_extracted_vals
} 

for(i in 2:ncol(ci_extracted_vals)){
  
  raster::extract(stack_daily_scaled_filtered[[i-1]],lake_station_coords) %>% as_tibble() -> ci_extracted_vals[1:nrow(ci_extracted_vals),i]
  
}
beep(2)

## Restructure 
ci_extracted_vals %>% 
  gather(colnames(ci_extracted_vals)[2:ncol(ci_extracted_vals)],
         key = 'date', value = 'ci_val') %>%
  mutate(date = ymd(date)) -> ci_vals 

joined_data <- right_join(df,ci_vals) %>% na.omit()

## Calculate chla
joined_data %>% 
  mutate(
    chla_tom = coef_tom*ci_val + const_tom,
    chla_who = ci_val*50000
  ) -> joined_data

## Calculate MAE
calc_mae <- function(actual,predicted){
  abs_error <- abs(predicted-actual)
  return(mean(abs_error))
}

calc_mae_log <- function(observed,predicted){
  abs_error <- abs(log10(predicted)-log10(observed))
  return(10^(mean(abs_error)))
}

calc_bias_log <- function(observed,predicted){
  num <- (log10(predicted)-log10(observed))
  return(10^mean(num))
}

mae_tom_daily <- calc_mae_log(observed = joined_data$field_val,predicted = joined_data$chla_tom)

bias_tom_daily <- calc_bias_log(observed = joined_data$field_val,predicted = joined_data$chla_tom)

# Examine for hypereutrophic range
head(joined_data)
joined_data %>% filter(field_val > 30 & field_val <= 90) -> joined_data_he

mae_tom_daily_he <- calc_mae_log(observed = joined_data_he$field_val,predicted = joined_data_he$chla_tom)

bias_tom_daily_he <- calc_bias_log(observed = joined_data_he$field_val,predicted = joined_data_he$chla_tom)

ggplot() +
  geom_point(data = ntl_validation,
             aes(x = `In situ Chl (ug/ml)`, y = (coef_tom*`CIcyano Test Set` + const_tom)),
             color = 'pink',
             alpha = 0.5) +
  geom_point(data = joined_data,
             aes(field_val,chla_tom),
             color = 'red',
             alpha = 0.5) +
  geom_abline(slope = 1,intercept = 0) +
  xlim(0,850) +
  ylim(0,850) +
  theme_minimal()

## Remove largest CI value because it has wayyyy too much leverage
joined_data %>%
  arrange(ci_val) %>%
  slice(1:(nrow(joined_data)-1)) -> joined_data_cut

# Repeat with weekly composites ----
sampling_dates_week_start <- rep(NA, length(sampling_dates))
for(i in 1:length(sampling_dates_week_start)) {
  # i <- 1
  sampling_dates_week_start[i] <- as_date(sampling_dates[i]) - wday(sampling_dates[i]) + 1
}

semi_join(
  tibble(week_date_start = dates_weekly_start,
         week_date_end = dates_weekly_end),
  tibble(week_date_start = as_date(sampling_dates_week_start))
) -> shared_dates_weekly

## Select corresponding raster images
tibble(week_date_start = dates_weekly_start) %>%
  mutate(index = row_number()) %>% 
  semi_join(shared_dates_weekly) %>% 
  pull(index) -> composite_index

stack_weekly_scaled[[composite_index]] -> stack_weekly_scaled_filtered

## Select corresponding chla samples
lake_chla %>% 
  filter((date - wday(date) + 1) %in% shared_dates_weekly$week_date_start) -> lake_chla_filtered_weekly

## Set up tibble
df_weekly <- lake_chla_filtered_weekly %>%
  gather(colnames(lake_chla_filtered_weekly)[2:ncol(lake_chla_filtered_weekly)],key = station_id,value = value) %>%
  mutate(station_id = str_extract(station_id,'.+(?=_CHLOROPHYLL)')) %>% 
  dplyr::select(date,station_id,value) %>%
  rename(field_val = value) %>%
  mutate(week_date_start = (date - wday(date) + 1))

ci_extracted_vals_weekly <- tibble(
  station_id = lake_stations
)

for(i in 1:nrow(shared_dates_weekly)){
  bind_cols(ci_extracted_vals_weekly,
            as_tibble_col(x = as.numeric(NA), 
                          column_name = as.character(
                            pull(shared_dates_weekly,week_date_start)[[i]]))) -> ci_extracted_vals_weekly
} 

for(i in 2:ncol(ci_extracted_vals_weekly)){
  
  raster::extract(stack_weekly_scaled_filtered[[i-1]],lake_station_coords) %>% as_tibble() -> ci_extracted_vals_weekly[1:nrow(ci_extracted_vals_weekly),i]
  
}
beep(2)

## Restructure 
ci_extracted_vals_weekly %>% 
  gather(colnames(ci_extracted_vals_weekly)[2:ncol(ci_extracted_vals_weekly)],
         key = 'week_date_start', value = 'ci_val') %>%
  mutate(week_date_start = ymd(week_date_start)) -> ci_vals_weekly

joined_data_weekly <- right_join(df_weekly,ci_vals_weekly) %>% na.omit()

joined_data_weekly %>%
  mutate(
    chla_tom = coef_tom*ci_val + const_tom,
    chla_who = ci_val*50000
  ) -> joined_data_weekly

mae_tom_weekly <- calc_mae(actual = joined_data_weekly$field_val,predicted = joined_data_weekly$chla_tom)
mae_who_weekly <-  calc_mae(actual = joined_data_weekly$field_val,predicted = joined_data_weekly$chla_who)

ggplot(joined_data_weekly) +
  geom_abline(slope = 1,intercept = 0,color = 'red') + 
  geom_point(aes(x = field_val,y = chla_tom), alpha = .5, size = 2) +
  xlim(c(0,300)) +
  ylim(c(0,300)) +
  xlab('Observed chlorophyll concentration') +
  ylab('Simulated chlorophyll concentration') +
  theme_minimal() +
  theme(text = element_text(size = 15))

# Plot for manuscript ----
ggplot() + 
  geom_abline(aes(slope = 1,intercept = 0,linetype = 'twodash'),color = 'lightblue4') + 
  geom_abline(aes(slope = 1,intercept = 0,linetype = 'solid'),color = 'black') + 
  geom_point(aes(x = log10(`In situ Chl (ug/ml)`), y = log10(coef_tom*`CIcyano Test Set` + const_tom),
                 color = 'lightblue3',), 
             data = ntl_validation, alpha = .5, size = 1.8) +
  geom_point(aes(x = log10(field_val),y = log10(chla_tom), color = 'skyblue4'),
             shape = 1,
             data = joined_data, size = 1.8) +
  # geom_hline(yintercept = log10(40), color = 'blue') +
  # geom_vline(xintercept = log10(40), color = 'blue') +
  # geom_hline(aes(yintercept = log10(30)), color = 'lightblue4',linetype = 'twodash') +
  geom_vline(aes(xintercept = log10(30)), color = 'lightblue4',linetype = 'twodash') +
  # geom_hline(aes(yintercept = log10(90)), color = 'lightblue4',linetype = 'twodash') +
  geom_vline(aes(xintercept = log10(90)), color = 'lightblue4',linetype = 'twodash') +
  # geom_hline(yintercept = log10(25), color = 'green') +
  # geom_vline(xintercept = log10(25), color = 'green') +
  # geom_hline(yintercept = log10(55), color = 'green') +
  # geom_vline(xintercept = log10(55), color = 'green') +
  xlab(expression('Log10 transformed observed chlorophyll concentration')) +
  ylab('Log10 transformed simulated chlorophyll concentration') +
  # xlim(c(0,3)) +
  # ylim(c(0,3)) +
  scale_x_continuous(limits = c(10^-1,3)) +
  scale_y_continuous(limits = c(10^-1,3)) +
  scale_color_identity(name = 'Data source:',
                       breaks = c('lightblue3','skyblue4'),
                       labels = c('Seegers et al. (2021)',
                                  'Lake Okeechobee'),
                       #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
                       guide = 'legend') +
  scale_fill_identity(breaks = c('transparent')) +
  scale_linetype_identity(name = 'Line types:',
                          breaks = c('solid','twodash'),
                          labels = c('1:1 Slope reference',
                                     'Hypereutrophic range'),
                          guide = 'legend') +
  guides(linetype = guide_legend(override.aes = list(linetype = c('solid', 'twodash'),color = c("black", "lightblue4"))),
         color = guide_legend(override.aes = list(shape = c(19,1), alpha = c(0.5,1)))) +
  theme_classic() +
  theme(text = element_text(size = 12))

ggsave(
  filename = 'one_to_one_plot_word.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 5,
  units = 'in'
)

ggplot() + 
  geom_abline(aes(slope = 1,intercept = 0,linetype = 'twodash'),color = 'lightblue4') + 
  geom_abline(aes(slope = 1,intercept = 0,linetype = 'solid'),color = 'black') + 
  geom_point(aes(x = `In situ Chl (ug/ml)`, y = coef_tom*`CIcyano Test Set` + const_tom,
                 color = 'lightblue3'), 
             data = ntl_validation, alpha = .4, size = 1.8) +
  geom_point(aes(x = field_val,y = chla_tom, color = 'skyblue4'),
             data = joined_data, shape = 1, size = 1.8) +
  # geom_hline(aes(yintercept = 30), color = 'lightblue4',linetype = 'twodash') +
  geom_vline(aes(xintercept = 30), color = 'lightblue4',linetype = 'twodash') +
  # geom_hline(aes(yintercept = 90), color = 'lightblue4',linetype = 'twodash') +
  geom_vline(aes(xintercept = 90), color = 'lightblue4',linetype = 'twodash') +
  xlab('*In situ* chlorophyll-*a* (µg L^(-1))') +
  ylab('Satellite-derived  
       chlorophyll-*a* (µg L^(-1))') +
  # xlim(c(0,3)) +
  # ylim(c(0,3)) +
  scale_x_continuous(limits = c(10^-1,850),trans = 'log10',
                     breaks = c(1,10,100,1000),
                     labels = c(1,10,100,1000)) +
  scale_y_continuous(limits = c(10^-1,850),trans = 'log10',
                     breaks = c(1,10,100,1000),
                     labels = c(1,10,100,1000)) +
  scale_color_identity(name = 'Data source:',
                       breaks = c('lightblue3','skyblue4'),
                       labels = c('Seegers et al. (2021)',
                                  'Lake Okeechobee'),
                       #(expression(paste("Taylor et al. 2012 \n", italic("(366 adults)"))))
                       guide = 'legend') +
  scale_linetype_identity(name = 'Line types:',
                          breaks = c('solid','twodash'),
                          labels = c('1:1 Slope reference',
                                     'Hypereutrophic range'),
                          guide = 'legend') +
  guides(linetype = guide_legend(override.aes = list(linetype = c('solid', 'twodash'),
                                                     color = c("black", "lightblue4"))),
         color = guide_legend(override.aes = list(shape = c(19,1), alpha = c(0.5,1)))) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        panel.border = element_rect(color = 'black', fill = 'transparent', size = 1),
        axis.line = element_blank())

ggsave(
  filename = 'one_to_one_plot_word_3.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5.25,
  height = 3,
  units = 'in'
)

# labs(
#   x = "Groups",
#   y = "No. of *bacteria X* isolates with corresponding types",
#   fill = "Var1"
# ) +
#   theme(axis.title.y = ggtext::element_markdown())


# quickly visualize
joined_data_weekly %>%
  ggplot(aes(x = ci_val,y = field_val)) +
  geom_point() +
  geom_smooth() +
  theme_minimal() +
  labs(title = 'Weekly Composites')

joined_data_cut %>%
  ggplot(aes(x = ci_val,y = field_val)) +
  geom_point() +
  geom_smooth() +
  theme_minimal() +
  labs(title = 'Individual Images')

# Save workspace
# save.image("~/Documents/R projects/olci_ci/chla_mod_ci_cyano_workspace.RData")

# Visualize the data ----
ggdensity(joined_data_cut$field_val)
ggqqplot(joined_data_cut$field_val)

# log10(y) -- REALLY BAD QQ PLOT
ggdensity(log10(joined_data$field_val))
ggqqplot(log10(joined_data_cut$field_val))

# Box-Cox -- this verifies the log transform shouldn't be used bc lambda does not contain zero (from stats 516 notes)
MASS::boxcox(field_val ~ ci_val, data = joined_data_cut, plotit = T)

# The following three transformations result in long-tail distributions and, therefore, least squares regressions should not be used 
## sqrt(y) - BEST DENSITY PLOT/QQPLOT COMBO
ggdensity(sqrt(joined_data_cut$field_val))
ggqqplot(sqrt(joined_data_cut$field_val))

## y^(1/3)
ggdensity((joined_data_cut$field_val)^(1/3))
ggqqplot((joined_data_cut$field_val)^(1/3))

## y^(1/4) -- THIS SEEMS BEST FOR DENSITY PLOT
ggdensity((joined_data_cut$field_val)^(1/4))
ggqqplot((joined_data_cut$field_val)^(1/4))

## 1/sqrt(y) - MAKES WORSE
ggdensity(1/(sqrt(joined_data_cut$field_val)))

## 1/y - MAKES WORSE
ggdensity(1/joined_data_cut$field_val)

# General visualizations
joined_data %>%
  ggplot(aes(x = ci_val,y = field_val)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

joined_data_cut %>%
  ggplot(aes(x = ci_val,y = field_val)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

joined_data_cut %>%
  ggplot(aes(x = ci_val,y = sqrt(field_val))) +
  geom_point() +
  geom_smooth(method = 'gam') +
  theme_minimal()

joined_data_cut %>%
  ggplot(aes(x = log10(ci_val),y = field_val)) +
  geom_point() +
  geom_smooth() +
  theme_minimal() # This gives the most consistent variance

joined_data_cut %>%
  ggplot(aes(x = log10(ci_val),y = log10(field_val))) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

# Model, no transformation ----
mod <- lm(field_val ~ ci_val, data = joined_data_cut) 
summary(mod)

par(mfrow = c(2,2))
plot(mod) ## Errors are NOT normally distributed

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod),x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

# Transformations ----
# Square root of field value -- Really horrible diagnostic plots
mod_2 <- lm(sqrt(field_val) ~ ci_val, data = joined_data_cut)
plot(mod_2)

# Log10 transformed ci -- Really horrible diagnostic plots
mod_log_ci <- lm(field_val ~ log10(ci_val), data = joined_data_cut)
summary(mod_log_ci)
plot(mod_log_ci)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_log_ci),x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

# Log transformed ci and field -- Better diagnostic plots but still not great
mod_log <- lm(log10(field_val) ~ log10(ci_val), data = joined_data_cut)
summary(mod_log)

plot(mod_log)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_log),
                 x = log10(joined_data_cut$field_val))) + 
  xlim(0,2.5) + 
  ylim(0,2.5) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

## Back Transformed
ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = 10^(predict(mod_log)),
                 x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

# Transformed field
mod_glm <- glm(sqrt(field_val) ~ ci_val, data = joined_data_cut)
summary(mod_glm)
plot(mod_glm)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red')

# Non-linear models ----

## Quadratic
mod_quad <- lm(field_val ~ poly(ci_val,2), data = joined_data_cut)
summary(mod_quad)

par(mfrow = c(2,2))
plot(mod_quad)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_quad),x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

mod_quad_2 <- lm(sqrt(field_val) ~ poly(ci_val,2), data = joined_data_cut)
summary(mod_quad_2)

par(mfrow = c(2,2))
plot(mod_quad_2)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = (predict(mod_quad_2))^2,x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

joined_data_cut %>%
  ggplot(aes(x =ci_val,y = sqrt(field_val))) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ poly(x,2)) +
  theme_minimal()

## Cubic -- DOESN'T WORK
mod_cube <- glm(field_val ~ ci_val + I(ci_val^2) + I(ci_val^3), data = joined_data_cut)
summary(mod_cube)

par(mfrow = c(2,2))
plot(mod_cube)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_cube),x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

mod_cube_2 <- glm(sqrt(field_val) ~ ci_val + I(ci_val^2) + I(ci_val^3), data = joined_data_cut)
summary(mod_cube)

par(mfrow = c(2,2))
plot(mod_cube_2)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y =( predict(mod_cube_2))^2,x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

## Quartic -- ALSO DOESN'T WORK
mod_quart <- lm(field_val ~ ci_val + I(ci_val^2) + I(ci_val^3) + I(ci_val^4), data = joined_data)
summary(mod_quart)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_quart),x = joined_data$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

# LOESS
mod_loess <- loess(field_val ~ ci_val, data = joined_data_cut, family = 'symmetric')
summary(mod_loess)

plot(mod_loess)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_loess),x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

joined_data_cut %>%
  ggplot(aes(x = ci_val,y = field_val)) + 
  geom_point() +
  geom_smooth(method = 'loess',span = .75) +
  theme_minimal()

## Oddly, the transformed loess doesn't seem to work as well as the un-transformed
mod_loess_2 <- loess(sqrt(field_val) ~ ci_val, data = joined_data_cut, family = 'symmetric')
summary(mod_loess_2)

plot(mod_loess_2)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = (predict(mod_loess_2))^2,x = joined_data_cut$field_val)) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

# Multi-linear regression
mod_multi <- lm(sqrt(field_val) ~ ci_val + month(date), data = joined_data_cut)
par(mfrow = c(2,2))
plot(mod_multi)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = (predict(mod_multi))^2,x = joined_data_cut$field_val)) +  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nmulti-linear regression mod (μg/L)') +
  theme_minimal()
  
# Piecewise ----
# Fitting a model from [0,0.0005]

joined_data_cut %>%
  filter(ci_val < 0.0005) -> joined_data_filtered

mod_piece_1 <- lm(field_val ~ ci_val, joined_data_filtered)
summary(mod_piece_1)

par(mfrow = c(2,2))
plot(mod_piece_1)

mae(mod_piece_1,joined_data_filtered)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_piece_1),x = joined_data_filtered$field_val)) +  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nfiltered regression model (μg/L)') +
  theme_minimal()

ggplot() + 
  geom_abline(intercept = mod_piece_1$coefficients[[1]], 
              slope = mod_piece_1$coefficients[[2]],
              color = 'red') +
  geom_point(data = joined_data_filtered, aes(x = ci_val, y = field_val)) +
  theme_minimal()

mod_piece_2 <- lm(sqrt(field_val) ~ ci_val, joined_data_filtered)
summary(mod_piece_2)

par(mfrow = c(2,2))
plot(mod_piece_2)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = (predict(mod_piece_2))^2,x = joined_data_filtered$field_val)) +  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nfiltered regression model with transformation (μg/L)') +
  theme_minimal()

# ???
ggplot() + 
  geom_abline(intercept = mod_piece_2$coefficients[[1]], 
              slope = mod_piece_2$coefficients[[2]],
              color = 'red') +
  geom_point(data = joined_data_filtered, aes(x = ci_val, y = sqrt(field_val))) +
  theme_minimal()

# Messing around -- VERY UNORGANIZED CODE ----
joined_data %>%
  # ggplot(aes(x = log10(ci_val),y = field_val)) +
  ggplot(aes(x = ci_val,y = field_val)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

joined_data %>%
  # ggplot(aes(x = log10(ci_val),y = field_val)) +
  ggplot(aes(x = ci_val,y = log10(field_val))) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

joined_data %>%
  arrange(ci_val) %>%
  slice(1:(nrow(joined_data)-1)) %>%
  # ggplot(aes(x = log10(ci_val),y = field_val)) +
  ggplot(aes(x = ci_val,y = log10(field_val))) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

loess(field_val ~ ci_val, data = (joined_data %>%
        arrange(ci_val) %>%
        slice(1:(nrow(joined_data)-1))), family = 'symmetric', 
      span = .02) -> mod_loess

summary(mod_quad)

mod_loess <- loess(field_val ~ ci_val, data = joined_data, family = 'symmetric')
summary(mod_loess)

ggplot() +
  geom_abline(intercept = 0,slope = 1,color = 'red') +
  geom_point(aes(y = predict(mod_loess),x = (joined_data %>%
                                               arrange(ci_val) %>%
                                               slice(1:(nrow(joined_data)-1)) %>%
                                               pull(field_val)
                                             ))) + 
  xlim(0,150) + 
  ylim(0,150) +
  labs(title = '1:1 plot for concentrations from \nlinear regression mod (μg/L)') +
  theme_minimal()

library(splines)
mod_spline <- lm(field_val ~ bs(ci_vals, knots = (5*10^-4)), data = joined_data)

attr(bs(joined_data$ci_val),'knots')

joined_data %>%
  arrange(ci_val) %>%
  slice(1:(nrow(joined_data)-1)) %>%
  # ggplot(aes(x = log10(ci_val),y = field_val)) +
  ggplot(aes(x = log10(ci_val),y = (field_val))) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

loess(field_val ~ log10(ci_val), data = joined_data_cut, family = 'symmetric', 
      span = .02) -> mod_loess


joined_data_cut %>%
  ggplot(aes(x =ci_val,y = sqrt(field_val))) +
  geom_point() +
  geom_smooth(method = 'loess', family = 'symmetric') +
  theme_minimal()

write_csv(joined_data,file = 'extracted_values.csv')
