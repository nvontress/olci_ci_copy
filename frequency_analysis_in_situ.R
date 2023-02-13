# The purpose of this script is to compare in situ data with the frequency analysis 

# Set up work space ----
# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(sf)
library(dbhydroR)
library(janitor)
library(raster)
library(exactextractr)
library(colormap)
library(lubridate)
library(relayer)
library(ggtext)
library(zoo)
library(AnalystHelper)
library(effectsize)
library(ggpubr)
library(ggpattern)
library(beepr)

# Set Working Directory
# setwd('/Volumes/Hard Drive/Data/') # for hard drive
# -- or --
setwd('/Users/natalie/Documents/select_data/') # for data on laptop

# Load THIS WORKSPACE -- Last updated 10/27 after updating legend labels
load("~/Documents/R projects/olci_ci_copy/frequency_analysis_in_situ_COPY_workspace.RData")

# Load workspace from frequency analysis - not necessary if you ran line 26
# load("~/Documents/R projects/olci_ci_copy/frequency_analysis_ci_cyano_COPY_workspace.RData")

# Load in situ data
# When I remove all the S308C, POLESOUT, LZ2, and PALMOUT your points start to align with Bridget’s national response.
# lake_stations <- c(#'S308C',
#                    "KISSR0.0","L001","L004","L005","L006","L007","L008",
#                    #"LZ2",
#                    "LZ25A","LZ30","LZ40"
#                    # ,"PALMOUT","POLESOUT"
#                    )

# Need internet connection for this step
get_wq(station_id = lake_stations,
       date_min = '2016-05-01',
       date_max = '2021-04-30',
       test_name = 'CHLOROPHYLL-A(LC)') %>% as_tibble() -> lake_chla

# read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
#   filter(station %in% lake_stations,
#          activity_t == 'Chemistry',
#          activity_s == 'Surface Water Grab') %>%
#   dplyr::select(station, geometry) %>%
#   rename(station_id = station) %>%
#   st_transform(crs = epsg) -> lake_station_coords

get_wq(station_id = 'S308C',
       date_min = '2016-05-01',
       date_max = '2021-04-30',
       test_name = 'CHLOROPHYLL-A(LC)') %>% as_tibble() -> s308c_chla

# s308c_coords <- read_sf('dbhydro_stations/DBHYDRO_SITE_STATION.shp') %>% clean_names() %>%
#   filter(station == 'S308C',
#          activity_t == 'Chemistry',
#          activity_s == 'Surface Water Grab') %>%
#   dplyr::select(station, geometry) %>%
#   rename(station_id = station) %>%
#   st_transform(crs = epsg)

# Read in hydrologic data
s308_usgs <- read_csv('NWISdata/20160501-20210430_s308.csv')
# beep(2)

# Generate value summary tibbles for in situ data ----
# generate_vals_tibble_in_situ <- function(chla = lake_chla,
#                                          station_coords = lake_station_coords,
#                                          bound) {
#   # chla <- lake_chla
#   # station_coords <- lake_station_coords
#   # bound <- nw_quad
#   
#   station_coords$station_id[st_intersects(station_coords,bound,sparse = F) %>% 
#                               which(T) %>% as_tibble() %>% pull(row)] -> stations
#   
#   chla %>%
#     gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = value) %>%
#     mutate(station_id = str_extract(station_id,'.+(?=_CHLOROPHYLL)')) %>% 
#     filter(station_id %in% stations,
#            !is.na(value)) %>%
#     mutate(station_id = factor(x = station_id, labels = stations), 
#            date = date(date),
#            ci_est = (value - 20)/4050,
#            class = cut(x = ci_est,breaks = breaks_weekly,labels = F)) %>%
#     group_by(date,class) %>% count() %>% spread(class,n) %>% return()
# }
# 
# generate_vals_tibble_in_situ(bound = lake_bound) -> lake_vals_chla
# generate_vals_tibble_in_situ(bound = nw_quad) -> nw_vals_chla
# generate_vals_tibble_in_situ(bound = ne_quad) -> ne_vals_chla
# generate_vals_tibble_in_situ(bound = sw_quad) -> sw_vals_chla
# generate_vals_tibble_in_situ(bound = se_quad) -> se_vals_chla

# Bloom/no bloom frequencies ----
generate_bloom_tib <- function(chla = lake_chla,
                               station_coords = lake_station_coords,
                               bound) {
  
  station_coords$station_id[st_intersects(station_coords,bound,sparse = F) %>% 
                              which(T) %>% as_tibble() %>% pull(row)] -> stations
    
  chla %>% 
    gather(colnames(lake_chla)[2:ncol(lake_chla)],key = station_id,value = value) %>%
    mutate(station_id = str_extract(station_id,'.+(?=_CHLOROPHYLL)')) %>% 
    filter(station_id %in% stations,
           !is.na(value)) %>%
    dplyr::select(date,station_id,value) %>%
    mutate(station_id = factor(x = station_id, labels = stations), 
           date = date(date),
           bloom = cut(x = value,breaks = c(0,40,max(value)),labels = F)) %>%
    group_by(date,bloom) %>% count() %>% spread(bloom,n) %>%
    rename(no_bloom = `1`,
           bloom = `2`) %>% 
  return()
  
}

generate_bloom_tib(bound = lake_bound) -> lake_vals_bloom
generate_bloom_tib(bound = nw_quad) -> nw_vals_bloom
generate_bloom_tib(bound = ne_quad) -> ne_vals_bloom
generate_bloom_tib(bound = sw_quad) -> sw_vals_bloom
generate_bloom_tib(bound = se_quad) -> se_vals_bloom

## Annual ----
annual_freq_bloom <- function(tib) {
  tib %>%
    group_by(WY(date)) %>%
    summarize(
      total_no_bloom = sum(no_bloom,na.rm = T),
      total_bloom = sum(bloom,na.rm = T),
      n = sum(total_bloom,total_no_bloom)
    ) %>%
    mutate(
      freq_no_bloom = total_no_bloom/n*100,
      freq_bloom = total_bloom/n*100
    ) %>% return()
}

annual_freq_lake_bloom <- annual_freq_bloom(lake_vals_bloom)
annual_freq_nw_quad_bloom <- annual_freq_bloom(nw_vals_bloom)
annual_freq_ne_quad_bloom <- annual_freq_bloom(ne_vals_bloom)
annual_freq_sw_quad_bloom <- annual_freq_bloom(sw_vals_bloom)
annual_freq_se_quad_bloom <- annual_freq_bloom(se_vals_bloom)

### Plot ----
# scales::show_col(
#   colormap(
#     colormap = colormaps$salinity,
#     nshades = 12
#   )
# )



# plot_annual_freq_bloom <- function(annual_bloom_freq_tib) {
#   # annual_bloom_freq_tib <- annual_freq_lake_bloom
#   annual_bloom_freq_tib %>%
#     dplyr::select(-total_no_bloom,-total_bloom,-n) %>%
#     rename(Year = `WY(date)`,
#            Bloom = freq_bloom,
#            `No bloom` = freq_no_bloom) %>%
#     mutate(`Water Year` = factor(Year,levels = unique(Year),ordered = T)) %>%
#     gather(`No bloom`,Bloom,key = 'Class',value = 'Frequency') %>%
#     ggplot(aes(x = `Water Year`,y = Frequency,fill = Class)) +
#     geom_bar(position = 'stack',stat = 'identity') +
#     scale_fill_manual(values = colors_2,name = 'Classification') +
#     scale_y_continuous(breaks = seq(0,100,10),
#                        minor_breaks = NULL) +
#     theme_minimal() +    
#     theme(text = element_text(size = 14),
#           axis.text.x = element_text(margin = margin(-10,0,5,0)),
#           axis.text.y = element_text(margin = margin(0,0,0,3)),
#           panel.grid.major.x = element_blank()) %>% return()
# }
# 
# plot_annual_freq_bloom(annual_freq_lake_bloom) 
# 
# bind_rows(mutate(annual_freq_nw_quad_bloom,quad = 'Northwest Quadrant'),
#           mutate(annual_freq_ne_quad_bloom,quad = 'Northeast Quadrant'),
#           mutate(annual_freq_sw_quad_bloom,quad = 'Southwest Quadrant'),
#           mutate(annual_freq_se_quad_bloom,quad = 'Southeast Quadrant')) %>%
#   plot_annual_freq_bloom() + facet_wrap(~quad)  +
#   theme(axis.text.x = element_text(margin = margin(0,0,0,0)))  

## Monthly ----
monthly_freq_bloom <- function(tib) {
  tib %>%
    group_by(month(date)) %>%
    summarize(
      total_no_bloom = sum(no_bloom,na.rm = T),
      total_bloom = sum(bloom,na.rm = T),
      n = sum(total_bloom,total_no_bloom)
    ) %>%
    mutate(
      freq_no_bloom = total_no_bloom/n*100,
      freq_bloom = total_bloom/n*100
    ) %>% return()
}

monthly_freq_lake_bloom <- monthly_freq_bloom(lake_vals_bloom)
monthly_freq_nw_quad_bloom <- monthly_freq_bloom(nw_vals_bloom)
monthly_freq_ne_quad_bloom <- monthly_freq_bloom(ne_vals_bloom)
monthly_freq_sw_quad_bloom <- monthly_freq_bloom(sw_vals_bloom)
monthly_freq_se_quad_bloom <- monthly_freq_bloom(se_vals_bloom)

### Plot ----
# plot_monthly_freq_bloom <- function(monthly_bloom_freq_tib) {
#   # monthly_bloom_freq_tib <- monthly_freq_lake_bloom
#   monthly_bloom_freq_tib %>%
#     dplyr::select(-total_no_bloom,-total_bloom,-n) %>%
#     rename(Month = `month(date)`,
#            Bloom = freq_bloom,
#            `No bloom` = freq_no_bloom) %>%
#     mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
#     gather(`No bloom`,Bloom,key = 'Class',value = 'Frequency') %>%
#     ggplot(aes(x = Month,y = Frequency,fill = Class)) +
#     geom_bar(position = 'stack',stat = 'identity') +
#     scale_fill_manual(values = colors_2,name = 'Classification') +
#     scale_y_continuous(breaks = seq(0,100,10),
#                        minor_breaks = NULL) +
#     theme_minimal()  + 
#     theme(text = element_text(size = 14),
#           axis.text.x = element_text(angle = 45,
#                                      margin = margin(-5,0,0,0)),
#           axis.text.y = element_text(margin = margin(0,0,0,5)),
#           panel.grid.major.x = element_blank()) %>% return()
# }
# 
# plot_monthly_freq_bloom(monthly_freq_lake_bloom) 
# 
# bind_rows(mutate(monthly_freq_nw_quad_bloom,quad = 'Northwest Quadrant'),
#           mutate(monthly_freq_ne_quad_bloom,quad = 'Northeast Quadrant'),
#           mutate(monthly_freq_sw_quad_bloom,quad = 'Southwest Quadrant'),
#           mutate(monthly_freq_se_quad_bloom,quad = 'Southeast Quadrant')) %>%
#   plot_monthly_freq_bloom() + facet_wrap(~quad) +
#   theme(axis.text.x = element_text(margin = margin(0,0,0,0)))  

# Hydrographs ----
month_labels_WY <- c('M','J','J','A','S','O','N','D','J','F','M','A')
month_labels <- c('J','F','M','A','M','J','J','A','S','O','N','D')
year_labels <- rep(month_labels_WY,5)

# dates_daily[which(WY(dates_daily) != 2016 & WY(dates_daily) != 2022)] -> dates_daily

s308_usgs %>%
  mutate(WY = WY(dateTime),
         month = month(dateTime)) %>%
  filter(WY != 2016 & WY != 2022) -> s308_usgs_filtered

# s308_usgs_filtered %>% group_by(WY,month) %>% summarize(mean_flow = mean(flow_cfs))

# Figuring out which mean to use for plotting - THIS FIG IS NOT FOR THE MANUSCRIPT
s308_usgs_filtered %>%
  group_by(date(dateTime)) %>% 
  summarize(daily_flow = mean(flow_cfs)) %>%
  mutate(avg_7d = zoo::rollmean(daily_flow, k = 7, fill = NA),
         avg_14d = zoo::rollmean(daily_flow, k = 14, fill = NA)) %>%
  ggplot() +
  geom_line(aes(x = `date(dateTime)`, y = daily_flow,color = 'black')) +
  geom_line(aes(x = `date(dateTime)`, y = avg_7d,color = 'red')) +
  geom_line(aes(x = `date(dateTime)`, y = avg_14d,color = 'blue')) +
  scale_color_identity(name = 'Flow rolling average',
                       breaks = c('black','red','blue'), 
                       labels = c('1-day','7-day','14-day'), 
                       guide = 'legend') +
  xlab('Date') +
  ylab('Average flow (cfs)') +
  theme_minimal()
  
# start_of_month_marks <- tibble(
#   day = '01',
#   month = as.character(month(c(5:12,rep(seq(1,12),4),1:4),label = T)),
#   year = as.character(c(rep(2016,8),rep(2017,12),
#            rep(2018,12),rep(2019,12),rep(2020,12),rep(2021,4))),
#   date = dmy(str_c(day,month,year))
#     ) %>%
#   dplyr::select(year,date) %>%
#   mutate(
#     year = as.numeric(year),
#       month = c(5:12,rep(seq(1,12),4),1:4)
#     )

s308_usgs_filtered %>%
  group_by(date(dateTime)) %>% 
  summarize(daily_flow = mean(flow_cfs)) %>%
  mutate(avg_7d = zoo::rollmean(daily_flow, k = 7, fill = NA)) %>%
  rename(date = `date(dateTime)`) %>%
  mutate(month = month(date),
         year = WY(date)) -> s308_usgs_filtered

# Hydrograph for manuscript ----
WY_labs <- c('WY 2017','WY 2018','WY 2019','WY 2020','WY 2021')
names(WY_labs) <- as.character(2017:2021)

roi_tile_vals %>%
    # mutate(total_bloom = total_2 + total_3,
    #        bloom_freq = total_bloom/(count_pixels(roi)*n)*100) %>%
  dplyr::select(month,year,total_pot_bloom_freq) %>%  # year references WY
  full_join(s308_usgs_filtered) %>%
  ggplot() +
  geom_hline(yintercept = 0,color = 'red',linetype = 'twodash') +
  geom_line(aes(x = date,y = avg_7d)) +
  geom_rug(aes(x = date,color = total_pot_bloom_freq),
           sides = 'b') +
  facet_grid(col = vars(year),scales = 'free_x',
             labeller = labeller(year = WY_labs)) +
  xlab('Month') + ylab('Flow (cfs)') +
  scale_x_date(breaks = seq(from=as.Date("2016-05-01"),
                            to=as.Date("2021-05-01"),
                            by="month"),
               minor_breaks = NULL,
               labels = c(year_labels,'M')) +
  scale_color_colormap(aesthetics = 'color', 
                       name = 'Total potential\nbloom frequency',
                       discrete = F,
                       colormap = colormaps$salinity #,
                       # limits = c(0,50)
                       ) +
  theme_minimal() +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.text.x = element_text(hjust = -0.3,vjust = 1),
        text = element_text(size = 12)) 

# ggsave(
#   filename = 'hydrograph_word_redo.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 9,
#   height = 5,
#   units = 'in'
# )

# Get quartiles of roi_tile_vals -- referenced in section 3.2.2 of manuscript
roi_tile_vals %>% pull(total_pot_bloom_freq) %>% quantile(probs = seq(0,1,.1)) %>% 
  as_tibble() %>% mutate(percentile = seq(0,100,10)) %>%
  ggplot(aes(x = percentile,y = value)) +
  geom_col() +
  geom_text(aes(label = round(value,digits = 1),vjust = -.7)) +
  scale_x_continuous(breaks = seq(0,100,10)) +
  theme_minimal()


# Hydrograph for pptx
roi_tile_vals %>%
  rename(
    month = `month`,
    year = `year`,
  ) %>%
  mutate(total_bloom = total_2 + total_3,
         bloom_freq = total_bloom/(count_pixels(roi)*n)*100) %>%
  dplyr::select(month,year,bloom_freq) %>% # year references WY
  full_join(s308_usgs_filtered) %>%
  ggplot() +
  geom_hline(yintercept = 0,color = 'red',linetype = 'twodash') +
  geom_line(aes(x = date,y = avg_7d)) +
  geom_rug(aes(x = date,color = bloom_freq),
           sides = 'b') +
  facet_grid(col = vars(year),scales = 'free_x',
             labeller = labeller(year = WY_labs)) +
  xlab('Month') + ylab('Flow (cfs)') +
  scale_x_date(breaks = seq(from=as.Date("2016-05-01"),
                            to=as.Date("2021-05-01"),
                            by="month"),
               minor_breaks = NULL,
               labels = c(year_labels,'M')) +
  scale_color_colormap(aesthetics = 'color', 
                       name = 'Total potential\nbloom frequency',
                       discrete = F,
                       colormap = colormaps$salinity,
                       limits = c(0,50)) +
  theme_minimal() +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.text.x = element_text(hjust = -0.3,vjust = 1),
        text = element_text(size = 16))

# ggsave(
#   filename = 'hydrograph_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 11.25,
#   height = 5.5,
#   units = 'in'
# )

# Plots for manuscript ----
colors_2 <- c(colormap(colormaps$salinity,12)[11], # lime
              colormap(colormaps$salinity,12)[5]) # green

colors_3 <- c(colormap(colormaps$salinity,12)[11], # lime
              colormap(colormaps$salinity,12)[9], # green
              colormap(colormaps$salinity,12)[5]) # teal

plot_annual_side_by_side_bloom <- function(annual_freq_rs, annual_freq_bl){
  # RS = remote sensing
  # bl = bloom
  # 
  # annual_freq_rs <- annual_freq_lake_weekly
  # annual_freq_bl <- annual_freq_lake_bloom
  
  annual_freq_rs %>%
    dplyr::select(-total_1,-total_2,-total_3,-n) %>%
    rename(
      Year = `WY(date)`,
      # `0` = freq_0,
      `1` = freq_1,
      `2` = freq_2,
      `3` = freq_3,
      # `NA` = freq_NA
    ) %>%
    # mutate(Year = factor(Year,labels = Year.abb[unique(Year)],ordered = T)) %>%
    gather(`1`,`2`,`3`,key = 'class_rs',value = 'freq') %>%
    mutate(source = 'Remote sensing') -> rs_tib
  
  annual_freq_bl %>%
    dplyr::select(-total_no_bloom,-total_bloom,-n) %>%
    rename(
      Year = `WY(date)`,
      `No bloom` = freq_no_bloom,
      Bloom = freq_bloom
    ) %>%
    # mutate(Year = factor(Year,labels = Year.abb[unique(Year)],ordered = T)) %>%
    gather(`No bloom`,Bloom,key = 'class_is',value = 'freq') %>%
    mutate(source = 'In situ') -> bl_tib
  
  w <- .35
  b <- .4
  a <- .9
  
  bind_rows(rs_tib,bl_tib) %>%
    mutate(class_rs = factor(class_rs,levels = c(3,2,1),labels = c('Bloom','Possible bloom','No bloom'),ordered = T),
           class_is = factor(class_is,levels = c('Bloom','No bloom'), ordered = T)) %>%
    spread(source,freq) %>%
    ggplot() +
    # Plot
    geom_col(aes(x = Year,y = `Remote sensing`,fill = class_rs),
             alpha = a,
             position = 'stack',
             width = w) +    
    geom_col_pattern(aes(x = Year + b,y = `In situ`,fill2 = class_is),
                     alpha = a,
                     position = 'stack',
                     width = w,
                     pattern = 'stripe', 
                     pattern_color = 'white',
                     pattern_fill = 'white',
                     pattern_alpha = 0.3,
                     pattern_density = 0.05)  %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
    # Scale
    scale_y_continuous(breaks = seq(0,100,20)) +
    scale_fill_manual(aesthetics = 'fill',
                      values = colors_3,
                      name = 'Remote sensing\nclassifications',
                      drop = F) +
    scale_fill_manual(aesthetics = 'fill2',
                      values = colors_2,
                      name = expression(italic('In situ')~'\nclassifications')) +
    ylab('Frequency') + xlab('Water Year') +
    # Format
    theme_minimal() + 
    theme(text = element_text(size = 12),
          axis.text.x = element_text(vjust = 1,
                                     hjust = 0),
          axis.text.y = element_text(margin = margin(0,0,0,0)),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.ticks = element_blank()) +
    guides(fill = guide_legend(override.aes = list(lwd = 0))) +
    # Clip
    coord_cartesian(expand = F) %>% return()
}

plot_annual_side_by_side_bloom(annual_freq_lake_weekly,annual_freq_lake_bloom)

plot_monthly_side_by_side_bloom <- function(monthly_freq_rs, monthly_freq_bl){
  # RS = remote sensing
  # bl = bloom
  
  # monthly_freq_rs <- monthly_freq_lake_weekly
  # monthly_freq_bl <- monthly_freq_lake_bloom
  
  monthly_freq_rs %>%
    dplyr::select(-total_1,-total_2,-total_3,-n) %>%
    rename(
      Month = `month(date)`,
      `1` = freq_1,
      `2` = freq_2,
      `3` = freq_3,
    ) %>%
    # mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
    gather(`1`,`2`,`3`,key = 'class_rs',value = 'freq') %>%
    mutate(source = 'Remote sensing') -> rs_tib
  
  monthly_freq_bl %>%
    dplyr::select(-total_no_bloom,-total_bloom,-n) %>%
    rename(
      Month = `month(date)`,
      `No bloom` = freq_no_bloom,
      Bloom = freq_bloom
    ) %>%
    # mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
    gather(`No bloom`,Bloom,key = 'class_is',value = 'freq') %>%
    mutate(source = 'In situ') -> bl_tib
  
  w <- .35
  b <- .4
  a <- .9
  
  bind_rows(rs_tib,bl_tib) %>%
    mutate(class_rs = factor(class_rs,levels = c(3,2,1),labels = c('Bloom','Possible bloom','No bloom'),ordered = T),
           class_is = factor(class_is,levels = c('Bloom','No bloom'), ordered = T))%>%
    spread(source,freq) %>%
    ggplot() +
    # Plot
    geom_col(aes(x = Month,y = `Remote sensing`,fill = class_rs),
             alpha = a,
             position = 'stack',
             width = w) +    
    geom_col_pattern(aes(x = Month + b,y = `In situ`,fill2 = class_is),
             alpha = a,
             position = 'stack',
             width = w,
             pattern = 'stripe', 
             pattern_color = 'white',
             pattern_fill = 'white',
             pattern_alpha = 0.3,
             pattern_density = 0.05)  %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
    # Scale
    scale_x_discrete(name = 'Month',
                     # breaks = seq(1,12,1),
                     labels = month_labels,
                     limits = seq(1,12,1)) +  
    scale_y_continuous(breaks = seq(0,100,20)) +
    scale_fill_manual(aesthetics = 'fill',
                      values = colors_3,
                      name = 'Remote sensing\nclassifications',
                      drop = F) +
    scale_fill_manual(aesthetics = 'fill2',
                      values = colors_2,
                      name = expression(italic('In situ')~'\nclassifications')) +
    ylab('Frequency') +
    # Format
    theme_minimal() + 
    theme(text = element_text(size = 12),
          axis.text.x = element_text(vjust = 1,
                                     hjust = 0),
          axis.text.y = element_text(margin = margin(0,0,0,0)),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.ticks = element_blank()) +
    guides(fill = guide_legend(override.aes = list(lwd = 0))) +
    # Clip
    coord_cartesian(expand = F) %>% return()
}

annual_lake_plot <- plot_annual_side_by_side_bloom(annual_freq_lake_weekly,annual_freq_lake_bloom)

ggsave(
  plot = annual_lake_plot,
  filename = 'annual_freq_lake_plot_redo_2.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5,
  height = 4,
  units = 'in',
  dpi = 300
)

monthly_lake_plot <- plot_monthly_side_by_side_bloom(monthly_freq_lake_weekly,monthly_freq_lake_bloom)

ggsave(
  plot = monthly_lake_plot,
  filename = 'monthly_freq_lake_plot_redo_2.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 5,
  height = 4,
  units = 'in',
  dpi = 300
)

# for AGU overview slides ----
annual_plot <- plot_annual_side_by_side_bloom(annual_freq_lake_weekly,annual_freq_lake_bloom)

monthly_plot <- plot_monthly_side_by_side_bloom(monthly_freq_lake_weekly,monthly_freq_lake_bloom)

ggarrange(annual_plot,
          monthly_plot,
          ncol = 1,
          nrow = 2,
          legend = 'right',
          legend.grob = get_legend(
            p = annual_plot,
            position = 'right'))

ggsave(filename = 'frequencies_plot_agu_overview.png',
       path = '/Users/natalie/Documents/thesis_writing/r_plots',
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300)

# Annual quadrants ----
bind_rows(mutate(annual_freq_nw_quad_weekly,quad = 'Northwest Quadrant'),
          mutate(annual_freq_ne_quad_weekly,quad = 'Northeast Quadrant'),
          mutate(annual_freq_sw_quad_weekly,quad = 'Southwest Quadrant'),
          mutate(annual_freq_se_quad_weekly,quad = 'Southeast Quadrant')) %>%
  mutate(quad = factor(quad,
                       levels = c('Northwest Quadrant','Northeast Quadrant',
                                  'Southwest Quadrant','Southeast Quadrant')
                       , ordered = T)
  ) ->  annual_freq_quads_weekly

bind_rows(mutate(annual_freq_nw_quad_bloom,quad = 'Northwest Quadrant'),
          mutate(annual_freq_ne_quad_bloom,quad = 'Northeast Quadrant'),
          mutate(annual_freq_sw_quad_bloom,quad = 'Southwest Quadrant'),
          mutate(annual_freq_se_quad_bloom,quad = 'Southeast Quadrant')) %>%
  mutate(quad = factor(quad,
                       levels = c('Northwest Quadrant','Northeast Quadrant',
                                  'Southwest Quadrant','Southeast Quadrant')
                       , ordered = T)
  ) ->  annual_freq_quads_bloom

annual_quads_plot <- plot_annual_side_by_side_bloom(annual_freq_quads_weekly,annual_freq_quads_bloom) + 
  facet_wrap(~quad) +
  theme(
    axis.text.x = element_text(hjust = .2),
    strip.background = element_rect(fill = 'transparent'),
    panel.spacing = unit(1,'line')
  )

ggsave(
  plot = annual_quads_plot,
  filename = 'annual_freq_quads_plot_redo_2.jpg',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 4,
  units = 'in',
  dpi = 300
)

# Monthly quadrants ----

bind_rows(mutate(monthly_freq_nw_quad_weekly,quad = 'Northwest Quadrant'),
          mutate(monthly_freq_ne_quad_weekly,quad = 'Northeast Quadrant'),
          mutate(monthly_freq_sw_quad_weekly,quad = 'Southwest Quadrant'),
          mutate(monthly_freq_se_quad_weekly,quad = 'Southeast Quadrant')) %>%
  mutate(quad = factor(quad,
                       levels = c('Northwest Quadrant','Northeast Quadrant',
                                  'Southwest Quadrant','Southeast Quadrant')
                         , ordered = T)
         ) -> monthly_freq_quads_weekly

bind_rows(mutate(monthly_freq_nw_quad_bloom,quad = 'Northwest Quadrant'),
          mutate(monthly_freq_ne_quad_bloom,quad = 'Northeast Quadrant'),
          mutate(monthly_freq_sw_quad_bloom,quad = 'Southwest Quadrant'),
          mutate(monthly_freq_se_quad_bloom,quad = 'Southeast Quadrant'))  %>%
  mutate(quad = factor(quad,
                       levels = c('Northwest Quadrant','Northeast Quadrant',
                                  'Southwest Quadrant','Southeast Quadrant')
                       , ordered = T)
  ) ->  monthly_freq_quads_bloom

monthly_quads_plot <- plot_monthly_side_by_side_bloom(monthly_freq_quads_weekly,monthly_freq_quads_bloom) + 
  facet_wrap(~quad) +
  theme(
    axis.text.x = element_text(vjust = 1,hjust = .15),
    strip.background = element_rect(fill = 'transparent'),
    panel.spacing = unit(1,'line')
  )

ggsave(
  plot = monthly_quads_plot,
  filename = 'monthly_freq_quads_plot_redo_2.jpg',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 4,
  units = 'in',
  dpi = 300
)

# Combined plots ----

annual_full_plot <- ggarrange(
  annual_lake_plot + 
    ggtitle('Entire Lake') + 
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = 0.1)),
  annual_quads_plot + theme_minimal() + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12,vjust = -1),
          axis.text.x = element_text(hjust = .2),
          panel.spacing = unit(1,'line'),
          strip.text = element_text(size = 10)),
  ncol = 1, 
  nrow = 2,
  heights = c(1,2),
  legend = 'right',
  common.legend = T
) %>%
  annotate_figure(left = text_grob('Frequency (%)', size = 12, rot = 90))

ggsave(
  plot = annual_full_plot,
  filename = 'annual_freq_plot_all.jpg',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 7,
  units = 'in',
  dpi = 300
)

monthly_full_plot <- ggarrange(
  monthly_lake_plot + 
    ggtitle('Entire Lake') + 
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = 0.1)),
  monthly_quads_plot + theme_minimal() + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 12,vjust = -1),
          axis.text.x = element_text(hjust = .2),
          panel.spacing = unit(1,'line'),
          strip.text = element_text(size = 10)),
  ncol = 1, 
  nrow = 2,
  heights = c(1,2),
  align = 'v',
  legend = 'right',
  common.legend = T
) %>%
  annotate_figure(left = text_grob('Frequency (%)', size = 12, rot = 90))

ggsave(
  plot = monthly_full_plot,
  filename = 'monthly_freq_plot_all.jpg',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 6.5,
  height = 7,
  units = 'in',
  dpi = 300
)

# Plots for presentations - need to remake ----

# Functions
# plot_annual_side_by_side_bloom_pptx <- function(annual_freq_rs, annual_freq_bl){
#   # RS = remote sensing
#   # bl = bloom
#   
#   # annual_freq_rs <- annual_freq_lake_weekly
#   # annual_freq_bl <- annual_freq_lake_bloom
#   
#   annual_freq_rs %>%
#     dplyr::select(-total_0,-total_1,-total_2,-total_3,-total_4,-total_5,-total_NA,-n) %>%
#     rename(
#       Year = `year(date)`,
#       `0` = freq_0,
#       `1` = freq_1,
#       `2` = freq_2,
#       `3` = freq_3,
#       `4` = freq_4,
#       `5` = freq_5,
#       `NA` = freq_NA
#     ) %>%
#     # mutate(Year = factor(Year,labels = Year.abb[unique(Year)],ordered = T)) %>%
#     gather(`NA`,`0`,`1`,`2`,`3`,`4`,`5`,key = 'class_rs',value = 'freq') %>%
#     mutate(source = 'Remote sensing') -> rs_tib
#   
#   annual_freq_bl %>%
#     dplyr::select(-total_no_bloom,-total_bloom,-n) %>%
#     rename(
#       Year = `year(date)`,
#       `No bloom` = freq_no_bloom,
#       Bloom = freq_bloom
#     ) %>%
#     # mutate(Year = factor(Year,labels = Year.abb[unique(Year)],ordered = T)) %>%
#     gather(`No bloom`,Bloom,key = 'class_is',value = 'freq') %>%
#     mutate(source = 'In situ') -> bl_tib
#   
#   w <- .35
#   b <- .4
#   a <- .9
#   
#   colors_6 <- c(colormap(colormaps$salinity,12)[11:10],
#                 colormap(colormaps$salinity,12)[6:4],
#                 colormap(colormaps$salinity,12)[1])
#   
#   bind_rows(rs_tib,bl_tib) %>%
#     mutate(class_rs = factor(class_rs,levels = c(5,4,3,2,1,0,'NA'),labels = c('Q5','Q4','Q3','Q2','Q1','Below MDL','Missing values'),ordered = T),
#            class_is = factor(class_is,levels = c('Bloom','No bloom'), ordered = T)) %>%
#     spread(source,freq) %>%
#     ggplot() +
#     # Plot
#     geom_col(aes(x = Year,y = `Remote sensing`,fill = class_rs),
#              alpha = a,
#              position = 'stack',
#              width = w) +    
#     geom_col(aes(x =(Year + b),y = `In situ`,fill2 = class_is),
#              alpha = a,
#              position = 'stack',
#              width = w)  %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
#     scale_y_continuous(breaks = seq(0,100,20)) +
#     scale_fill_manual(aesthetics = 'fill',
#                       values = c(colors_6,'black'),
#                       name = 'Remote sensing\nclassifications',
#                       drop = F) +
#     scale_fill_manual(aesthetics = 'fill2',
#                       values = colors_2,
#                       name = expression(italic('In situ')~'\nclassifications')) +
#     ylab('Frequency') +
#     # Format
#     theme_minimal() + 
#     theme(text = element_text(size = 16),
#           legend.text = element_markdown(),
#           axis.text.x = element_text(hjust = 0.15),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.major.x = element_blank(),
#           axis.ticks = element_blank()) +
#     # Clip
#     coord_cartesian(expand = F) %>% return()
# }
# 
# plot_monthly_side_by_side_bloom_pptx <- function(monthly_freq_rs, monthly_freq_bl){
#   # RS = remote sensing
#   # bl = bloom
#   
#   # monthly_freq_rs <- monthly_freq_lake_weekly
#   # monthly_freq_bl <- monthly_freq_lake_bloom
#   
#   monthly_freq_rs %>%
#     dplyr::select(-total_0,-total_1,-total_2,-total_3,-total_4,-total_5,-total_NA,-n) %>%
#     rename(
#       Month = `month(date)`,
#       `0` = freq_0,
#       `1` = freq_1,
#       `2` = freq_2,
#       `3` = freq_3,
#       `4` = freq_4,
#       `5` = freq_5,
#       `NA` = freq_NA
#     ) %>%
#     # mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
#     gather(`NA`,`0`,`1`,`2`,`3`,`4`,`5`,key = 'class_rs',value = 'freq') %>%
#     mutate(source = 'Remote sensing') -> rs_tib
#   
#   monthly_freq_bl %>%
#     dplyr::select(-total_no_bloom,-total_bloom,-n) %>%
#     rename(
#       Month = `month(date)`,
#       `No bloom` = freq_no_bloom,
#       Bloom = freq_bloom
#     ) %>%
#     # mutate(Month = factor(Month,labels = month.abb[unique(Month)],ordered = T)) %>%
#     gather(`No bloom`,Bloom,key = 'class_is',value = 'freq') %>%
#     mutate(source = 'In situ') -> bl_tib
#   
#   w <- .35
#   b <- .4
#   a <- .9
#   
#   colors_6 <- c(colormap(colormaps$salinity,12)[11:10],
#                 colormap(colormaps$salinity,12)[6:4],
#                 colormap(colormaps$salinity,12)[1])
#   
#   bind_rows(rs_tib,bl_tib) %>%
#     mutate(class_rs = factor(class_rs,levels = c(5,4,3,2,1,0,'NA'),labels = c('Q5','Q4','Q3','Q2','Q1','Below MDL','Missing values'),ordered = T),
#            class_is = factor(class_is,levels = c('Bloom','No bloom'), ordered = T)) %>%
#     spread(source,freq) %>%
#     ggplot() +
#     # Plot
#     geom_col(aes(x = Month,y = `Remote sensing`,fill = class_rs),
#              alpha = a,
#              position = 'stack',
#              width = w) +    
#     geom_col(aes(x = Month + b,y = `In situ`,fill2 = class_is),
#              alpha = a,
#              position = 'stack',
#              width = w)  %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
#     # Scale
#     scale_x_discrete(name = 'Month',
#                      # breaks = seq(1,12,1),
#                      labels = month_labels,
#                      limits = seq(1,12,1)) +  
#     scale_y_continuous(breaks = seq(0,100,20)) +
#     scale_fill_manual(aesthetics = 'fill',
#                       values = c(colors_6,'black'),
#                       name = 'Remote sensing\nclassifications',
#                       drop = F) +
#     scale_fill_manual(aesthetics = 'fill2',
#                       values = colors_2,
#                       name = expression(italic('In situ')~'\nclassifications')) +
#     ylab('Frequency') +
#     # Format
#     theme_minimal() + 
#     theme(text = element_text(size = 16),
#           axis.text.x = element_text(vjust = 1,
#                                      hjust = 0),
#           axis.text.y = element_text(margin = margin(0,0,0,0)),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.major.x = element_blank(),
#           axis.ticks = element_blank()) +
#     # Clip
#     coord_cartesian(expand = F) %>% return()
# }
# 
# # Plot and save
# plot_annual_side_by_side_bloom_pptx(annual_freq_lake_weekly,annual_freq_lake_bloom)
# 
# ggsave(
#   filename = 'annual_freq_lake_plot_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 6.5,
#   height = 5,
#   units = 'in'
# )
# 
# plot_monthly_side_by_side_bloom_pptx(monthly_freq_lake_weekly,monthly_freq_lake_bloom)
# 
# ggsave(
#   filename = 'monthly_freq_lake_plot_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 6.5,
#   height = 5,
#   units = 'in'
# )
# 
# Annual quadrants ----
# plot_annual_side_by_side_bloom_pptx(annual_freq_quads_weekly,annual_freq_quads_bloom) + 
#   facet_wrap(~quad) +
#   theme(
#     axis.text.x = element_text(hjust = .2),
#     strip.background = element_rect(fill = 'transparent'),
#     panel.spacing = unit(1,'line')
#   )
# 
# ggsave(
#   filename = 'annual_freq_quads_plot_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 11.5,
#   height = 5.5,
#   units = 'in'
# )

# Monthly quadrants ----

# plot_monthly_side_by_side_bloom_pptx(monthly_freq_quads_weekly,monthly_freq_quads_bloom) + 
#   facet_wrap(~quad) +
#   theme(
#     axis.text.x = element_text(vjust = 1,hjust = .15),
#     strip.background = element_rect(fill = 'transparent'),
#     panel.spacing = unit(1,'line')
#   )
# 
# ggsave(
#   filename = 'monthly_freq_quads_plot_pptx.png',
#   path = '/Users/natalie/Documents/thesis_writing/r_plots',
#   width = 11.5,
#   height = 5.5,
#   units = 'in'
# )


# Generate side-by-side tables for reference when writing ----
annual_freq_lake_weekly %>% 
  dplyr::select(`WY(date)`, freq_2, freq_3) %>%
  mutate(freq_range = freq_2 + freq_3) %>%
  left_join(dplyr::select(annual_freq_lake_bloom,`WY(date)`, freq_bloom)) -> annual_summary_tib_lake

annual_freq_quads_weekly %>%
  dplyr::select(`WY(date)`, freq_2, freq_3,quad) %>%
  mutate(freq_range = freq_2 + freq_3) %>%
  left_join(dplyr::select(annual_freq_quads_bloom,`WY(date)`,quad,freq_bloom)) %>% 
  mutate(quad = factor(quad, levels = c('Northwest Quadrant','Northeast Quadrant','Southwest Quadrant','Southeast Quadrant'),
                       labels = c('nw_quad','ne_quad','sw_quad','se_quad'),ordered = T)) -> annual_summary_tib_quads

monthly_freq_lake_weekly %>% 
  dplyr::select(`month(date)`, freq_2, freq_3) %>%
  mutate(freq_range = freq_2 + freq_3) %>%
  left_join(dplyr::select(monthly_freq_lake_bloom,`month(date)`, freq_bloom)) -> monthly_summary_tib_lake

monthly_freq_quads_weekly %>% 
  dplyr::select(`month(date)`, freq_2, freq_3,quad) %>%
  mutate(freq_range = freq_2 + freq_3) %>%
  left_join(dplyr::select(monthly_freq_quads_bloom,`month(date)`,quad,freq_bloom)) %>% 
  mutate(quad = factor(quad, levels = c('Northwest Quadrant','Northeast Quadrant','Southwest Quadrant','Southeast Quadrant'),
                       labels = c('nw_quad','ne_quad','sw_quad','se_quad'),ordered = T))  -> monthly_summary_tib_quads

# ALL ANNUAL FREQUENCIES
annual_summary_tib_lake %>%
  mutate(sb = 'lake') %>%
  bind_rows(annual_summary_tib_quads %>% rename(sb = quad)) -> annual_summary_tib

# ALL MONTHLY FREQUENCIES
monthly_summary_tib_lake %>%
  mutate(sb = 'lake') %>%
  bind_rows(monthly_summary_tib_quads %>% rename(sb = quad)) -> monthly_summary_tib

# Correlations ----
tibble(
  spatial_boundary = rep(c('lake','nw_quad','ne_quad','sw_quad','se_quad'),2),
  time_group = c(rep('annual',5),rep('monthly',5)),
  tau = NA,
  p_value = NA,
  cohens_d = NA
) -> cor_tibble

cor_from_summary_tibs <- function(spatial_boundary,time_group) {
  if(time_group == 'annual') {
    tib <- annual_summary_tib
  } else {tib <- monthly_summary_tib}
  
  tib %>% filter(sb == spatial_boundary) %>% pull(freq_range) -> freq_sim
  tib %>% filter(sb == spatial_boundary) %>% pull(freq_bloom) -> freq_obs
  
  cor_vals <- cor.test(x = freq_sim,y = freq_obs,method = 'kendall',alternative = 'two.sided')
  
  effect_size <- cohens_d(x = freq_sim,y = freq_obs,paired = T,alternative = 'two.sided')
  
  return(c(cor_vals$estimate[[1]],cor_vals$p.value,effect_size$Cohens_d))
}

for(i in 1:nrow(cor_tibble)) {
  vals_i <- cor_from_summary_tibs(spatial_boundary = cor_tibble$spatial_boundary[i],
                                  time_group = cor_tibble$time_group[i])
  cor_tibble$tau[i] <- vals_i[1]
  cor_tibble$p_value[i] <- vals_i[2]
  cor_tibble$cohens_d[i] <- vals_i[3]
}

annual_summary_tib %>% filter(sb == 'lake') %>% pull(freq_range)

cohens_d(x = freq_range,
         y = freq_bloom,
         data = filter(annual_summary_tib,sb == 'lake'),
         alternative = 'two.sided',
         paired = T) 

cohens_d(x = (annual_summary_tib %>% filter(sb == 'lake') %>% pull(freq_range)),
            y = (annual_summary_tib %>% filter(sb == 'lake') %>% pull(freq_bloom)),
            alternative = 'two.sided', paired = T) ->d

boxplot((annual_summary_tib %>% filter(sb == 'lake') %>% pull(freq_range)),
        (annual_summary_tib %>% filter(sb == 'lake') %>% pull(freq_bloom)))

boxplot((monthly_summary_tib %>% filter(sb == 'lake') %>% pull(freq_range)),
        (monthly_summary_tib %>% filter(sb == 'lake') %>% pull(freq_bloom)))

# Save workspace -- last saved 4/28 after re-calculating correlations (no change lol) and re-making figs
# save.image("~/Documents/R projects/olci_ci_copy/frequency_analysis_in_situ_workspace_COPY.RData")

# Look at S-308 chla values - combined plot ----

hydro <- roi_tile_vals %>%
  dplyr::select(month,year,total_pot_bloom_freq) %>%  # year references WY
  full_join(s308_usgs_filtered) %>%
  ggplot() +
  geom_hline(yintercept = 0,color = 'gray70') +
  geom_line(aes(x = date,y = avg_7d)) +
  geom_rug(aes(x = date,color = total_pot_bloom_freq),
           length = unit(0.04, "npc"),
           sides = 'b') +
  facet_grid(col = vars(year),scales = 'free_x',
             labeller = labeller(year = WY_labs)) +
  xlab('Month') + ylab('Flow (cfs)') +
  scale_x_date(breaks = seq(from=as.Date("2016-05-01"),
                            to=as.Date("2021-05-01"),
                            by="month"),
               minor_breaks = NULL,
               labels = c(year_labels,'M')) +
  scale_color_colormap(aesthetics = 'color', 
                       name = 'Total potential\nbloom frequency',
                       discrete = F,
                       colormap = colormaps$salinity,
                       labels = label_percent(accuracy = 1,
                                              scale = 1)) +
  theme_minimal() +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.text.x = element_text(hjust = -0.3,vjust = 1),
        strip.text.x = element_blank(),
        text = element_text(size = 12)) 


ggarrange(
  ggplot(data = (s308c_chla %>%
           mutate(date = date(date),
                  month = month(date),
                  year = WY(date)) %>% 
             rename(chla = `S308C_CHLOROPHYLL-A(LC)_ug/L`)%>%
             full_join(s308_usgs_filtered))) +
    geom_line(aes(x = date,y = avg_7d*0.001),color = 'transparent') +
    geom_hline(yintercept = 40,linetype = 'longdash',
               color = 'seagreen') +
    geom_hline(yintercept = 0,color = 'gray70') +
    geom_point(aes(x = date,y = chla),alpha = 0.8) + 
    facet_grid(col = vars(year),scales = 'free_x',
               labeller = labeller(year = WY_labs)) + 
    scale_x_date(breaks = seq(from=as.Date("2016-05-01"),
                              to=as.Date("2021-05-01"),
                              by="month"),
                 minor_breaks = NULL,
                 labels = c(year_labels,'M')) +
    scale_y_continuous(name = expression(paste('Clorophyll-',italic('a'),' (µg/L)')),
                       sec.axis = sec_axis(trans = ~. / 0.01, 
                                           name = '')) +
    theme_minimal() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          text = element_text(size = 12)),
  
  hydro,
  
  ncol = 1, nrow = 2,common.legend = T,legend = 'right',heights = c(1,2),align = 'v'
)

s308c_chla %>%
  mutate(date = date(date),
         month = month(date),
         year = WY(date)) %>% 
  rename(chla = `S308C_CHLOROPHYLL-A(LC)_ug/L`)%>%
  left_join(s308_usgs_filtered) %>% filter(chla >= 40) 

s308c_chla %>%
  mutate(date = date(date),
         month = month(date),
         year = WY(date)) %>% 
  rename(chla = `S308C_CHLOROPHYLL-A(LC)_ug/L`)%>%
  left_join(s308_usgs_filtered) %>% View()

ggsave(
  filename = 'hydrograph_word_with_chla_4.png',
  path = '/Users/natalie/Documents/thesis_writing/r_plots',
  width = 9,
  height = 5,
  units = 'in'
)

# I think everything from here down was experimental ----
roi_tile_vals %>%
  rename(
    month = `month(date)`,
    year = `WY(date)`,
  ) %>%
  mutate(total_bloom = total_2 + total_3,
         bloom_freq = total_bloom/(count_pixels(roi)*n)*100) %>%
  dplyr::select(month,year,bloom_freq) %>%  # year references WY
  full_join(s308_usgs_filtered) %>%
  full_join(s308c_chla %>%
              mutate(date = date(date),
                     month = month(date),
                     year = WY(date)) %>% 
              rename(chla = `S308C_CHLOROPHYLL-A(LC)_ug/L`)) %>% 
  gather(avg_7d,chla,key = 'data', value = 'value') %>%
  ggplot() + 
  geom_hline(yintercept = 0,color = 'red',linetype = 'twodash') +
  geom_hline(yintercept = 40, color = 'seagreen',linetype = 'longdash') + 
  geom_point(aes(x = date, 
                 y = chla),
             size = 2,
             alpha = 0.7, color = 'blue') +
  geom_line(aes(x = date,y = avg_7d*.1),
            color = 'black', size = .75) + 
  geom_rug(aes(x = date,color = bloom_freq),
           sides = 'b') +

  facet_grid(col = vars(year),scales = 'free_x',
             labeller = labeller(year = WY_labs)) +
  scale_y_continuous(name = 'Clorophyll-a concentation (µg/L)',
                     sec.axis = sec_axis(trans = ~. / 0.1, 
                                         name = 'Discharge to SLC (cfs)')) +
  xlab('Month') + ylab('Flow (cfs)') +
  scale_x_date(breaks = seq(from=as.Date("2016-05-01"),
                            to=as.Date("2021-05-01"),
                            by="month"),
               minor_breaks = NULL,
               labels = c(year_labels,'M')) +
  scale_color_colormap(aesthetics = 'color', 
                       name = 'Total potential\nbloom frequency',
                       discrete = F,
                       colormap = colormaps$salinity,
                       limits = c(0,50)) +
  theme_minimal() +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.text.x = element_text(hjust = -0.3,vjust = 1),
        text = element_text(size = 12))



roi_tile_vals %>%
  rename(
    month = `month(date)`,
    year = `WY(date)`,
  ) %>%
  mutate(total_bloom = total_2 + total_3,
         bloom_freq = total_bloom/(count_pixels(roi)*n)*100) %>%
  dplyr::select(month,year,bloom_freq) %>%  # year references WY
full_join(s308_usgs_filtered) %>%
  ggplot() +
  geom_hline(yintercept = 0,color = 'red',linetype = 'twodash') +
  geom_line(aes(x = date,y = avg_7d)) +
  geom_rug(aes(x = date,color = bloom_freq),
           sides = 'b') +
  facet_grid(col = vars(year),scales = 'free_x',
             labeller = labeller(year = WY_labs)) +
  xlab('Month') + ylab('Flow (cfs)') +
  scale_x_date(breaks = seq(from=as.Date("2016-05-01"),
                            to=as.Date("2021-05-01"),
                            by="month"),
               minor_breaks = NULL,
               labels = c(year_labels,'M')) +
  scale_color_colormap(aesthetics = 'color', 
                       name = 'Total potential\nbloom frequency',
                       discrete = F,
                       colormap = colormaps$salinity,
                       limits = c(0,50)) +
  theme_minimal() +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.text.x = element_text(hjust = -0.3,vjust = 1),
        text = element_text(size = 12))

# plot to visually verify
annual_summary_tib %>% 
  ggplot() +  
  geom_abline(slope = 1,intercept = 0,color = 'red') +
  geom_smooth(aes(x = freq_range,y = freq_bloom),method = 'lm') +
  geom_point(aes(x = freq_range,y = freq_bloom)) +
  facet_wrap(~sb) +
  theme_bw() +
  xlab('frequencies from satellite') +
  ylab('frequencies from field data') +
  xlim(0,30) + ylim(0,30) 

monthly_summary_tib %>% 
  ggplot() +  
  geom_abline(slope = 1,intercept = 0,color = 'red') +
  geom_smooth(aes(x = freq_range,y = freq_bloom),method = 'lm') +
  geom_point(aes(x = freq_range,y = freq_bloom)) +
  facet_wrap(~sb) +
  theme_bw() +
  xlab('frequencies from satellite') +
  ylab('frequencies from field data') +
  xlim(0,40) + ylim(0,40) 

