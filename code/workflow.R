# workflow.R


# Libraries ---------------------------------------------------------------

# Load all meta-data and functions
source("code/functions.R")

# Remove scientific notation from results
# options(scipen=999)


# Reference analysis ------------------------------------------------------

# Combine the three reference ts, run results, and save
# sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
#                  mutate(sst_NW_Atl, site = "NW_Atl"),
#                  mutate(sst_Med, site = "Med"))
# system.time(
#   sst_ALL_results <- plyr::ddply(sst_ALL, c("site"), single_analysis, .parallel = T,
#                                  full_seq = T, clim_metric = T, count_miss = T, windows = T)
# ) # 72 seconds
# saveRDS(sst_ALL_results, "data/sst_ALL_results.Rda")


# Random analysis ---------------------------------------------------------

# Calculate the full analysis on 100 random pixels
# doMC::registerDoMC(cores = 25) # NB: 50 cores uses too much RAM
# set.seed(666)
# system.time(
#   random_results <- plyr::ldply(1:100, random_analysis, .parallel = T, .id = "random")
# ) # 417 seconds
# saveRDS(random_results, "data/random_results.Rda")


# Global analysis ---------------------------------------------------------

# Wrapper function with a chatty output as it chugs along
# global_analysis_single <- function(file_sub){
#   OISST_slice <- OISST_files[file_sub]
#   lon_row_pad <- str_pad(file_sub, width = 4, pad = "0", side = "left")
#   print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))
#   slice_res <- global_analysis(OISST_slice)
#   saveRDS(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rda"))
#   print(paste0("Finished run on step ",lon_row_pad," at ",Sys.time()))
#   rm(slice_res); gc()
# }
# plyr::l_ply(1:1440, global_analysis_single, .parallel = T) # This took 37.5 hours to run

# The nightly running of the MHW Tracker seems to have interfered with the
# calculation of several lon slices
# full_files <- paste0("slice_",str_pad(1:1440, width = 4, pad = "0", side = "left"),".Rda")
# missing_files <- setdiff(full_files, dir("data/global"))
# missing_index <- as.numeric(sapply(str_split(sapply(str_split(missing_files, "_"),
#                                                     "[[", 2), "[.]"), "[[", 1))
# plyr::l_ply(missing_index, global_analysis_single, .parallel = T)


# Combine global results --------------------------------------------------

# Load and combine each longitude slice of results
# NB: Uses too much RAM to load everything in one shot...
# system.time(
#   global_results_1 <- plyr::ldply(dir("data/global", full.names = T)[1:500], readRDS, .parallel = T)
# ) # 184 seconds
# system.time(
#   global_results_2 <- plyr::ldply(dir("data/global", full.names = T)[501:1000], readRDS, .parallel = T)
# ) # 312 seconds
# global_results_1_2 <- rbind(global_results_1, global_results_2)
# rm(global_results_1, global_results_2); gc()
# system.time(
#   global_results_3 <- plyr::ldply(dir("data/global", full.names = T)[1001:1440], readRDS, .parallel = T)
# ) # xxx seconds
# global_results <- rbind(global_results_1_2, global_results_3)
# rm(global_results_1_2, global_results_3); gc()
#
# # Save
# saveRDS(global_results, "data/global_results.Rda")


# Process results ---------------------------------------------------------

# Calculate the simple linear slopes for the different tests at each pixel

# Set cores
# doMC::registerDoMC(cores = 50)

# Load global data
# global_summary <- map_dfr(readRDS, dir("data/global", full.names = T))

# Calculate summary slopes
# global_summary_slope <- plyr::ddply(global_summary, .variables = c("lat"),
#                                   .fun = global_slope, .parallel = T)
# saveRDS(global_summary_slope, "data/global_summary_slope.Rda")

# Load focus data
# global_focus <- readRDS("data/global_focus.Rda")

# Calculate focus slopes
# global_focus_slope <- plyr::ddply(global_focus, .variables = c("lat"),
#                                   .fun = global_slope, .parallel = T)
# saveRDS(global_focus_slope, "data/global_focus_slope.Rda")


# Global relationships ----------------------------------------------------

# In this section we look at the relationships between certain global variables

# Set cores
# doMC::registerDoMC(cores = 50)

# Global decadal trends
# load("data/global_dec_trend.Rdata")

# The effect on single events
# load("data/global_effect_event.Rdata")

# The slope results for single events
# load("data/global_effect_event_slope.Rdata")

# Merge and filter data
# event_slope_dec_trend <- left_join(global_effect_event_slope, global_dec_trend, by = c("lon", "lat")) %>%
#   filter(test == "length") %>%
#   na.omit() %>%
#   droplevels()

# Linear model of relationship
# event_slope_dec_trend_model <- event_slope_dec_trend %>%
#   group_by(test, metric) %>%
#   do(model = broom::glance(lm(slope ~ dec_trend, data = .))) %>%
#   unnest()

# A regression scatterplot between decadal trend and duration/max.int. and ts length slope
# event_dec_scatterplot <- ggplot(data = filter(event_slope_dec_trend,
#                                               metric %in% c("duration", "intensity_max")),
#                                 aes(x = dec_trend, y = slope)) +
#   geom_point(aes(colour = metric)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~metric, scales = "free_y", ncol = 1)
# event_dec_scatterplot

# The spread of the results per year measured
# effect_event_spread <- global_effect_event %>%
#   filter(test == "length") %>%
#   droplevels() %>%
#   spread(key = index_vals, value = val) %>%
#   mutate(val_spread = `30`-`10`)

# Plot showing value spread between 30 and 10 years
# duration_spread_plot <- ggplot(filter(effect_event_spread, metric == "intensity_max"),
#                                aes(x = lon, y = lat)) +
#   geom_raster(aes(fill = val_spread)) +
#   geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
#   scale_fill_gradient2(low = "blue", high = "red") +
#   coord_equal(expand = F) +
#   # labs(fill = "Linear trend (°C/dec)") +
#   theme_void() +
#   theme(legend.position = "bottom",
#         legend.key.width = unit(2, "cm"))
# duration_spread_plot

# Join in the event metrics detected from full 30 year time series
# event_slope_prop <- left_join(event_slope_dec_trend, effect_event_spread,
#                               by = c("lon", "lat", "test", "metric")) %>%
#   na.omit() %>%
#   mutate(slope_prop_10 = round(slope/`10`, 3),
#          slope_prop_20 = round(slope/`20`, 3),
#          slope_prop_30 = round(slope/`30`, 3))

# Get just slope ten and change it for easy plotting
# event_slope_prop_10 <- event_slope_prop %>%
#   dplyr::select(lat:metric, slope_prop_10) %>%
#   dplyr::rename(slope = slope_prop_10)

# # Global map showing patterns of the slope as a proportion of the overall change
# fig_4 <- global_effect_event_slope_plot(test_sub = "length", metric_sub = "intensity_max",
#                                         df = event_slope_prop_10, prop = T)
# # fig_4
# ggsave(plot = fig_4, filename = "LaTeX/fig_4.pdf", height = 6, width = 10)
# ggsave(plot = fig_4, filename = "LaTeX/fig_4.png", height = 6, width = 10)
# ggsave(plot = fig_4, filename = "LaTeX/fig_4.jpg", height = 6, width = 10)
#
# fig_5 <- global_effect_event_slope_plot(test_sub = "length", metric_sub = "duration",
#                                         df = event_slope_prop_10, prop = T)
# # fig_5
# ggsave(plot = fig_5, filename = "LaTeX/fig_5.pdf", height = 6, width = 10)
# ggsave(plot = fig_5, filename = "LaTeX/fig_5.png", height = 6, width = 10)
# ggsave(plot = fig_5, filename = "LaTeX/fig_5.jpg", height = 6, width = 10)


# Boxplot showing spread of proportion values


# Linear model of relationship between full event metric and rate of change from shortening
# event_slope_prop_model <- event_slope_prop %>%
#   group_by(test, metric) %>%
#   filter(metric %in% c("duration", "intensity_max")) %>%
#   do(model = broom::glance(lm(slope ~ dec_trend, data = .))) %>%
#   unnest()

# A regression scatterplot between decadal trend and duration/max.int. and ts length slope
# event_slope_prop_scatterplot <- ggplot(data = filter(event_slope_prop,
#                                                      metric %in% c("duration", "intensity_max")),
#                                        aes(x = val, y = slope)) +
#   geom_point(aes(colour = metric)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~metric, scales = "free", ncol = 1)
# event_slope_prop_scatterplot


# More thoughts -----------------------------------------------------------

# It's not reasonable to run a decadal trend test in addition to the missing
# and length tests because this value is being added in a completely
# controlled way.
# Meaning that controlling for it is the same as simply not adding it
# in the first place.
# To that end we want to provide a post-hoc correction.
# Because we have seen that max. intensity is a function of the
# decadal trend and the peak date of the event we should be able to
# correct for intensities by subtracting the decadal trend multiplied
# by the peak date as it relates to the length of the time series.
# For example, if the decadal trend is 0.3C/dec, and the peak date of an event
# is in the 25th year of a 30 year time series (0.8333 of the length),
# then the the impact of the decadal trend on the max. intensity should be
# 0.3*0.83 = 0.25C
# So one would subtract this value in order to correct the results to
# match a de-trended time series.
# The problem then becomes, why would anyone actually want to do this?
# No. Rather we must provide a post-hoc fix for the potential impact of
# a decadal trend on MHWs in conjunction with short time series.
# And then on top of that add in the correction for missing data.
# Mercifully the correction for missing data is very simple and should
# play nice with the other issues.
# This is because the linear interpolation of NA gaps will provide
# data that are matching the temperatures around them so the problem of
# having more missing data on one end of a time series over the other
# should not be pronounced.

# To correct for the effect of time series length on event metrics we need
# two things:
# The length of the time series
# The decadal trend in the region
# It may be that rather than correcting the value, we should provide a CI instead

# The relationship may be better aided by the known variance in the time series

# Must see what the R2 + SE is between decadal trend and change in duration/max.int.

# Where on the x axis things go wrong is the main question to be answers

# Show the difference in the moving 30 year clim vs. the preferred 30 year clim as an appendix figure

# It may end up being best to offer advise based on the change in MHW count/days
# Also the proportion shift in duration and max int based on something bio relevant in literature
