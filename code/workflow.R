# workflow.R


# Libraries ---------------------------------------------------------------

# Load all meta-data and functions
source("code/functions.R")

# Remove scientific notation from results
# options(scipen=999)


# Reference analysis ------------------------------------------------------

# Combine the three reference time series, run analysis, and save
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
#   random_results <- plyr::ldply(1:100, random_analysis, .parallel = T)
# ) # 417 seconds
# saveRDS(random_results, "data/random_results_100.Rda")

# Calculate the full analysis on 1000 random pixels
# doMC::registerDoMC(cores = 25) # NB: 50 cores uses too much RAM
# set.seed(666)
# system.time(
#   random_results <- plyr::ldply(1:1000, random_analysis, .parallel = T)
# ) # 417 seconds
# saveRDS(random_results, "data/random_results_1000.Rda")


# Global analysis ---------------------------------------------------------

# Wrapper function with a chatty output as it chugs along
# global_analysis_single <- function(file_sub, par_op = F){
#   OISST_slice <- OISST_files[file_sub]
#   lon_row_pad <- str_pad(file_sub, width = 4, pad = "0", side = "left")
#   print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))
#   slice_res <- global_analysis(OISST_slice, par_op = par_op)
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
# plyr::l_ply(missing_index, global_analysis_single, .parallel = F, par_op = T)


# Combine global results --------------------------------------------------

# Load and combine each longitude slice of results
# NB: Uses too much RAM to load everything in one shot...
# NB: Best to run this in a terminal and not RStudio
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

# Save
# saveRDS(global_results, "data/global_results.Rda")


# Global trends -----------------------------------------------------------
# Calculate the simple linear trends for the different tests at each pixel

# Set cores
# doMC::registerDoMC(cores = 50)

# Load global data
# system.time(
#   global_results <- readRDS("data/global_results.Rda")
# ) # 665 seconds

# Calculate trends and save
# NB: Not run as it takes too long
# system.time(
#   global_trend <- plyr::ddply(global_results, .variables = c("lat"),
#                               .fun = pixel_trend, .parallel = T)
# ) # xxx seconds
# saveRDS(global_trend, "data/global_trend.Rda")

# Filter down to more immediately necessary results
# system.time(
#   global_mean_perc <- global_results %>%
#     filter(id %in% c("mean_perc", "n_diff"),
#            var %in% c("count", "focus_count",
#                       "duration", "focus_duration",
#                       "intensity_max", "focus_intensity_max"))
# ) # 91 seconds

# Calculate trends and save
# NB: Not run in favour of the following chunk
# system.time(
#   global_mean_perc_trend <- plyr::ddply(global_mean_perc, .variables = c("lat"),
#                                         .fun = pixel_trend, .parallel = T)
# ) # 59 seconds for one lon slice,
# saveRDS(global_mean_perc_trend, "data/global_mean_perc_trend.Rda")

# system.time(
#   global_focus_trend <- plyr::ldply(dir("data/global", full.names = T),
#                                         .fun = focus_trend, .parallel = T)
# ) # 41 seconds for one lon slice, 44 minutes for all
# saveRDS(global_focus_trend, "data/global_focus_trend.Rda")


# Figures -----------------------------------------------------------------

# All of the figures that use the above results are made in "code/figures.R"


# Results -----------------------------------------------------------------
# Code that generates numeric results referred to in the text outside of figures

# Upper and lower quantile ranges for sub-optimal tests
# NB: `full_results` created in Figure 2 section of code/figures.R
bad_pixels <- c("140.375_0.625", "-73.625_-77.125")
quant_subopt <- full_results %>%
  filter(var %in% c("count", "duration", "intensity_max",
                    "focus_count", "focus_duration", "focus_intensity_max"),
         id %in% c("n_diff", "mean_perc", "sum_perc", "mean_perc"),
         !(site %in% bad_pixels)) %>%
  group_by(test, index_vals, var, id) %>%
  summarise(lower = quantile(val, 0.05),
            upper = quantile(val, 0.95)) %>%
  ungroup()
quant_length <- filter(quant_subopt, test == "length")
quant_miss <- filter(quant_subopt, test == "missing")
quant_trend <- filter(quant_subopt, test == "trend")



# Appendix 1 --------------------------------------------------------------

# The effect of the sub-optimal tests on seas/thresh


# Appendix 2 --------------------------------------------------------------
# The difference between the proper 30 year base period and all other 30 year base periods

# Combine the three reference time series, run analysis, and save
# sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
#                  mutate(sst_NW_Atl, site = "NW_Atl"),
#                  mutate(sst_Med, site = "Med"))
# system.time(
#   sst_ALL_results <- plyr::ddply(sst_ALL, c("site"), base_period_analysis, .parallel = T, clim_metric = T)
# ) # 3 seconds
# saveRDS(sst_ALL_results, "data/sst_ALL_bp_results.Rda")

# Calculate the base period analysis on 100 random pixels
# doMC::registerDoMC(cores = 25) # NB: 50 cores uses too much RAM
# set.seed(666)
# system.time(
#   random_results <- plyr::ldply(1:100, random_analysis, .parallel = T, base_period = T)
# ) # 53 seconds
# saveRDS(random_results, "data/random_bp_results_100.Rda")


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
