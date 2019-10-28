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
#                                  full_seq = T, clim_metric = F, count_miss = T, windows = T)
# ) # 65 seconds
# saveRDS(sst_ALL_results, "data/sst_ALL_results.Rda")


# Random analysis ---------------------------------------------------------

# Calculate the full analysis on 100 random pixels
# doMC::registerDoMC(cores = 50)
# set.seed(666)
# system.time(
#   random_results <- plyr::ldply(1:100, random_analysis, .parallel = T)
# ) # 262 seconds
# saveRDS(random_results, "data/random_results_100.Rda")

# Calculate the full analysis on 1000 random pixels
# doMC::registerDoMC(cores = 25) # 50 appears to be too much
# set.seed(666)
# system.time(
#   random_results <- plyr::ldply(1:1000, random_analysis, .parallel = T)
# ) # ~54 minutes
# saveRDS(random_results, "data/random_results_1000.Rda")

# A couple of cores slipped during the processing of this file
# random_results <- readRDS("data/random_results_1000.Rda")
# length(unique(unite(random_results, lon, lat, col = "site")$site)) # 920 unique time series
# It is necessary to add 80 more time series to this file
# doMC::registerDoMC(cores = 40)
# set.seed(999)
# system.time(
#   random_results_patch <- plyr::ldply(1:80, random_analysis, .parallel = T)
# ) # 226 seconds
# random_results_patched <- rbind(random_results, random_results_patch)
# saveRDS(random_results_patched, "data/random_results_1000.Rda")


# Why do some MHWs dissapear from wider windows? --------------------------

# -112.125 -28.875 # A pixel negatively affected by window widening
# which(c(seq(0.125, 179.875, by = 0.25), seq(-179.875, -0.125, by = 0.25)) == -112.125)
# sst <- load_noice_OISST(OISST_files[992]) %>%
#   filter(lat == -28.875)

# Detrend the selected ts
# sst_flat <- detrend(sst)

# Calculate MHWs from detrended ts
# sst_flat_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2018-12-31")))

# Pull out the largest event in the ts
# focus_event <- sst_flat_MHW$event %>%
#   filter(date_start >= "2009-01-01") %>%
#   filter(intensity_cumulative == max(intensity_cumulative)) %>%
#   select(event_no, date_start:date_end, duration, intensity_cumulative, intensity_max) %>%
#   mutate(intensity_cumulative = round(intensity_cumulative, 2),
#          intensity_max = round(intensity_max, 2))

# Quickly visualise the largest heatwave in the last 10 years of data
# heatwaveR::event_line(sst_flat_MHW, start_date = "2009-01-01", metric = "intensity_cumulative")

# Normal window width
# window_5_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2018-12-31")))
# heatwaveR::event_line(window_5_MHW, start_date = "2009-01-01", metric = "intensity_cumulative")

# 10 day window
  # Already here we see why the event falls away
  # The focus MHW was just staying above the down slope of the seasonal dive into winter
  # When the window half width is expanded the seasonal decline becomes less steep and the
  # observed temperature is no longer above the 90th percentile
# window_10_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2018-12-31"), windowHalfWidth = 10))
# heatwaveR::event_line(window_10_MHW, start_date = "2009-01-01", metric = "intensity_cumulative")

# 20 day window
# window_20_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2018-12-31"), windowHalfWidth = 20))
# heatwaveR::event_line(window_20_MHW, start_date = "2009-01-01", metric = "intensity_cumulative")

# 30 day window
# window_30_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2018-12-31"), windowHalfWidth = 30))
# heatwaveR::event_line(window_30_MHW, start_date = "2009-01-01", metric = "intensity_cumulative")

# Now let's have a peak at each step along the way, just for laughs
# ts2clm_window <- function(window_choice, df = sst_flat){
#   res <- ts2clm(df, climatologyPeriod = c("1982-01-01", "2018-12-31"), windowHalfWidth = window_choice) %>%
#     mutate(site_label = paste0("window_",window_choice))
#   return(res)
# }

# Calculate clims
# sst_clim <- plyr::ldply(seq(5, 30, by = 5), ts2clm_window, .parallel = T)

# Climatologies doy
# sst_clim_only <- sst_clim %>%
#   select(-t, -temp) %>%
#   unique()

# Calculate events
# sst_event <- sst_clim %>%
#   group_by(site_label) %>%
#   group_modify(~detect_event(.x)$event)

# Find largest event in most recent ten years of data
# focus_event <- sst_event %>%
#   filter(date_start >= "2009-01-01") %>%
#   group_by(site_label) %>%
#   filter(intensity_cumulative == max(intensity_cumulative)) %>%
#   ungroup()

# Merge with results for better plotting
# sst_focus <- left_join(sst_clim,
#                        focus_event[,c("site_label", "date_start", "date_peak", "date_end")], by = "site_label") %>%
#   mutate(site_label = factor(site_label, levels = c("window_5", "window_10", "window_15",
#                                                     "window_20", "window_25", "window_30")))

# trend_fig <- fig_1_plot(sst_focus, spread = 150)
# trend_fig

# Look at differences between the seas/thresh for each window
# sst_clim_only %>%
#   select(-doy) %>%
#   gather(key = "var", value = "val", seas, thresh) %>%
#   group_by(site_label, var) %>%
#   summarise_if(.predicate = is.numeric, .funs = c("min", "median", "mean", "max")) %>%
#   ungroup() %>%
#   gather(key = "stat", value = "val", -site_label, - var) %>%
#   mutate(site_label = factor(site_label, levels = c("window_5", "window_10", "window_15",
#                                                     "window_20", "window_25", "window_30"))) %>%
#   arrange(site_label) %>%
#   ggplot(aes(x = stat, y = val, colour = site_label)) +
#   geom_point() +
#   scale_colour_brewer() +
#   facet_wrap(~var)

# Now let's look at all of the 100 random results to see how this shakes out
# random_results <- readRDS("data/random_results_100.Rda")
# unique(random_results$test)
# all_clims <- random_results %>%
#   filter(test %in% c("length", "window_10", "window_20", "window_30"),
#          index_vals == 30,
#          var %in% c("seas", "thresh"),
#          id %in% c("min", "median", "mean", "max", "sd")) %>%
#   ggplot(aes(x = id, y = val, fill = test)) +
#   geom_boxplot() +
#   scale_fill_brewer(palette = "YlOrRd") +
#   facet_wrap(~var)
# all_clims


# Why are the buil-in time series so anomalous? ---------------------------

# We'll use the WA time series and the pixel just adjacent to it
# sst_WA_flat <- detrend(sst_WA) %>%
#   mutate(site = "WA") %>%
#   select(site, t, temp)

# which(c(seq(0.125, 179.875, by = 0.25), seq(-179.875, -0.125, by = 0.25)) == 112.375)
# sst_flat <- load_noice_OISST(OISST_files[450]) %>%
#   filter(lat > -30, lat < -20) %>%
#   unite(lon, lat, col = "site") %>%
#   group_by(site) %>%
#   group_modify(~detrend(.x)) %>%
#   data.frame()

# sst_ALL <- rbind(sst_flat, sst_WA_flat)
# unique(sst_ALL$site)

# Plot all time series together
# sst_ALL %>%
#   filter(t >= "2010-01-01", t <= "2012-12-31") %>%
#   ggplot(aes(x = t, y = temp)) +
#   geom_line(aes(group = site, colour = site)) +
#   scale_colour_viridis_d()

# Run full analyses on both
# result_ALL <- plyr::ddply(sst_ALL, c("site"), single_analysis, full_seq = T, .parallel = T)

# fig_box_plot(result_ALL, tests = "base", result_choice = "10_years")
# fig_box_plot(result_ALL, tests = "base", result_choice = "focus")

# These boxplots don't show us much, taking a different approach
# result_ALL %>%
#   filter(test == "trend", var == "duration", id == "sum_perc") %>%
#   ggplot(aes(x = index_vals, y = val)) +
#   geom_line(aes(group = site, colour = site)) +
#   scale_colour_viridis_d()
# In the figure above we may see that the closer we appraoch the centre of the WA MHW the less of an effect the
# increasing decadal trend is having on the overall number of MHWs detected
# We may deduce that this is because the WA MHW was so intense that it artificially raising up the 90th percentile
# so high that even with added decadal warming it is not enough to increase the other MHWs

# Now we want to look at how the count of overall events are affected
# result_ALL %>%
#   filter(test == "trend", var == "count", id == "n_perc") %>%
#   ggplot(aes(x = index_vals, y = val)) +
#   geom_line(aes(group = site, colour = site)) +
#   scale_colour_viridis_d()

# And there you have it, the built-in time series are just super wacky, their is nothing incorrect with the results


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
# plyr::l_ply(1:500, global_analysis_single, .parallel = T) # This took xxx hours to run

# The nightly running of the MHW Tracker seems to have interfered with the
# calculation of several lon slices
# full_files <- paste0("slice_",str_pad(1:1440, width = 4, pad = "0", side = "left"),".Rda")
# missing_files <- setdiff(full_files, dir("data/global"))
# missing_index <- as.numeric(sapply(str_split(sapply(str_split(missing_files, "_"),
#                                                     "[[", 2), "[.]"), "[[", 1))
# plyr::l_ply(missing_index, global_analysis_single, .parallel = F, par_op = T)


# Global trends -----------------------------------------------------------
# Calculate the simple linear trends for the different tests at each pixel

# Set cores
doMC::registerDoMC(cores = 50)

# Calculate trends and save
system.time(
  global_var_trend <- plyr::ldply(dir("data/global", full.names = T),
                                  .fun = var_trend, .parallel = T)
) # 60 seconds for one lon slice, ~70 minutes for all
saveRDS(global_var_trend, "data/global_var_trend.Rda")


# Figures -----------------------------------------------------------------

# All of the figures that use the above results are made in "code/figures.R"


# Results -----------------------------------------------------------------
# Code that generates numeric results referred to in the text outside of figures

# Upper and lower quantile ranges for sub-optimal tests
# NB: `full_results` created in Figure 2 section of code/figures.R
# bad_pixels <- c("140.375_0.625", "-73.625_-77.125")
# quant_subopt <- full_results %>%
#   filter(var %in% c("count", "duration", "intensity_max",
#                     "focus_count", "focus_duration", "focus_intensity_max"),
#          id %in% c("n_diff", "mean_perc", "sum_perc", "mean_perc"),
#          !(site %in% bad_pixels)) %>%
#   group_by(test, index_vals, var, id) %>%
#   summarise(lower = quantile(val, 0.05),
#             upper = quantile(val, 0.95)) %>%
#   ungroup()
# quant_length <- filter(quant_subopt, test == "length")
# quant_miss <- filter(quant_subopt, test == "missing")
# quant_trend <- filter(quant_subopt, test == "trend")


# Best practices ----------------------------------------------------------
# Code that produces the tables and supports the statements made in the Best Practices section

# Load the random 1000 data
# system.time(
#   random_results <- readRDS("data/random_results_1000.Rda") %>%
#     unite("site", c(lon, lat))
# ) # 68 seconds, 15 seconds without the "site" column

# Create table showing the rates of change int the results
# Find the 5th, 25th, 50th,75th, and 95th quantile at each step
# Fit linear models to those and provide those trends + R2 values
# Also show where the inflection points may be where the trends change
# var_choice <- data.frame(var = c("count", "duration", "intensity_max", "focus_count", "focus_duration", "focus_intensity_max"),
#                          id = c("n_perc", "sum_perc", "mean_perc", "mean_perc", "sum_perc", "mean_perc"),
#                          stringsAsFactors = F)
# random_quant <- random_results %>%
#   right_join(var_choice, by = c("var", "id")) %>%
#   mutate(test = as.character(test)) %>%
#   filter(test %in% c("length", "missing", "interp", "trend")) %>%
#   group_by(test, index_vals, var, id) %>%
#   summarise(q05 = quantile(val, 0.05),
#             q25 = quantile(val, 0.25),
#             q50 = quantile(val, 0.50),
#             q75 = quantile(val, 0.75),
#             q95 = quantile(val, 0.95),
#             iqr50 = q75-q25,
#             iqr90 = q95-q05) %>%
#   ungroup()

# Custom function for trend calculation
# Moved to functions.R

# Function that allows for plyr::ldply to iteratively look for best fits
  # NB: It assumes that random_quant above has been created and is in the environment
# testers...
# end_int = 3
# start_val = 30
# start_val = 0.25
# test_sub = "length"
# test_sub = "missing"
# rev_trend = F
# trend_test <- function(end_int = 3, test_sub = "length", start_val = 30, rev_trend = F){
#   df <- random_quant %>%
#     filter(test == test_sub)
#   df_index_vals <- unique(df$index_vals)
#   if(start_val == 30 & rev_trend == F){
#     end_val <- start_val-end_int
#     by_val <- -1
#   } else if(start_val == 30 & rev_trend == T){
#     end_val <- start_val+end_int
#     by_val <- 1
#   }else{
#     end_val <- start_val+(end_int*0.01)
#     by_val <- 0.01
#   }
#   df_model <- df %>%
#     filter(index_vals %in% seq(start_val, end_val, by = by_val)) %>%
#     gather(key = "stat", value = "val", -c(test:id)) %>%
#     group_by(test, var, id, stat) %>%
#     group_modify(~trend_stats(.x)) %>%
#     mutate(start_val = start_val,
#            end_val = end_val)
#   return(df_model)
# }

# Run the linear models at each possible step to deduce where any inflections points may be
# This is determined by tracking the change in R2 values, with lower values being bad
# quant_missing <- plyr::ldply(3:50, trend_test, .parallel = T, test_sub = "missing", start_val = 0)
# quant_missing_A <- plyr::ldply(3:25, trend_test, .parallel = T, test_sub = "missing", start_val = 0)
# quant_missing_B <- plyr::ldply(3:24, trend_test, .parallel = T, test_sub = "missing", start_val = 0.26)
# quant_interp <- plyr::ldply(3:50, trend_test, .parallel = T, test_sub = "interp", start_val = 0)
# quant_trend <- plyr::ldply(3:30, trend_test, .parallel = T, test_sub = "trend", start_val = 0)
# quant_length_A <- plyr::ldply(3:20, trend_test, .parallel = T, test_sub = "length")
# quant_length_B <- plyr::ldply(3:7, trend_test, .parallel = T, test_sub = "length", rev_trend = T)
# quant_ALL <- rbind(quant_missing_A, quant_missing_B, quant_interp, quant_trend, quant_length_A, quant_length_B)

## Test visuals to determine that the trends above are lekker
# First create a line plot of the results
# quant_ALL %>%
#   filter(test == "missing") %>%
#   ggplot(aes(x = end_val, y = r2)) +
#   geom_point(aes(colour = var)) +
#   geom_line(aes(colour = var)) +
#   facet_grid(stat~test, scales = "free_x")

# Filter out the trends that cover the correct ranges
# trend_filter <- data.frame(test = c("length", "length", "missing", "missing", "interp", "trend"),
#                            start_val = c(30, 30, 0, 0.26, 0, 0),
#                            end_val = c(10, 37, 0.25, 0.5, 0.5, 0.3))
# quant_filter <- quant_ALL %>%
#   right_join(trend_filter, by = c("test", "start_val", "end_val")) %>%
#   mutate(slope = ifelse(test == "length" & end_val == 10, -slope, slope),
#          end_point = end_val*slope) %>%
#   filter(stat != "iqr50", stat != "iqr90")

# Then project them onto the real data
# ggplot(random_quant) +
#   # 90 CI crossbars
#   # Need different lines for tests due to the different x-axis interval sizes
#   geom_crossbar(data = filter(random_quant, test == "length"),
#                 aes(x = index_vals, y = 0, ymin = q05, ymax = q95),
#                 fatten = 0, fill = "grey70", colour = NA, width = 1) +
#   geom_crossbar(data = filter(random_quant, test != "length"),
#                 aes(x = index_vals, y = 0, ymin = q05, ymax = q95),
#                 fatten = 0, fill = "grey70", colour = NA, width = 0.01) +
#   # IQR Crossbars
#   geom_crossbar(data = filter(random_quant, test == "length"),
#                 aes(x = index_vals, y = 0, ymin = q25, ymax = q75),
#                 fatten = 0, fill = "grey50", width = 1) +
#   geom_crossbar(data = filter(random_quant, test != "length"),
#                 aes(x = index_vals, y = 0, ymin = q25, ymax = q75),
#                 fatten = 0, fill = "grey50", width = 0.01) +
#   # Median segments
#   geom_crossbar(data = filter(random_quant, test == "length"),
#                 aes(x = index_vals, y = 0, ymin = q50, ymax = q50),
#                 fatten = 0, fill = NA, colour = "black", width = 1) +
#   geom_crossbar(data = filter(random_quant, test != "length"),
#                 aes(x = index_vals, y = 0, ymin = q50, ymax = q50),
#                 fatten = 0, fill = NA, colour = "black", width = 0.01) +
#   geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
#   geom_segment(data = quant_filter, aes(x = start_val, y = intercept, xend = end_val, yend = end_point, colour = stat)) +
#   scale_colour_brewer(palette = "Dark2") +
#   scale_x_continuous(expand = c(0, 0)) +
#   facet_grid(var ~ test, scales = "free", switch = "both") +
#   theme(legend.position = "top")

# These all look okay, but length is a little suspicious
# The slopes from 10 to 30 years are inverted, and the slopes from 30 to 37 years are too massive
# Best to pull these things out 1-by-1
# slope_length <- random_quant %>%
#   filter(test == "length", var == "intensity_max", index_vals %in% 30:37) %>%
#   dplyr::select(test:id, q95) %>%
#   dplyr::rename(val = q95) %>%
#   mutate(index_vals = abs(index_vals-30))
# slope_length <- trend_stats(slope_length) %>%
#   mutate(start_val = 30,
#          end_val = 37,
#          start_point = 0,
#          end_point = slope*7)


# Visualise manual test runs
# random_quant %>%
#   filter(test == "length", var == "intensity_max") %>%
#   ggplot() +
#   # 90 CI crossbars
#   # Need different lines for tests due to the different x-axis interval sizes
#   geom_crossbar(aes(x = index_vals, y = 0, ymin = q05, ymax = q95),
#                 fatten = 0, fill = "grey70", colour = NA, width = 1) +
#   # IQR Crossbars
#   geom_crossbar(aes(x = index_vals, y = 0, ymin = q25, ymax = q75),
#                 fatten = 0, fill = "grey50", width = 1) +
#   # Median segments
#   geom_crossbar(aes(x = index_vals, y = 0, ymin = q50, ymax = q50),
#                 fatten = 0, fill = NA, colour = "black", width = 1) +
#   geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
#   geom_segment(data = slope_length, aes(x = start_val, y = start_point, xend = end_val, yend = end_point)) +
#   scale_colour_brewer(palette = "Dark2") +
#   scale_x_continuous(expand = c(0, 0)) +
#   facet_grid(var ~ test, scales = "free", switch = "both") +
#   theme(legend.position = "top")

# Okay, that is looking much better now
# The issue was getting at how the linear model wanted to interact with the length test data
# Having now brought the statistics to Mohamad everything is swell

# The function for correct trend calculatoin was moved to the functions script

# Calculate the table for the Best Practices section
slope_final <- random_quant %>%
  gather(key = "stat", value = "val", -c(test:id)) %>%
  mutate(test2 = test,
         var2 = var) %>%
  group_by(test2, var2, id, stat) %>%
  group_modify(~trend_correct(.x)) %>%
  dplyr::rename(test = test2,
                var = var2) %>%
  # unite(var, id, col = "var_id", sep = " - ") %>%
  ungroup() %>%
  select(-id) %>%
  mutate(slope = round(slope, 2),
         R2 = paste0("(",R2,")")) %>%
  unite(slope, R2, col = "slope_R2", sep = " ") %>%
  select(-intercept, -p) %>%
  # gather(slope, R2, p, key = "var", value = "val") %>%
  spread(stat, slope_R2) %>%
  # select(test, var, range, q50, iqr50, iqr90) %>%
  select(test, var, range, q05, q25, q50, q75, q95) %>%
  mutate(test = factor(test, levels = c("length", "missing", "interp", "trend"))) %>%
  arrange(test)

# Smack it together to round this puppy out
best_table_average <- slope_final %>%
  filter(!grepl("focus", var))
saveRDS(best_table_average, "data/best_table_average.Rda")

best_table_focus <- slope_final %>%
  filter(grepl("focus", var))
saveRDS(best_table_focus, "data/best_table_focus.Rda")


# Discussion --------------------------------------------------------------

# Support the claim about interannual variability relating to changes in results with 30+ year base periods

# Support the claim about how many consecutive missing days is acceptable
# Most claims about missing data need to be backed up

# I don't think it is possible to be making the claims about max. intensity corrections that I am making

# Supplementary 1 ---------------------------------------------------------

# The effect of the sub-optimal tests on seas/thresh


# Supplementary 2 ---------------------------------------------------------

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


# Supplementary 3 ---------------------------------------------------------

# The global patterns in missing data are unremarkable and generally consistent across the oceans


# Supplementary 4 ---------------------------------------------------------

# The global patterns in added decadal trends generally show that MHW metrics increase


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
