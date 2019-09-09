# workflow.R


# Libraries ---------------------------------------------------------------

# Load all meta-data and functions
source("code/functions.R")

# Remove scientific notation from results
# options(scipen=999)


# Analyses ----------------------------------------------------------------

# Where on the x axis things go wrong is the main question to be answers

# Create a map at which the year after which a threshold is exceeded in the change in the statistic in question

# Include the maps showing the results with the trends left in as supplementary material

# Write a paragraph in the discussion that talks about why using p-values for deciding what not to use is a bad idea/problematic

# Show the difference in the moving 30 year clim vs. the preferred 30 year clim as an appendix figure

# It may end up being best to offer advise based on the change in MHW count/days
# Also the proportion shift in duration and max int based on something bio relevant in literature

# Need to run the above code on the three base ts and see what those results look like
# If that looks good then it is time to go global

# Combine the three reference ts, run results, and save
# sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
#                  mutate(sst_NW_Atl, site = "NW_Atl"),
#                  mutate(sst_Med, site = "Med"))
# system.time(
  # sst_ALL_results <- plyr::dlply(sst_ALL, c("site"), single_analysis,
                                 # full_seq = T, clim_metric = T, count_miss = T, .parallel = T)
# ) # 54 seconds
# saveRDS(sst_ALL_results, "data/sst_ALL_results.Rds")

## Global run
# Run sequentially so that each lon slice can be saved en route
# i <- 6
# for(i in 1:length(OISST_files)){
# # for(i in 273){
#
#   # Determine file
#   OISST_slice <- OISST_files[i]
#   lon_row_pad <- str_pad(i, width = 4, pad = "0", side = "left")
#   print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))
#
#   # Calculate tests etc.
#   # system.time(
#   slice_res <- global_analysis(OISST_slice)
#   # ) # ~ 3 - 4 minutes
#   saveRDS(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rds"))
#   print(paste0("Finished run on step ",lon_row_pad," at ",Sys.time()))
#
#   # Clear up some RAM
#   rm(slice_res); gc()
# } # ~ 2.5 minutes each

# It appears as though we are experiencing some sort of core slippage when running
# one file on multiple cores
# So now we are going to try running multiple files on one core each
# It may also be that the MHW Trackers scheduled run is knocking this script out of order

global_analysis_single <- function(file_sub){
  OISST_slice <- OISST_files[file_sub]
  lon_row_pad <- str_pad(file_sub, width = 4, pad = "0", side = "left")
  print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))
  slice_res <- global_analysis(OISST_slice)
  saveRDS(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rds"))
  print(paste0("Finished run on step ",lon_row_pad," at ",Sys.time()))
  rm(slice_res); gc()
}

plyr::l_ply(1:1440, global_analysis_single, .parallel = T)

# I project that this may take 3 days...

# This took ~xxx hours to run


# Test visuals ------------------------------------------------------------

# Overall mean values
# ggplot(sst_test$Med$summary, aes(x = index_vals, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line(aes(colour = var), size = 2, alpha = 0.7) +
#   geom_point(data = filter(sst_summary, difference == TRUE),
#              colour = "red", size = 1, alpha = 0.4) +
#   # scale_colour_manual(values = c("black", "red")) +
#   facet_grid(~test, scales = "free")# +
#   # coord_cartesian(ylim = c(-2, 2))

# Percent change away from control valuesfor base three tests
# ggplot(filter(sst_test$Med$summary, test %in% c("length", "missing", "trend")),
#               aes(x = index_vals, y = mean_perc)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line(aes(colour = var), size = 2, alpha = 0.7) +
#   # geom_point(data = filter(sst_summary, difference == TRUE),
#              # colour = "red", size = 1, alpha = 0.4) +
#   # scale_colour_manual(values = c("black", "red")) +
#   facet_wrap(~test, scales = "free") +
#   coord_cartesian(ylim = c(-2, 2))

# Percent change away from control values with interpolation
# ggplot(filter(sst_test$Med$summary, test %in% c("missing", "interp")),
#        aes(x = index_vals, y = mean_perc)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line(aes(colour = var), size = 2, alpha = 0.7) +
#   # geom_point(data = filter(sst_summary, difference == TRUE),
#              # colour = "red", size = 1, alpha = 0.4) +
#   # scale_colour_manual(values = c("black", "red")) +
#   facet_wrap(~test, scales = "free") +
#   coord_cartesian(ylim = c(-1, 1))

# Change in count of MHWs
# ggplot(filter(sst_test$Med$summary, test %in% c("missing", "interp")),
#        aes(x = index_vals, y = n_diff)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line() +
#   facet_grid(~test, scales = "free")

# Percent change away from control values with interpolation
# ggplot(filter(sst_test$Med$summary, test %in% c("length", "window_10", "window_20", "window_30")),
#        aes(x = index_vals, y = mean_perc)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line(aes(colour = var), size = 2, alpha = 0.7) +
#   # geom_point(data = filter(sst_summary, difference == TRUE),
#   #            colour = "red", size = 1, alpha = 0.4) +
#   # scale_colour_manual(values = c("black", "red")) +
#   facet_wrap(~test, scales = "free") +
#   coord_cartesian(ylim = c(-1, 1))

# Change in count of MHWs
# ggplot(filter(sst_test$WA$summary, test %in% c("length", "window_10", "window_20", "window_30")),
#        aes(x = index_vals, y = n)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line() +
#   facet_grid(~test, scales = "free")

# Change in percent of MHW days
# ggplot(filter(sst_test$Med$summary, var == "duration"),
#        aes(x = index_vals, y = sum_perc)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line() +
#   facet_grid(~test, scales = "free")

# Change in focus event
# ggplot(filter(sst_focus, var == "focus_duration"),
# ggplot(sst_test$Med$focus,
#        aes(x = index_vals, y = val_perc)) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_line() +
#   facet_grid(var~test, scales = "free")


# Visuals -----------------------------------------------------------------

# Prep event data for pretty plotting
# effect_event_pretty <- effect_event %>%
#   filter(metric %in% c("count", "duration", "intensity_max"),
#          !index_vals %in% seq(1, 9),
#          index_vals <= 0.5 | index_vals >= 10) %>%
#   mutate(panel_label = case_when(metric == "count" & test == "length" ~ "A",
#                                  metric == "count" & test == "missing" ~ "B",
#                                  metric == "count" & test == "trended" ~ "C",
#                                  metric == "duration" & test == "length" ~ "D",
#                                  metric == "duration" & test == "missing" ~ "E",
#                                  metric == "duration" & test == "trended" ~ "F",
#                                  metric == "intensity_max" & test == "length" ~ "H",
#                                  metric == "intensity_max" & test == "missing" ~ "I",
#                                  metric == "intensity_max" & test == "trended" ~ "J"),
#          metric = case_when(metric == "intensity_max" ~ "max. intensity (°C)",
#                             metric == "duration" ~ "duration (days)",
#                             metric == "count" ~ "count (event)"),
#          test = case_when(test == "length" ~ "length (years)",
#                           test == "missing" ~ "missing data (proportion)" ,
#                           test == "trended" ~ "added trend (°C/dec)"),
#          test = as.factor(test),
#          test = factor(test, levels = levels(test)[c(2,3,1)]),
#          site = as.character(site)) %>%
#   group_by(test) %>%
#   mutate(panel_label_x = min(index_vals)) %>%
#   group_by(metric) %>%
#   mutate(panel_label_y = max(val)) %>%
#   ungroup()
# effect_event_pretty$site[effect_event_pretty$site == "NW_Atl"] <- "NWA"
# effect_event_pretty$site <- factor(effect_event_pretty$site, levels = c("WA", "NWA", "Med"))
# effect_event_pretty$index_vals[effect_event_pretty$test == "missing"] <- effect_event_pretty$index_vals[effect_event_pretty$test == "missing"]*100

### Visualise
## Climatologies
# Sub-optimal data
# ggplot(effect_clim, aes(x = index_vals)) +
#   # geom_ribbon(aes(ymin = min, ymax = max, fill = site), alpha = 0.2) +
#   geom_line(aes(y = mean, colour = site)) +
#   # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
#   facet_grid(metric~test, scales = "free")
# # Fixed data
# ggplot(effect_clim_fix, aes(x = index_vals)) +
#   # geom_ribbon(aes(ymin = min, ymax = max, fill = site), alpha = 0.2) +
#   geom_line(aes(y = mean, colour = site)) +
#   # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
#   facet_grid(metric~test, scales = "free")

## Event metrics
# Sub-optimal data
# plot_event_effect <- ggplot(effect_event_pretty, aes(x = index_vals)) +
#   # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
#   # geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
#   # stat_smooth(aes(y = val, colour = site), geom = "line",
#               # method = "lm", alpha = 0.5, size = 1) +
#   geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1) +
#   geom_text(aes(label = panel_label, y = panel_label_y, x = panel_label_x)) +
#   # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
#   scale_colour_brewer(palette = "Dark2") +
#   facet_grid(metric~test, scales = "free", switch = "both") +
#   labs(x = NULL, y = NULL, colour = "Site") +
#   theme(legend.position = "bottom")
# plot_event_effect
# ggsave(plot_event_effect, filename = "output/effect_event.pdf", height = 5, width = 10)
# ggsave(plot_event_effect, filename = "output/effect_event.png", height = 5, width = 10)
# ggsave(plot_event_effect, filename = "LaTeX/fig_3.pdf", height = 5, width = 10)
# ggsave(plot_event_effect, filename = "LaTeX/fig_3.png", height = 5, width = 10)
# ggsave(plot_event_effect, filename = "LaTeX/fig_3.jpg", height = 5, width = 10)

# Fixed data
# ggplot(effect_event_fix, aes(x = index_vals)) +
#   # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
#   geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
#   stat_smooth(aes(y = val, colour = site), geom = "line",
#               method = "lm", alpha = 0.5, size = 1) +
#   geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1.2) +
#   # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
#   facet_grid(metric~test, scales = "free", switch = "both") +
#   labs(x = NULL, y = NULL, colour = "Site") +
#   theme(legend.position = "bottom")

# Missing data only + fix


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
