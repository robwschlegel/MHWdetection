# code/functions.R
# This script contains all of the functions used throughout 'code/workflow.R'


# Libraries ---------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))

# library(jsonlite, lib.loc = "~/R-packages/")
# library(dplyr, lib.loc = "~/R-packages/") # Development version for group_modify()
library(tidyverse, lib.loc = "~/R-packages/")
library(ggridges)
# library(broom)
library(heatwaveR, lib.loc = "~/R-packages/")
# cat(paste0("heatwaveR version = ",packageDescription("heatwaveR")$Version))
library(lubridate) # This is intentionally activated after data.table
# library(fasttime)
library(ggpubr)
# library(boot)
# library(FNN)
# library(mgcv)
library(doMC); doMC::registerDoMC(cores = 50)
library(tidync, lib.loc = "~/R-packages/")
library(rgdal, lib.loc = "~/R-packages/")
library(pgirmess)
# library(egg)
# library(rcompanion)


# Meta-data ---------------------------------------------------------------

# The MHW category colour palette
MHW_colours <- c(
  "I Moderate" = "#ffc866",
  "II Strong" = "#ff6900",
  "III Severe" = "#9e0000",
  "IV Extreme" = "#2d0000"
)

# Location of NOAA OISST files
OISST_files <- dir(path = "~/data/OISST", full.names = T)

# Loation of MHW category files
category_files <- as.character(dir(path = "~/data/cat_clim", pattern = "cat.clim",
                                   full.names = TRUE, recursive = TRUE))

# Reference site coordinates
sst_ALL_coords <- data.frame(site = c("WA", "NW_Atl", "Med"),
                             lon = c(112.5, -67, 9),
                             lat = c(-29.5, 43, 43.5))

# The base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, col = "grey80", plot = FALSE)) %>%
  dplyr::rename(lon = long) %>%
  # filter(lat >= 25.6) %>%
  mutate(group = ifelse(lon > 180, group+9999, group),
         lon = ifelse(lon > 180, lon-360, lon))


# De-trending -------------------------------------------------------------

detrend <- function(df){
  resids <- broom::augment(lm(temp ~ t, df))
  res <- df %>%
    mutate(temp = round((temp - resids$.fitted),2))
  return(res)
}


# Control test variables --------------------------------------------------

# Length of time series
control_length <- function(year_begin, df){
  res <- df %>%
    filter(year(t) >= year_begin) %>%
    mutate(index_vals = 2018-year_begin+1,
           test = "length",
           fix = "none") %>%
    dplyr::select(test, fix, index_vals, t, temp)
  return(res)
}

# Missing data
control_missing <- function(prop, df){
  # NB: Don't allow sampling of first and last value to ensure
  # all time series are the same length
  ts_length <- nrow(df)
  miss_index <- sample(seq(2, ts_length-1, 1), ts_length*prop, replace = F)
  res <- df %>%
    mutate(row_index = 1:n(),
           temp = replace(temp, which(row_index %in% miss_index), NA),
           test = "missing") %>%
    mutate(index_vals = prop,
           fix = "none") %>%
    dplyr::select(test, fix, index_vals, t, temp)
  return(res)
}

# Linear decadal trend
control_trend <- function(rate, df){
  daily_step <- rate/3652.5
  res <- df %>%
    mutate(row_index = 1:n(),
           temp = temp + (row_index*daily_step),
           test = "trend") %>%
    mutate(index_vals = rate,
           fix = "none") %>%
    dplyr::select(test, fix, index_vals, t, temp)
  return(res)
}


# Count consecutive missing days ------------------------------------------

con_miss <- function(df){
  ex1 <- rle(is.na(df$temp))
  ind1 <- rep(seq_along(ex1$lengths), ex1$lengths)
  s1 <- split(1:nrow(df), ind1)
  if(length(s1) == 1){
    res <- tibble(duration = 1, count = 0)
  } else{
    res <- do.call(rbind,
                   lapply(s1[ex1$values == TRUE], function(x)
                     data.frame(index_start = min(x), index_end = max(x))
                   )) %>%
      mutate(duration = index_end - index_start + 1) %>%
      group_by(duration) %>%
      summarise(count = n())
  }
  return(res)
}


# Calculate clims and metrics ---------------------------------------------

# testers...
# df <- filter(sst_length, index_vals == 30)
# df <- filter(sst_interp, index_vals == 0.2)
# df <- filter(sst_window_10, index_vals == 10)
# set_window = 5
# set_window = 30
# set_pad = F
# min_date = "2009-01-01"
# year_start = 0
# year_end = 0
# focus_dates = focus_event
clim_metric_calc <- function(df, set_window = 5, set_pad = F, min_date = "2009-01-01",
                             year_start = 0, year_end = 0, focus_dates){

  # First and last years for full clim period
  if(year_start == 0)  year_start <- min(lubridate::year(df$t))
  if(year_end == 0) year_end <- max(lubridate::year(df$t))

  # base calculation
  res <- ts2clm(df, windowHalfWidth = set_window, maxPadLength = set_pad, var = TRUE,
                climatologyPeriod = c(paste0(year_start,"-01-01"), paste0(year_end,"-12-31"))) %>%
    filter(t >= min_date) %>%
    detect_event()

  # The metrics for the event(s) occurring during the largest control event
  res_focus <- res$event %>%
    filter(date_peak >= focus_dates$date_start,
           date_peak <= focus_dates$date_end) %>%
    summarise(count = n(),
              duration = sum(duration),
              intensity_max = max(intensity_max),
              intensity_cumulative = sum(intensity_cumulative)) %>%
    gather(var, val) %>%
    mutate(id = "focus_event",
           var = paste0("focus_",var),
           val = if_else(is.infinite(val), 0, val)) %>%
    select(id, var, val)

  # Extract desired clim values
  res_clim <- res$climatology %>%
    dplyr::select(doy, seas, thresh, var) %>%
    unique() %>%
    arrange(doy) %>%
    gather(var, val, -doy) %>%
    dplyr::rename(id = doy) %>%
    mutate(id = paste0("doy_",id)) %>%
    select(id, var, val)

  # Extract desired metric values
  res_metric <- res$event %>%
    dplyr::select(event_no, duration,
                  intensity_cumulative, intensity_max) %>%
    gather(var, val, -event_no) %>%
    dplyr::rename(id = event_no) %>%
    mutate(id = paste0("event_no_",id))%>%
    select(id, var, val)

  # Combine and exit
  res_all <- rbind(res_clim, res_metric, res_focus) %>%
    mutate(val = round(val, 3))
  return(res_all)
}


# Significance tests ------------------------------------------------------

# tester...
# df <- sst_clim_metric %>%
#   filter(test == "missing", var == "duration") %>%
#   ungroup() %>%
#   select(-test, -var)

# Kruskal-Wallis post-hoc
kruskal_post_hoc <- function(df){
  df$index_vals <- factor(df$index_vals)
  if("30" %in% levels(df$index_vals)){
    df$index_vals <- relevel(df$index_vals, "30")
  }
  res <- kruskalmc(df$val, df$index_vals)$dif.com %>%
    mutate(comp_index = row.names(.)) %>%
    separate(comp_index, into = c("control", "index_vals"),
             sep = '-', convert = T) %>%
    filter(control == levels(df$index_vals)[1]) %>%
    # filter(grepl(levels(df$index_vals)[1], control)) %>%
    select(-obs.dif, -critical.dif, -control) #%>%
    # mutate(difference = replace_na(difference, FALSE))
  return(res)
}

# Tukey post-hoc
tukey_post_hoc <- function(df){
  res <- TukeyHSD(aov(val ~ as.factor(index_vals), df)) %>%
    broom::tidy(.) %>%
    separate(comparison, into = c("comp", "control"), sep = '-') %>%
    mutate(p.value = round(adj.p.value, 4)) %>%
    filter(control == levels(df$index_vals)[1]) %>%
    arrange(adj.p.value)
  return(aov_tukey)
}


# Summary stats -----------------------------------------------------------

# Boot strap mean for confidence intervals
  # NB: CI's don't appear to be very useful/informative
  # This is because the different tests do not have a large effect on
  # the mean and distribution of the clims/metrics
# boot_mean <- function(data, indices) {
#   d <- data[indices] # allows boot to select sample
#   return(mean(d))
# }

# Summary stats for broad results
# tester...
# df <- sst_clim_metric %>%
# ungroup() %>%
# filter(test == "missing") %>%
# select(-test)
summary_stats <- function(df){

  # Determine the control group
  if(30 %in% df$index_vals){
    control_val <- 30
  } else{
    control_val <- 0
  }

  # Count of values
  res_count <- df %>%
    group_by(index_vals) %>%
    count(var) %>%
    filter(var == "duration") %>%
    select(-var) %>%
    unique() %>%
    data.frame()

  # Summary stats for each index_val
  res_base <- df %>%
    group_by(index_vals, var) %>%
    summarise(min = min(val, na.rm = T),
              # lower = boot.ci(boot(data = val, statistic = boot_mean,
              #                      R = 1000), type = "basic")$basic[4],
              median = median(val, na.rm = T),
              mean = mean(val, na.rm = T),
              # upper = boot.ci(boot(data = val, statistic = boot_mean,
              #                      R = 1000), type = "basic")$basic[5],
              max = max(val, na.rm = T),
              sum = sum(val, na.rm = T),
              # range = max-min,
              sd = sd(val, na.rm = T)) %>%
    ungroup() %>%
    left_join(res_count, by = "index_vals") %>%
    mutate(index_vals = as.character(index_vals)) %>%
    mutate_if(is.numeric, round, 3) %>%
    mutate(index_vals = as.numeric(index_vals))

  # Extract control row
  res_control <- filter(res_base, index_vals == control_val)

  # Find proportions and exit
  res_perc <- res_base %>%
    mutate(min_perc = (min-res_control$min)/abs(res_control$min),
           median_perc = (median-res_control$median)/abs(res_control$median),
           mean_perc = (mean-res_control$mean)/abs(res_control$mean),
           max_perc = (max-res_control$max)/abs(res_control$max),
           sum_perc = (sum-res_control$sum)/abs(res_control$sum),
           sd_perc = (sd-res_control$sd)/abs(res_control$sd),
           n_diff = n-res_control$n)
  return(res_perc)
}

# Summary stats for focus events only
# df <- sst_clim_metric %>%
# ungroup() %>%
# filter(id == "focus_event", test == "missing") %>%
# select(-test)
summary_stats_focus <- function(df){

  # Determine the control group
  if(30 %in% df$index_vals){
    control_val <- 30
  } else{
    control_val <- 0
  }

  # Extract control row
  res_control <- filter(df, index_vals == control_val)

  # Find proportions and exit
  res_perc <- df %>%
    select(-id) %>%
    mutate(val_perc = round((val-res_control$val)/abs(res_control$val), 2))
    # mutate(count = (focus_count-res_control$)/abs(res_control$min),
           # median_perc = (median-res_control$median)/abs(res_control$median),
           # mean_perc = (mean-res_control$mean)/abs(res_control$mean),
           # max_perc = (max-res_control$max)/abs(res_control$max),
           # sum_perc = (sum-res_control$sum)/abs(res_control$sum),
           # sd_perc = (sd-res_control$sd)/abs(res_control$sd),
           # n_diff = n-res_control$n)
  return(res_perc)
}


# Model linear relationships ----------------------------------------------

## Currently not calculating R2 for first global pass
# Potentially will calculate these on the global results afterwards

## Ideas:
# Fit linear models to everything and extract R2 values
# Look at relationship between change in seas/thresh and other metrics over time
# Look at this relationship for the percentage values, too
# Look at relationship between change in count of events and the percentage change of other summary stats

## Clever way of calculating R2 for any number of paired columns
# specify columns to regress
# y_col <- colnames(sst_summary)#"disp"
# y_col <- c("index_vals", "n", "n_diff", "mean")
# x_col <- colnames(sst_summary)[-c(1:3)]
#
# test <- expand.grid(y = y_col, x = x_col, stringsAsFactors = F) %>%
#   mutate(formula = paste(y,"~",x)) %>%
#   group_by(formula) %>%
#   mutate(r_sq = summary(lm(formula, data = filter(sst_summary, test == "missing", var == "seas")))$r.squared,
#          r_sq = round(r_sq, 2)) %>%
#   ungroup()


# Global functions --------------------------------------------------------

# The following wrapper functions use the above functions,
# but allow them to be used on the NOAA OISST NetCDF file structure

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

# This single function runs through and outputs all of the desired tests
single_analysis <- function(df){

  # Calculate the secular trend
  dec_trend <- round(as.numeric(broom::tidy(lm(temp ~ t, df))[2,2]*3652.5), 3)

  # Create base anomaly time series
    # NB: Rather just calcuate this for the globe by itself
    # It only takes a couple of hours to do this
  # sst_anom <- df %>%
  #   mutate(temp = round(temp-mean(temp), 2))
  # sst_anom_MHW <- detect_event(ts2clm(sst_anom, climatologyPeriod = c("1982-01-01", "2018-12-31")))

  # Create flat time series
  # sst_flat <- detrend(df)

  # Calculate MHWs from detrended
  sst_flat_MHW <- detect_event(ts2clm(detrend(df), climatologyPeriod = c("1982-01-01", "2018-12-31")))

  # Pull out the largest event in the ts
  focus_event <- sst_flat_MHW$event %>%
    filter(intensity_cumulative == max(intensity_cumulative)) %>%
    select(event_no, date_start:date_end, duration, intensity_cumulative, intensity_max)

  # Sub-optimise data
  ## Length
  sst_length <- plyr::ldply(1982:2009, control_length, df = sst_flat)
  ## Missing data
  sst_missing <- plyr::ldply(seq(0.00, 0.50, 0.01), control_missing, df = sst_flat)
  ## Decadal trend
  sst_trend <- plyr::ldply(seq(0.00, 0.30, 0.01), control_trend, df = sst_flat)

  # Calculate MHWs in most recent 10 years of data and return the desired clims and metrics
  system.time(
    sst_base_res <- rbind(sst_length, sst_missing, sst_trend) %>%
      group_by(test, index_vals) %>%
      group_modify(~clim_metric_calc(.x, focus_dates = focus_event)) %>%
      ungroup()
  ) # 22 seconds

  # Run the tests while also interpolating all gaps
  system.time(
    sst_interp_res <- sst_missing %>%
      group_by(test, index_vals) %>%
      group_modify(~clim_metric_calc(.x, focus_dates = focus_event, set_pad = 9999)) %>%
      ungroup() %>%
      mutate(test = "interp")
  ) # 11 seconds

  ### NB: It appears that increasing the window half width has a direct effect on reducing the MHW count
  ### This means that this test must almost certainly not be performed
  ### Need to first test this on several ts first before axing it

  # Increase rolling mean window to 10
  system.time(
    sst_window_10_res <- sst_length %>%
      group_by(test, index_vals) %>%
      group_modify(~clim_metric_calc(.x, focus_dates = focus_event, set_window = 10)) %>%
      ungroup() %>%
      mutate(test = "window_10")
  ) # 5 seconds

  # Increase rolling mean window to 20
  system.time(
    sst_window_20_res <- sst_length %>%
      group_by(test, index_vals) %>%
      group_modify(~clim_metric_calc(.x, focus_dates = focus_event, set_window = 20)) %>%
      ungroup() %>%
      mutate(test = "window_20")
  ) # 5 seconds

  # Increase rolling mean window to 30
  system.time(
    sst_window_30_res <- sst_length %>%
      group_by(test, index_vals) %>%
      group_modify(~clim_metric_calc(.x, focus_dates = focus_event, set_window = 30)) %>%
      ungroup() %>%
      mutate(test = "window_30")
  ) # 5 seconds

  # Combine for ease of use
  sst_clim_metric <- rbind(sst_base_res, sst_interp_res, sst_window_10_res, sst_window_20_res,  sst_window_30_res)

  # Run ANOVA/Tukey on MHW results for three different tests
  system.time(
    sst_signif <- sst_clim_metric %>%
      filter(id != "focus_event") %>%
      group_by(test, var) %>%
      group_modify(~kruskal_post_hoc(.x))
  ) # 1 second

  # Create summary statistics of MHW results
  system.time(
    sst_summary <- sst_clim_metric %>%
      filter(id != "focus_event") %>%
      group_by(test) %>%
      group_modify(~summary_stats(.x)) %>%
      left_join(sst_signif, by = c("test", "index_vals", "var")) %>%
      mutate(difference = replace_na(difference, FALSE))
  ) # 1 second

  # Effect on focus MHW
  sst_focus <- sst_clim_metric %>%
    filter(id == "focus_event") %>%
    group_by(test) %>%
    group_modify(~summary_stats_focus(.x))

  # Count consecutive missing days
    # NB: This would be useful information to have but it takes to long to calculate
  # system.time(
  #   sst_missing_count <- sst_missing %>%
  #     group_by(index_vals) %>%
  #     group_modify(~con_miss(.x))
  # ) # 18 seconds

  # Wrap it up
  res <- list(dec_trend = dec_trend,
              clim_metric = sst_clim_metric,
              summary = sst_summary,
              focus = sst_focus)
  return(res)
}

# The function that runs all of the tests on a single pixel/time series
# nc_file <- OISST_files[100]
global_analysis <- function(nc_file){

  sst <- tidync(nc_file) %>%
    hyper_tibble() %>%
    dplyr::rename(t = time, temp = sst) %>%
    mutate(t = as.Date(t, origin = "1970-01-01")) %>%
    filter(t <= "2018-12-31") %>%
    na.omit() %>%
    # Some ice pixels don't start on 1982-01-01, or go until 2018-12-31
    # So for now we are simply removing them
    group_by(lon, lat) %>%
    filter(max(t) == "2018-12-31",
           min(t) == "1982-01-01") %>%
    select(lon, lat, t, temp) %>%
    mutate(site = as.character(lat), rep = "1")

  # Also filter out pixels where there is 50% or more ice cover during the time series

  # Calculate the secular trend
  dec_trend <- round(as.numeric(broom::tidy(lm(temp ~ t, sst))[2,2]*3652.5), 3)

  # Remove random data
  # set.seed(666)
  sst_knockout <- plyr::ldply(.data = c(0.10, 0.25, 0.50),
                              .fun = random_knockout_global, df = sst)

  # Manually force certain key dates to be present
  sst_knockout$temp[sst_knockout$t == "1982-01-01"] <- sst$temp[sst$t == "1982-01-01"]
  sst_knockout$temp[sst_knockout$t == "1989-01-01"] <- sst$temp[sst$t == "1989-01-01"]
  sst_knockout$temp[sst_knockout$t == "1999-01-01"] <- sst$temp[sst$t == "1999-01-01"]
  sst_knockout$temp[sst_knockout$t == "2009-01-01"] <- sst$temp[sst$t == "2009-01-01"]
  sst_knockout$temp[sst_knockout$t == "2018-12-31"] <- sst$temp[sst$t == "2018-12-31"]

  # Make base calculations
  sst_res_base <- clim_event_cat_calc(sst) %>%
    mutate(index_vals = 0, test = as.factor("missing"))
  sst_res_base_fix <- sst_res_base %>%
    mutate(test = as.factor("missing_fix"))

  # Make missing data calculations
  sst_res_missing <- plyr::ddply(sst_knockout, c("index_vals"),
                                 clim_event_cat_calc) %>%
    mutate(test = as.factor("missing"))

  # Make missing data fix calculations
  sst_res_missing_fix <- plyr::ddply(sst_knockout, c("index_vals"),
                                     clim_event_cat_calc, fix = "missing") %>%
    mutate(test = as.factor("missing_fix"))

  # Make shortened time series calculations
  sst_res_length <- plyr::ldply(.data = seq(1989, 2009, by = 10),
                                .fun = shrinking_results_global, df = sst) %>%
    mutate(test = as.factor("length"))

  # Combine for upcoming tests
  sst_res <- rbind(sst_res_base, sst_res_missing,
                   sst_res_base_fix, sst_res_missing_fix,
                   sst_res_length) %>%
    select(test, index_vals, clim, event, cat) %>%
    mutate(site = as.character(lat), rep = "1")


  # Find the focus event
  focus_event <-  sst_res_base %>%
    select(-clim, -cat) %>%
    unnest(event) %>%
    filter(date_start >= "2009-01-01") %>%
    filter(intensity_cumulative == max(intensity_cumulative)) %>%
    mutate(site = as.character(lat)) %>%
    select(site, date_start:date_end) %>%
    dplyr::rename(date_start_control = date_start,
                  date_peak_control = date_peak,
                  date_end_control = date_end)

  # Quantify changes caused by the two tests
  effect_clim_res <- effect_clim_func(sst_res) %>%
    select(-site)
  effect_event_res <- effect_event_func(sst_res, focus_event) %>%
    select(-site)
  effect_cat_res <- effect_cat_func(sst_res, focus_event) %>%
    select(-site)

  # Wrap it up
  res <- list(lat = lat, lon = lon, dec_trend = dec_trend,
              effect_clim = effect_clim_res,
              effect_event = effect_event_res,
              effect_cat = effect_cat_res)
  return(res)
}

# One function to rule them all
# nc_file <- OISST_files[20]
global_analysis <- function(nc_file){

  # test_run <- global_analysis_sub()
  # system.time(
    suppressWarnings( # Ignore the coercing factor to character messages
      res_pixel <- plyr::llply(.data = seq(1, 720),
                               .fun = global_analysis_sub,
                               nc_file = nc_file, .parallel = T)
    )
  # )  # ~7 seconds for one, ~78 seconds for 720
  # save(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rdata"))
  return(res_pixel)
}

# The function that unpacks global resuls as desired
# testers...
# global_file <- global_files[1]
# sub_level <- 3
# sub_level <- 5
# sub_level <- 8
global_unpack_sub <- function(global_file, sub_level){
  # Load the one slice
  load(global_file)
  slice_full <- data.frame()
  for(i in 1:720){
    slice_step_1 <- slice_res[[i]]
    if(is.null(slice_step_1)){
     # slice_full <- slice_full
    } else if(nrow(as.data.frame(slice_step_1[sub_level])) > 0) {
      slice_step_2 <- as.data.frame(slice_step_1[c(1, 2, sub_level)])
      slice_full <- rbind(slice_full, slice_step_2)
    }
  }
  return(slice_full)
}

# The function that crawls through all of the global results
# unpacks them and then stitches them together before saving
global_unpack <- function(){
  # Point to files
  global_files <- dir("data/global", full.names = T)

  ## Run unpacker, save, and clear one at a time due to RAM restrictions on laptop
  # Decadal trends
  print("Unpacking decadal trends")
  global_dec_trend <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                                  .parallel = T, sub_level = 3)
  save(global_dec_trend, file = "data/global_dec_trend.Rdata")
  rm(global_dec_trend); gc()
  # KS climatology results
  print("Unpacking KS climatology results")
  global_KS_clim <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                               .parallel = T, sub_level = 4)
  colnames(global_KS_clim) <- gsub("KS_clim.", "", colnames(global_KS_clim))
  save(global_KS_clim, file = "data/global_KS_clim.Rdata")
  rm(global_KS_clim); gc()
  # KS event results
  print("Unpacking KS event results")
  global_KS_event <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                               .parallel = T, sub_level = 5)
  colnames(global_KS_event) <- gsub("KS_event.", "", colnames(global_KS_event))
  save(global_KS_event, file = "data/global_KS_event.Rdata")
  rm(global_KS_event); gc()
  # KS category results
  print("Unpacking KS category results")
  global_KS_cat <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                               .parallel = T, sub_level = 6)
  colnames(global_KS_cat) <- gsub("KS_cat.", "", colnames(global_KS_cat))
  save(global_KS_cat, file = "data/global_KS_cat.Rdata")
  rm(global_KS_cat); gc()
  # Single event climatology results
  print("Unpacking event climatology results")
  global_effect_clim <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                               .parallel = T, sub_level = 7)
  colnames(global_effect_clim) <- gsub("effect_clim.", "", colnames(global_effect_clim))
  save(global_effect_clim, file = "data/global_effect_clim.Rdata")
  rm(global_effect_clim); gc()
  # Single evnt metric results
  print("Unpacking event metric results")
  global_effect_event <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                               .parallel = T, sub_level = 8)
  colnames(global_effect_event) <- gsub("effect_event.", "", colnames(global_effect_event))
  save(global_effect_event, file = "data/global_effect_event.Rdata")
  rm(global_effect_event); gc()
  # Single event category results
  print("Unpacking event category results")
  global_effect_cat <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                               .parallel = T, sub_level = 9)
  colnames(global_effect_cat) <- gsub("effect_cat.", "", colnames(global_effect_cat))
  save(global_effect_cat, file = "data/global_effect_cat.Rdata")
  rm(global_effect_cat); gc()
}

# Function for calculating R2 values for sub-optimal results
global_slope_sub <- function(df){
  if(nrow(df) >2){
    round(broom::tidy(lm(val ~ index_vals, data = df))$estimate[2], 3)
  } else{
    return(NA)
  }
}

# This function expects to be given only one latitude slice at a time
# tester...
# df <- global_effect_event %>%
  # filter(lon == lon[20], lat == lat[129])
  # filter(lon == lon[20], lat == lat[129], test == "length", metric == "duration")
global_slope <- function(df){
  suppressWarnings( # Suppress perfect slope warnings
  df_slope <- df %>%
    group_by(lon, test, metric) %>%
    nest() %>%
    mutate(slope = purrr::map(data, global_slope_sub)) %>%
    select(-data) %>%
    unnest()
  )
}

# This function is the same as above but less picky about the output
# tester...
# df <- event_slope_dec_trend %>%
  # filter(lon == lon[20], lat == lat[129])
  # filter(lon == lon[20], lat == lat[129], test == "length", metric == "duration")
global_model <- function(df){
  suppressWarnings( # Suppress perfect slope warnings
    df_slope <- df %>%
      group_by(lon, test, metric) %>%
      # nest() %>%
      # mutate(slope = purrr::map(data, global_slope_sub)) %>%
      do(model = broom::augment(lm(slope ~ dec_trend, data = .))) #%>%
      # select(-data) %>%
      # unnest()
  )
}

# Length fix --------------------------------------------------------------

# Below are the functions used to perform post-hoc corrections on events
# detected in shorter time series

# To correct for the effect of time series length on event metrics we need
# two things:
  # The length of the time series
  # The decadal trend in the region
# It may be that rather than correcting the value, we should provide a CI instead

# The relationship may be better aided by the known variance in the time series

# Must see what the R2 + SE is between decadal trend and change in duration/max.int.


# Figure functions --------------------------------------------------------

# Wrapper for creating plugs for smoother fig 2 - 4 plotting
control_plug <- function(test_plug, site_plug){
  plug <- data.frame(test = test_plug, site = site_plug, metric = unique(sst_ALL_plot_long$metric),
                     index_vals = 30, p.value.min = 1.0, p.value.mean = 1.0,
                     p.value.max = 1.0, p.value.sd = 0)
}


# The code that creates the figure 1 panels from detect_event output
# The function expects to be given the dates that should be plotted
# testers...
# df <- sst_ALL_clim
# y_label <-  "Temperature (°C)"
# spread <- 183
fig_1_plot <- function(df, spread, y_label = "Temperature (°C)"){
  # Create category breaks and select slice of data.frame
  clim_cat <- df %>%
    group_by(site_label) %>%
    dplyr::mutate(diff = thresh - seas,
                  thresh_2x = thresh + diff,
                  thresh_3x = thresh_2x + diff,
                  thresh_4x = thresh_3x + diff) %>%
    dplyr::filter(t >= date_peak-spread, t <= date_peak+spread)

  # Peak event
  peak_event <- clim_cat %>%
    group_by(site_label) %>%
    filter(event_no == clim_cat$event_no[clim_cat$t == date_peak]) %>%
    mutate(date_start = min(t),
           date_end = max(t),
           start_point = thresh[t == min(t)],
           peak_point = thresh[t == date_peak],
           end_point = thresh[t == max(t)])

  # Set line colours
  lineColCat <- c(
    "Temperature" = "black",
    "Climatology" = "skyblue",
    "Threshold" = "navy",
    "2x Threshold" = "gray20",
    "3x Threshold" = "gray30",
    "4x Threshold" = "gray40"
  )

  # Set category fill colours
  fillColCat <- c(
    "I Moderate" = "#ffc866",
    "II Strong" = "#ff6900",
    "III Severe" = "#9e0000",
    "IV Extreme" = "#2d0000"
  )

  ggplot(data = clim_cat, aes(x = t, y = temp)) +
    geom_flame(aes(y2 = thresh, fill = "I Moderate"), n = 5, n_gap = 2) +
    geom_flame(aes(y2 = thresh_2x, fill = "II Strong")) +
    geom_flame(aes(y2 = thresh_3x, fill = "III Severe")) +
    geom_flame(aes(y2 = thresh_4x, fill = "IV Extreme")) +
    geom_line(aes(y = thresh_2x, col = "2x Threshold"), size = 0.7, linetype = "dashed") +
    geom_line(aes(y = thresh_3x, col = "3x Threshold"), size = 0.7, linetype = "dotdash") +
    geom_line(aes(y = thresh_4x, col = "4x Threshold"), size = 0.7, linetype = "dotted") +
    geom_line(aes(y = seas, col = "Climatology"), size = 0.7) +
    geom_line(aes(y = thresh, col = "Threshold"), size = 0.7) +
    geom_line(aes(y = temp, col = "Temperature"), size = 0.6) +
    # geom_segment(data = peak_event, arrow = arrow(),
    #              aes(x = date_peak[1], xend = date_peak[1],
    #                  y = max(peak_event$temp) + 2,
    #                  yend = max(peak_event$temp))) +
    geom_segment(data = peak_event, colour = "springgreen",
                 aes(x = date_start-1, xend = date_start-1,
                     y = start_point-1,
                     yend = start_point+1)) +
    geom_segment(data = peak_event, colour = "springgreen",
                 aes(x = date_end+1, xend = date_end+1,
                     y = end_point-1,
                     yend = end_point+1)) +
    geom_segment(data = peak_event, colour = "forestgreen",
                 aes(x = date_peak, xend = date_peak,
                     y = peak_point-1,
                     yend = peak_point+1)) +
    # geom_rug(data = peak_event, sides = "b", colour = "red3", size = 2,
    #          aes(x = c(min(peak_event$t), max(peak_event$t)), y = min(clim_cat$temp))) +
    scale_colour_manual(name = NULL, values = lineColCat,
                        breaks = c("Climatology", "Threshold", "Temperature",
                                   "2x Threshold", "3x Threshold", "4x Threshold")) +
    scale_fill_manual(name = NULL, values = fillColCat,
                      breaks = c("I Moderate", "II Strong",
                                 "III Severe", "IV Extreme")) +
    scale_x_date(date_labels = "%b %Y", expand = c(0, 0)) +
    # scale_x_date(date_labels = "%b %Y", expand = c(0, 0),
    #              breaks = c(clim_cat$t[round(nrow(clim_cat)*0.25)],
    #                         clim_cat$t[round(nrow(clim_cat)*0.5)],
    #                         clim_cat$t[round(nrow(clim_cat)*0.75)])) +
    guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid",
                                                                  "dashed", "dotdash", "dotted")))) +
    labs(y = y_label, x = NULL) +
    facet_wrap(~site_label, ncol = 1, scales = "free") +
    theme(legend.position = "top")
}

# The code that creates figure 2 - 4
# These show the p-values for the KS tests from different sub-optimal treatments
# testers...
# df <- sst_ALL_plot_long
# site_sub <- "WA"
# test_sub <- "length"
fig_line_plot <- function(test_sub, df = sst_ALL_plot_long){

  # Create x -axis label
  if(test_sub == "length"){
    x_label <- "Time series length (years)"
  } else if(test_sub == "missing" | test_sub == "missing_fix"){
    x_label <- "Missing data (%)"
  } else if(test_sub == "trended"){
    x_label <- "Added linear trend (°C/dec)"
  } else{
    stop("Typos make baby pandas cry...")
  }

  # Filter data
  df_sub <- df %>%
    filter(test == test_sub)

  # Correct percentage missing index_vals
  if(test_sub == "missing" | test_sub == "missing_fix"){
    df_sub$index_vals <- df_sub$index_vals*100
  }

  # Set label names and colours
  colour_choice <- c("skyblue", "navy",
                     "springgreen", "forestgreen",
                     "#ffc866", "#ff6900", "#9e0000", "#2d0000")
  label_choice <- c("Climatology", "Threshold",
                   "Duration (days)", "Max. intensity (°C)",
                   "Prop. moderate", "Prop. strong", "Prop. severe", "Prop. extreme")
  # Create figure
  fig_plot <- ggplot(df_sub, aes(x = index_vals, y = p.value.mean, colour = metric)) +
    geom_line() +
    # geom_point() +
    geom_ribbon(alpha = 0.1, colour = NA,
                aes(ymin = p.value.mean-p.value.sd,
                    ymax = p.value.mean+p.value.sd, fill = metric)) +
    geom_point(data = filter(df_sub, p.value.mean <= 0.05), shape = 15, colour = "red", size = 2) +
    geom_hline(yintercept = 0.05, colour = "red", linetype = "dashed") +
    scale_colour_manual(values = colour_choice,
                        labels = label_choice) +
    scale_fill_manual(values = colour_choice,
                      labels = label_choice) +
    # scale_y_continuous(limits = c(0, 1)) +
    coord_cartesian(ylim = c(0, 1), expand = F) +
    # geom_errorbarh(aes(xmin = year_long, xmax = year_short)) +
    # guides(colour = guide_legend(override.aes = list(shape = 15, linetype = NA, size = 3))) +
    labs(y = "Mean p-value +- SD", x = x_label, colour = NULL, fill = NULL) +
    # facet_grid(site~test, scales = "free_x", switch = "x") +
    theme(legend.position = "top") +
    facet_wrap(~site_label)
  # fig_plot
  return(fig_plot)
}


# Function that plots the effect on a single focus MHW



# Function for easily plotting subsets from the global slope results
# testers...
# test_sub <- "length"
# metric_sub <- "intensity_max"
# metric_sub <- "duration"
global_effect_event_slope_plot <- function(test_sub, metric_sub,
                                           df = global_effect_event_slope,
                                           prop = FALSE) {

  # Prepare Viridis colour palette
  if(metric_sub == "duration") {
    vir_op <- "C"
    col_split <- c("purple", "forestgreen")
  } else {
    vir_op <- "A"
    col_split <- c("blue", "red")
  }

  # Filter base data
  base_sub <- df %>%
    filter(test == test_sub, metric == metric_sub) #%>%
    # Convert from proportion to percent
    # mutate(slope = slope*100)

  if(prop){
    type_sub <- "prop"
    base_sub$slope <- base_sub$slope*100
  } else{
    type_sub <- "slope"
  }

  # Find quantiles
  slope_quantiles <- quantile(base_sub$slope, na.rm = T,
                              probs = c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1.0))

  # Correct base data to quantiles as the tails are very long
  base_quantile <- base_sub %>%
    mutate(slope = case_when(slope > slope_quantiles[6] ~ slope_quantiles[6],
                             slope < slope_quantiles[2] ~ slope_quantiles[2],
                             slope <= slope_quantiles[6] | slope >= slope_quantiles[2] ~ slope))

  # The map
  slope_map <- ggplot(base_quantile, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
    # scale_fill_viridis_c(option = vir_op) +
    # scale_fill_gradientn(colors = scales::viridis_pal(option = vir_op)(9),
    # limits = c(as.numeric(slope_quantiles[2]),
    # as.numeric(slope_quantiles[4])),
    # breaks = c(as.numeric(slope_quantiles[2:6]))) +
    scale_fill_gradient2(low = col_split[1], high = col_split[2],
                         breaks = c(as.numeric(slope_quantiles[2:6]))) +
    coord_equal(expand = F) +
    # labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "cm"))
  if(metric_sub == "intensity_max" & prop == F){
    slope_map <- slope_map + labs(fill = "Change in max. intensity (°C)\nper year")
  } else if(metric_sub == "duration" & prop == F){
    slope_map <- slope_map + labs(fill = "Change in duration (days)\nper year")
  } else if(metric_sub == "intensity_max" & prop == T){
    slope_map <- slope_map + labs(fill = "Percent change in max. intensity (°C)\n per year from 10 year value")
  } else if(metric_sub == "duration" & prop == T){
    slope_map <- slope_map + labs(fill = "Percent change in duration (days)\n per year from 10 year value")
  }
  # slope_map

  # The density polygon
  # slope_density <- ggplot(base_sub, aes(x = slope)) +
  #   geom_density(aes(fill = slope)) +
  #   coord_flip() +
  #   scale_x_continuous(expand = c(0,0))
  # slope_density

  # The ridgwlinw plot
  # slope_ridge <- ggplot(base_sub, aes(x = slope, y = metric)) +
  #   stat_density_ridges(aes(fill = factor(..quantile..)),
  #                       geom = "density_ridges_gradient", calc_ecdf = TRUE,
  #                       quantiles = 4, quantile_lines = TRUE) +
  #   viridis::scale_fill_viridis(discrete = TRUE, name = "Quartiles", alpha = 0.7) +
  #   coord_flip(expand = F) +
  #   theme(axis.text.x = element_blank())
  # slope_ridge

  ggsave(slope_map,
         filename = paste0("output/",test_sub,"_",metric_sub,"_",type_sub,"_plot.png"), height = 6, width = 10)
  return(slope_map)
}


# Old figure functions ----------------------------------------------------

# Expects a one row data.frame with a 'lon' and 'lat' column
# df <- focus_WA
map_point <- function(df){
  category_data <- readRDS(category_files[grepl(pattern = as.character(df$date_peak),
                                                x = category_files)])
  map_out <- ggplot(data = df, aes(x = lon, y = lat)) +
    geom_tile(data = category_data, aes(fill = category)) +
    borders(fill = "grey80", colour = "black") +
    geom_point(shape = 21, colour = "white", fill = "hotpink", size = 2) +
    scale_fill_manual("Category",
                      values = c("#ffc866", "#ff6900", "#9e0000", "#2d0000"),
                      labels = c("I Moderate", "II Strong",
                                 "III Severe", "IV Extreme")) +
    coord_cartesian(xlim = c(df$lon[1]-20, df$lon[1]+20),
                    ylim = c(df$lat[1]-20, df$lat[1]+20)) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "bottom")
  return(map_out)
}

# This function expects the output of the clim, event, cat pipe
# df <- filter(sst_ALL_res, site == "NW_Atl")
table_summary <- function(site_1){

  df <- filter(sst_ALL_res, site == site_1)

  # Seasonal min/mean/max
  # Threshold min/mean/max
  clim_long <- df %>%
    select(clims) %>%
    unnest() %>%
    select(doy, seas, thresh) %>%
    unique() %>%
    select(seas, thresh) %>%
    gather(key = "metric", value = "val")

  event_long <- df %>%
    select(events) %>%
    unnest() %>%
    filter(row_number() %% 2 == 0) %>%
    unnest() %>%
    select(duration, intensity_mean, intensity_max, intensity_cumulative) %>%
    gather(key = "metric", value = "val")

  event_count <- df %>%
    select(cats) %>%
    unnest() %>%
    select(category) %>%
    group_by(category) %>%
    summarise(val = n()) %>%
    mutate(metric = "count") %>%
    spread(key = category, value = val)
  if(!"IV Extreme" %in% colnames(event_count)){
    event_count <- cbind(event_count, tibble('IV Extreme' = 0))
  }
  event_count <- event_count %>%
    select(metric, 'I Moderate', 'II Strong', 'III Severe', 'IV Extreme') %>%
    dplyr::rename(' ' = metric)

  summary_res <- rbind(clim_long, event_long) %>%
    group_by(metric) %>%
    summarise_all(.funs = c("min", "mean", "max", "sd")) %>%
    mutate_if(is.numeric, round, 2) %>%
    arrange(metric)

  tbl_1 <- gridExtra::tableGrob(summary_res, rows = NULL)
  tbl_2 <- gridExtra::tableGrob(event_count, rows = NULL)

  tbl_all <- gridExtra::grid.arrange(tbl_1, tbl_2,
                                     nrow = 2, as.table = TRUE)

  return(tbl_all)
  # res_list <- list(summary_res = summary_res,
                   # event_count = event_count)
  # return(res_list)
}

# Wrapper for time series + clims + event rug
ts_clim_rug <- function(site_1){
  ts_data <- filter(sst_ALL, site == site_1)
  clim_data <- filter(sst_ALL_clim, site == site_1)
  event_data <- filter(sst_ALL_event, site == site_1)
  fig <- ggplot(data = ts_data, aes(x = t, y = temp)) +
    geom_line(colour = "grey20") +
    geom_line(data = clim_data, aes(y = seas),
              linetype = "dashed", colour = "steelblue3") +
    geom_line(data = clim_data, linetype = "dotted", colour = "tomato3",
              aes(x = t, y = thresh)) +
    geom_rug(data = event_data, sides = "b", colour = "red3", size = 2,
             aes(x = date_peak, y = min(ts_data$temp))) +
    labs(x = NULL, y = "Temperature (°C)")
  return(fig)
}

# Wrapper for clim only line plot
clim_line <-function(site_1){
  clim_data <- filter(sst_ALL_clim_only, site == site_1)
  ggplot(data = clim_data, aes(x = doy)) +
    geom_line(aes(y = seas), colour = "steelblue3") +
    geom_line(aes(y = thresh), colour = "tomato3") +
    labs(x = NULL, y = "Temperature (°C)")
}

