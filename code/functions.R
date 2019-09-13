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
library(tidync, lib.loc = "~/R-packages/")
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

# Location of global result files
global_files <- dir(path = "data/global", full.names = T)

# Loation of MHW category files
# category_files <- as.character(dir(path = "~/data/cat_clim", pattern = "cat.clim",
#                                    full.names = TRUE, recursive = TRUE))

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
           test = "length") %>%
    dplyr::select(test, index_vals, t, temp)
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
           test = "missing",
           index_vals = prop) %>%
    dplyr::select(test, index_vals, t, temp)
  return(res)
}

# Linear decadal trend
control_trend <- function(rate, df){
  daily_step <- rate/3652.5
  res <- df %>%
    mutate(row_index = 1:n(),
           temp = temp + (row_index*daily_step),
           test = "trend",
           index_vals = rate) %>%
    dplyr::select(test, index_vals, t, temp)
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
clim_metric_focus_calc <- function(df, set_window = 5, set_pad = F, min_date = "2009-01-01",
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
    select(var, id, val)

  # Extract desired clim values
  res_clim <- res$climatology %>%
    dplyr::select(doy, seas, thresh, var) %>%
    unique() %>%
    arrange(doy) %>%
    gather(var, val, -doy) %>%
    dplyr::rename(id = doy) %>%
    mutate(id = paste0("doy_",id)) %>%
    select(var, id, val)

  # Extract desired metric values
  res_metric <- res$event %>%
    dplyr::select(event_no, duration,
                  intensity_cumulative, intensity_max) %>%
    gather(var, val, -event_no) %>%
    dplyr::rename(id = event_no) %>%
    mutate(id = paste0("event_no_",id))%>%
    select(var, id, val)

  # Combine and exit
  res_all <- rbind(res_clim, res_metric, res_focus) %>%
    mutate(val = round(val, 3))
  return(res_all)
}


# Window manipulation -----------------------------------------------------

window_test <- function(df, focus_event){

  # Increase rolling mean window to 10
  sst_window_10_res <- df %>%
    group_by(test, index_vals) %>%
    group_modify(~clim_metric_focus_calc(.x, focus_dates = focus_event, set_window = 10)) %>%
    ungroup() %>%
    mutate(test = "window_10")

  # Increase rolling mean window to 20
  sst_window_20_res <- df %>%
    group_by(test, index_vals) %>%
    group_modify(~clim_metric_focus_calc(.x, focus_dates = focus_event, set_window = 20)) %>%
    ungroup() %>%
    mutate(test = "window_20")

  # Increase rolling mean window to 30
  sst_window_30_res <- df %>%
    group_by(test, index_vals) %>%
    group_modify(~clim_metric_focus_calc(.x, focus_dates = focus_event, set_window = 30)) %>%
    ungroup() %>%
    mutate(test = "window_30")

  # Exit
  res <- rbind(sst_window_10_res, sst_window_20_res, sst_window_30_res)
  return(res)
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
    dplyr::rename(val = difference) %>%
    mutate(id = "difference",
           val = as.numeric(val)) %>%
    select(index_vals, id, val)
  return(res)
}

# Tukey post-hoc
  # NB: Not used as the data are not normal
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

# Summary stats for broad results
# tester...
# df <- sst_clim_metric %>%
# ungroup() %>%
# filter(test == "length") %>%
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
    mutate(var = "count",
           n_diff = n-(nrow(filter(df,
                                   index_vals == control_val,
                                   var == "duration")))) %>%
    gather(id, val, -index_vals, -var)

  # Summary stats for each index_val
  res_base <- df %>%
    group_by(index_vals, var) %>%
    summarise(min = min(val, na.rm = T),
              median = median(val, na.rm = T),
              mean = mean(val, na.rm = T),
              max = max(val, na.rm = T),
              sum = sum(val, na.rm = T),
              sd = sd(val, na.rm = T),
              sd = replace_na(sd, 0)) %>%
    gather(id, val, -index_vals, -var)

  # Extract control row
  res_control <- filter(res_base, index_vals == control_val) %>%
    ungroup() %>%
    dplyr::rename(cont_val = val) %>%
    select(-index_vals)

  # Find percentages of change, attach counts, and exit
  res_perc <- left_join(res_base, res_control, by = c("var", "id")) %>%
    mutate(perc = (val-cont_val)/abs(cont_val),
           perc = replace_na(perc, 0),
           perc_name = paste0(id,"_perc")) %>%
    select(index_vals, var, perc_name, perc) %>%
    dplyr::rename(id = perc_name, val = perc)

  # Combine, round, and exit
  res <- rbind(res_base, res_perc, res_count) %>%
    mutate(val = round(val, 3))
  return(res)
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


# Single full analysis ----------------------------------------------------

# This single function runs through and outputs all of the desired tests as a list
# testers...
# df <- sst_Med
# full_seq = T
# clim_metric = T
# count_miss = T
# windows = T
single_analysis <- function(df, full_seq = F, clim_metric = F, count_miss = F, windows = F){

  # Calculate the secular trend
  dec_trend <- data.frame(test = "length", index_vals = 37, var = "dec_trend", id = "slope",
                          val = round(as.numeric(broom::tidy(lm(temp ~ t, df))[2,2]*3652.5), 3),
                          stringsAsFactors = F)

  # Detrend the ts
  sst_flat <- detrend(df)

  # Calculate MHWs from detrended ts
  sst_flat_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2018-12-31")))

  # Pull out the largest event in the ts
  focus_event <- sst_flat_MHW$event %>%
    filter(date_start >= "2009-01-01") %>%
    filter(intensity_cumulative == max(intensity_cumulative)) %>%
    select(event_no, date_start:date_end, duration, intensity_cumulative, intensity_max) %>%
    mutate(intensity_cumulative = round(intensity_cumulative, 2),
           intensity_max = round(intensity_max, 2))
  if(nrow(focus_event) > 1){
    focus_event <- slice(focus_event, nrow(focus_event))
  }

  # Create vectors for sub-optimal data
  if(full_seq){
    length_seq <- 1982:2009
    missing_seq <- seq(0.00, 0.50, 0.01)
    trend_seq <- seq(0.00, 0.30, 0.01)
  } else{
    length_seq <- seq(1984, 2009, 5)
    missing_seq <- seq(0.00, 0.50, 0.1)
    trend_seq <- seq(0.00, 0.30, 0.1)
  }

  # Sub-optimise data
  ## Length
  sst_length <- plyr::ldply(length_seq, control_length, df = sst_flat)
  ## Missing data
  sst_missing <- plyr::ldply(missing_seq, control_missing, df = sst_flat)
  ## Decadal trend
  sst_trend <- plyr::ldply(trend_seq, control_trend, df = sst_flat)

  # Calculate MHWs in most recent 10 years of data and return the desired clims and metrics
  sst_base_res <- rbind(sst_length, sst_missing, sst_trend) %>%
    group_by(test, index_vals) %>%
    group_modify(~clim_metric_focus_calc(.x, focus_dates = focus_event)) %>%
    ungroup()

  # Run the tests while also interpolating all gaps
  sst_interp_res <- sst_missing %>%
    group_by(test, index_vals) %>%
    group_modify(~clim_metric_focus_calc(.x, focus_dates = focus_event, set_pad = 9999)) %>%
    ungroup() %>%
    mutate(test = "interp")

  # Combine for ease of use
  sst_clim_metric <- data.frame(rbind(sst_base_res, sst_interp_res))

  if(windows){
    sst_windows <- window_test(sst_length, focus_event)
    sst_clim_metric <- rbind(sst_clim_metric, sst_windows)
  }

  # Run ANOVA/Tukey on MHW results for three different tests
  sst_signif <- sst_clim_metric %>%
    filter(id != "focus_event") %>%
    group_by(test, var) %>%
    group_modify(~kruskal_post_hoc(.x))

  # Create summary statistics of MHW results
  sst_summary <- sst_clim_metric %>%
    group_by(test) %>%
    group_modify(~summary_stats(.x)) %>%
    rbind(sst_signif) %>%
    data.frame()

  # Include clim/metric data if requested
  if(clim_metric){
    sst_summary <- rbind(sst_summary, sst_clim_metric)
  }

  # Count consecutive missing days
  # NB: This would be useful information to have for each pixel
  # but it takes too long to calculate
  if(count_miss){
    sst_missing_count <- sst_missing %>%
      group_by(test, index_vals) %>%
      group_modify(~con_miss(.x)) %>%
      dplyr::rename(id = duration,
                    val = count) %>%
      mutate(var = "con_miss",
             id = as.character(id)) %>%
      data.frame()
    sst_summary <- rbind(sst_summary, sst_missing_count)
  }

  # Add decadal trend to end of data.frame and exit
  sst_summary <- rbind(sst_summary, dec_trend) %>%
    mutate(test = as.factor(test)) %>%
    data.frame()
  return(sst_summary)
}


# Global functions --------------------------------------------------------

# The function that runs all of the tests on a single pixel/time series
# nc_file <- OISST_files[1]
global_analysis <- function(nc_file){

  sst <- tidync(nc_file) %>%
    hyper_tibble() %>%
    dplyr::rename(t = time, temp = sst) %>%
    mutate(t = as.Date(t, origin = "1970-01-01")) %>%
    filter(t <= "2018-12-31",
           round(temp, 1) != -1.8) %>%
    na.omit() %>%
    # Filter out pixels that don't cover the whole time series
    # Filter out pixels with any ice cover during the time series
    group_by(lon, lat) %>%
    filter(n() == 13514) %>% # A full time series is 13,514 days long
    ungroup() %>%
    select(lon, lat, t, temp) %>%
    mutate(lon = ifelse(lon > 180, lon-360, lon)) %>%
    data.frame()

  # Run the analysis on each lon/lat pixel
  # system.time(
  res <- plyr::ddply(sst, c("lon", "lat"), single_analysis, .parallel = F, .progress = "text")
  # ) # 121 seconds in parallel, xxx seconds not in parallel
  # ~5 seconds for one pixel, ~30 minutes not in parallel
  # Times will vary by perhaps as much as 50% due to change in pixel count per longitude step
  return(res)
}

# The function that unpacks global resuls as desired
# testers...
# global_file <- global_files[1]
# sub_level <- 1
# sub_level <- 2
# sub_level <- 3
global_unpack_sub <- function(global_file, sub_level){

  # Load the one slice
  slice_res <- readRDS(global_file)

  # Pull it apart as desired
  slice_full <- data.frame()
  for(i in 1:length(slice_res)){
    slice_step_1 <- slice_res[[i]]
    slice_step_coords <- as.vector(strsplit(names(slice_res)[i], "[.]"))
    slice_step_lon <- paste0(slice_step_coords[[1]][1],".",slice_step_coords[[1]][2])
    slice_step_lat <- paste0(slice_step_coords[[1]][3],".",slice_step_coords[[1]][4])
    if(is.null(slice_step_1)){
    } else if(length(slice_step_1[sub_level]) > 0) {
      names(slice_step_1) <- NULL # Prevents column names from being altered
      slice_step_2 <- data.frame(slice_step_1[sub_level],
                                 lon = slice_step_lon,
                                 lat = slice_step_lat)
      if(ncol(slice_step_2) == 3) colnames(slice_step_2)[1] <- "dec_trend"
      slice_full <- rbind(slice_full, slice_step_2)
    }
  }

  # Sort it out and exit
  slice_full <- slice_full %>%
    select(lon, lat, everything()) %>%
    mutate(lon = as.numeric(as.character(lon)),
           lat = as.numeric(as.character(lat)))
  return(slice_full)
}

# The function that crawls through all of the global results
# unpacks them and then stitches them together before saving
global_unpack <- function(){
  # Point to files
  global_files <- dir("data/global", full.names = T)

  ## Run unpacker, save, and clear one at a time to minimise RAM use
  # Decadal trends
  print(paste0("Began unpacking decadal trends at ",Sys.time()))
  global_dec_trend <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                                  .parallel = T, sub_level = 1)
  saveRDS(global_dec_trend, "data/global_dec_trend.Rda")
  rm(global_dec_trend); gc()
  print(paste0("Finished unpacking decadal trends at ",Sys.time()))

  # Summary stats
  print(paste0("Began unpacking summary stats at ",Sys.time()))
  global_summary <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                                .parallel = T, sub_level = 2)
  # colnames(global_effect_event) <- gsub("effect_event.", "", colnames(global_effect_event))
  saveRDS(global_summary, "data/global_summary.Rda")
  rm(global_summary); gc()
  print(paste0("Finished unpacking summary stats at ",Sys.time()))

  # Focus event
  print(paste0("Began unpacking focus event at ",Sys.time()))
  global_focus <- plyr::ldply(.data = global_files, .fun = global_unpack_sub,
                              .parallel = T, sub_level = 3)
  # colnames(global_effect_event) <- gsub("effect_event.", "", colnames(global_effect_event))
  saveRDS(global_focus, "data/global_focus.Rda")
  rm(global_focus); gc()
  print(paste0("Finished unpacking focus event at ",Sys.time()))
}

# Function for calculating R2 values for sub-optimal results
global_slope_sub <- function(df){
  if(nrow(df) > 2){
    round(broom::tidy(lm(val ~ index_vals, data = df))$estimate[2], 3)
  } else{
    return(NA)
  }
}

# This function expects to be given only one latitude slice at a time
# tester...
# df <- global_focus %>%
# filter(lon == lon[20], lat == lat[129])
# filter(lon == lon[20], lat == lat[129], test == "length", var == "focus_duration")
# filter(test == "length", var == "focus_duration")
# df <- global_focus %>%
# filter(lon == lon[20], lat == lat[129])
# filter(lon == lon[20], lat == lat[129], test == "length", var == "duration")
global_slope <- function(df){
  suppressWarnings( # Suppress perfect slope warnings
  df_slope <- df %>%
    group_by(lon, test, var) %>%
    group_modify(~global_slope_sub(.x))
    # nest() %>%
    # mutate(slope = purrr::map(data, global_slope_sub)) %>%
    # select(-data) %>%
    # unnest()
  )
}


# Figure functions --------------------------------------------------------

# Function for rounding decimal places of plot labels
fmt_dcimals <- function(decimals = 0){
  function(x) as.character(round(x, decimals))
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
    # dplyr::mutate(diff = thresh - seas,
                  # thresh_2x = thresh + diff,
                  # thresh_3x = thresh_2x + diff,
                  # thresh_4x = thresh_3x + diff) %>%
    dplyr::filter(t >= date_peak-spread, t <= date_peak+spread)

  # Peak event
  peak_event <- clim_cat %>%
    group_by(site_label) %>%
    filter(t >= date_start-1, t <= date_end+1) %>%
    mutate(start_point = thresh[t == min(t)],
           peak_point = thresh[t == date_peak],
           end_point = thresh[t == max(t)])

  # Set line colours
  lineColCat <- c(
    "Temperature" = "black",
    "Climatology" = "skyblue",
    "Threshold" = "navy")

  # Set category fill colours
  fillColCat <- c(
    "Focus MHW" = "red",
    "Other MHWs" = "salmon"
  )

  ggplot(data = clim_cat, aes(x = t, y = temp)) +
    geom_flame(aes(y2 = thresh, fill = "Other MHWs"), n = 5, n_gap = 2) +
    geom_flame(data = peak_event, aes(y2 = thresh, fill = "Focus MHW"), n = 5, n_gap = 2) +
    geom_line(aes(y = seas, col = "Climatology"), size = 0.7) +
    geom_line(aes(y = thresh, col = "Threshold"), size = 0.7) +
    geom_line(aes(y = temp, col = "Temperature"), size = 0.6) +
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
    scale_colour_manual(name = NULL, values = lineColCat,
                        breaks = c("Climatology", "Threshold", "Temperature")) +
    scale_fill_manual(name = NULL, values = fillColCat,
                      breaks = c("Focus MHW", "Other MHWs")) +
    scale_x_date(date_labels = "%b %Y", expand = c(0, 0)) +
    labs(y = y_label, x = NULL) +
    facet_wrap(~site_label, ncol = 1, scales = "free") +
    theme(legend.position = "top")
}

# The code that creates Figure 2
# This shows the percentage change in sea/thresh and dur/max.int
# from control for the three base tests
# testers...
# df <- sst_ALL_results
fig_line_plot <- function(df = sst_ALL_results){

  # Prep data.frame for pretty plotting
  df_prep <- df %>%
    filter(test %in% c("length", "missing", "trend"),
           var %in% c("seas", "thresh", "duration", "intensity_max")) %>%
    mutate(val = val*100,
           index_vals = ifelse(test == "missing", index_vals*100, index_vals),
           site_label = case_when(site == "WA" ~ "A",
                                  site == "NWA" ~ "B",
                                  site == "Med" ~ "C"),
           var_label = case_when(var == "duration" ~ "Duration (days)",
                                 var == "intensity_max" ~ "Max intensity (°C)",
                                 var == "seas" ~ "Seasonal clim. (°C)",
                                 var == "thresh" ~ "Threshold clim. (°C)"),
           test_label = case_when(test == "length" ~ "Time series length (years)",
                                  test == "missing" ~ "Missing data (%)",
                                  test == "trend" ~ "Trend (°C/dec)"),
           var_label = factor(var_label,
                              levels = c("Seasonal clim. (°C)",
                                         "Threshold clim. (°C)",
                                         "Duration (days)",
                                         "Max intensity (°C)")))

  # Mean values
  df_mean <- df_prep %>%
    filter(id == "mean_perc")

  # Significance points
  df_sig <- df_prep %>%
    filter(id == "difference", val == 1)

  # Create figure
  fig_plot <- ggplot(df_mean, aes(x = index_vals, y = val)) +
    geom_hline(aes(yintercept = 0), colour = "grey") +
    geom_line(aes(colour = site), size = 2, alpha = 0.7) +
    geom_point(data = df_sig, colour = "red", size = 1, alpha = 0.4) +
    scale_colour_brewer(palette = "Dark2") +
    facet_grid(var_label~test_label, scales = "free", switch = "x") +
    labs(y = "Change from control value (%)", x = NULL, colour = "Site") +
    # coord_cartesian(ylim = c(-3, 3)) +
    theme(legend.position = "top")
  fig_plot
  return(fig_plot)
}

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

  # The ridgeline plot
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

# Wrapper for creating plugs for smoother fig 2 - 4 plotting
# control_plug <- function(test_plug, site_plug){
#   plug <- data.frame(test = test_plug, site = site_plug, metric = unique(sst_ALL_plot_long$metric),
#                      index_vals = 30, p.value.min = 1.0, p.value.mean = 1.0,
#                      p.value.max = 1.0, p.value.sd = 0)
# }

# Expects a one row data.frame with a 'lon' and 'lat' column
# df <- focus_WA
# map_point <- function(df){
#   category_data <- readRDS(category_files[grepl(pattern = as.character(df$date_peak),
#                                                 x = category_files)])
#   map_out <- ggplot(data = df, aes(x = lon, y = lat)) +
#     geom_tile(data = category_data, aes(fill = category)) +
#     borders(fill = "grey80", colour = "black") +
#     geom_point(shape = 21, colour = "white", fill = "hotpink", size = 2) +
#     scale_fill_manual("Category",
#                       values = c("#ffc866", "#ff6900", "#9e0000", "#2d0000"),
#                       labels = c("I Moderate", "II Strong",
#                                  "III Severe", "IV Extreme")) +
#     coord_cartesian(xlim = c(df$lon[1]-20, df$lon[1]+20),
#                     ylim = c(df$lat[1]-20, df$lat[1]+20)) +
#     labs(x = NULL, y = NULL) +
#     theme(legend.position = "bottom")
#   return(map_out)
# }

# This function expects the output of the clim, event, cat pipe
# df <- filter(sst_ALL_res, site == "NW_Atl")
# table_summary <- function(site_1){
#
#   df <- filter(sst_ALL_res, site == site_1)
#
#   # Seasonal min/mean/max
#   # Threshold min/mean/max
#   clim_long <- df %>%
#     select(clims) %>%
#     unnest() %>%
#     select(doy, seas, thresh) %>%
#     unique() %>%
#     select(seas, thresh) %>%
#     gather(key = "metric", value = "val")
#
#   event_long <- df %>%
#     select(events) %>%
#     unnest() %>%
#     filter(row_number() %% 2 == 0) %>%
#     unnest() %>%
#     select(duration, intensity_mean, intensity_max, intensity_cumulative) %>%
#     gather(key = "metric", value = "val")
#
#   event_count <- df %>%
#     select(cats) %>%
#     unnest() %>%
#     select(category) %>%
#     group_by(category) %>%
#     summarise(val = n()) %>%
#     mutate(metric = "count") %>%
#     spread(key = category, value = val)
#   if(!"IV Extreme" %in% colnames(event_count)){
#     event_count <- cbind(event_count, tibble('IV Extreme' = 0))
#   }
#   event_count <- event_count %>%
#     select(metric, 'I Moderate', 'II Strong', 'III Severe', 'IV Extreme') %>%
#     dplyr::rename(' ' = metric)
#
#   summary_res <- rbind(clim_long, event_long) %>%
#     group_by(metric) %>%
#     summarise_all(.funs = c("min", "mean", "max", "sd")) %>%
#     mutate_if(is.numeric, round, 2) %>%
#     arrange(metric)
#
#   tbl_1 <- gridExtra::tableGrob(summary_res, rows = NULL)
#   tbl_2 <- gridExtra::tableGrob(event_count, rows = NULL)
#
#   tbl_all <- gridExtra::grid.arrange(tbl_1, tbl_2,
#                                      nrow = 2, as.table = TRUE)
#
#   return(tbl_all)
#   # res_list <- list(summary_res = summary_res,
#                    # event_count = event_count)
#   # return(res_list)
# }

# Wrapper for time series + clims + event rug
# ts_clim_rug <- function(site_1){
#   ts_data <- filter(sst_ALL, site == site_1)
#   clim_data <- filter(sst_ALL_clim, site == site_1)
#   event_data <- filter(sst_ALL_event, site == site_1)
#   fig <- ggplot(data = ts_data, aes(x = t, y = temp)) +
#     geom_line(colour = "grey20") +
#     geom_line(data = clim_data, aes(y = seas),
#               linetype = "dashed", colour = "steelblue3") +
#     geom_line(data = clim_data, linetype = "dotted", colour = "tomato3",
#               aes(x = t, y = thresh)) +
#     geom_rug(data = event_data, sides = "b", colour = "red3", size = 2,
#              aes(x = date_peak, y = min(ts_data$temp))) +
#     labs(x = NULL, y = "Temperature (°C)")
#   return(fig)
# }

# Wrapper for clim only line plot
# clim_line <-function(site_1){
#   clim_data <- filter(sst_ALL_clim_only, site == site_1)
#   ggplot(data = clim_data, aes(x = doy)) +
#     geom_line(aes(y = seas), colour = "steelblue3") +
#     geom_line(aes(y = thresh), colour = "tomato3") +
#     labs(x = NULL, y = "Temperature (°C)")
# }

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
