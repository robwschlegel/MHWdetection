# code/functions.R
# This script contains all of the functions used throughout 'code/workflow.R'


# Libraries ---------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))

# library(jsonlite, lib.loc = "~/R-packages/")
# library(dplyr, lib.loc = "~/R-packages/") # Development version for group_modify()
library(tidyverse, lib.loc = "~/R-packages/")
# library(ggridges)
# library(broom)
library(heatwaveR, lib.loc = "~/R-packages/")
library(data.table, lib.loc = "~/R-packages/")
# cat(paste0("heatwaveR version = ",packageDescription("heatwaveR")$Version))
library(lubridate) # This is intentionally activated after data.table
# library(fasttime)
library(ggpubr)
library(boot)
# library(FNN)
# library(mgcv)
library(doMC); registerDoMC(cores = 50)
library(tidync, lib.loc = "~/R-packages/")
# library(rgdal)
# library(pgirmess)
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
OISST_files <- dir(path = "~/data/OISST", full.names = T, pattern = "avhrr")

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


# Load OISST --------------------------------------------------------------

# Load an OISST NetCDF file and exclude pixels with any ice (-1.8C) cover
# During iterations of this methodology it was found that there are pixels near
# the ice edge that don't ever quite reach -1.8C but remain in a "slush" state
# This was determined to be an artefact of the OISST process and so the new
# "ice" limit was set to -1.6C

# testers...
# An Antarctic pixel that doesn't have a temp = -1.8 -"-78.375	-76.125"
# which(c(seq(0.125, 179.875, by = 0.25), seq(-179.875, -0.125, by = 0.25)) == -78.375) # 1127
# nc_file <- OISST_files[1127]
load_noice_OISST <- function(nc_file){
  # suppressWarnings( # (2019-10-16) tidync started giving a meaningless warning message
  res <- tidync(nc_file) %>%
    hyper_filter(lat = between(lat, -70, 80)) %>% # No need to load the very poleward data
    hyper_tibble() %>%
    dplyr::rename(t = time, temp = sst) %>%
    mutate(t = as.Date(t, origin = "1970-01-01")) %>%
    filter(t <= "2018-12-31",
           round(temp, 1) > -1.6) %>%
    na.omit() %>%
    # Filter out pixels that don't cover the whole time series
    # Filter out pixels with any ice/slush cover during the time series
    group_by(lon, lat) %>%
    filter(n() == 13514) %>% # A full time series is 13,514 days long
    ungroup() %>%
    select(lon, lat, t, temp) %>%
    mutate(lon = ifelse(lon > 180, lon-360, lon)) %>%
    data.frame()
  # )
  return(res)
}

# test <- res %>%
#   filter(lat == -76.125)


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
# df <- filter(sst_length, index_vals == 25)
# df <- filter(sst_interp, index_vals == 0.2)
# df <- filter(sst_window_10, index_vals == 10)
# df <- filter(sst_missing, index_vals == 0.37)
# df <- sst_missing[sst_missing$index_vals == 0.36,]
# df <- filter(sst_trend, index_vals == 0.20)
# set_window = 5
# set_window = 20
# set_pad = F
# min_date = "2009-01-01"
# year_start = 0
# year_end = 0
# focus_dates = focus_event
clim_metric_focus_calc <- function(df, set_window = 5, set_pad = F, min_date = "2009-01-01",
                                   year_start = 0, year_end = 0, focus_dates){

  # First and last years for full clim period
  if(year_start == 0) year_start <- min(lubridate::year(df$t))
  if(year_end == 0) year_end <- max(lubridate::year(df$t))

  # base calculation
  res <- ts2clm(df, windowHalfWidth = set_window, maxPadLength = set_pad, var = TRUE,
                climatologyPeriod = c(paste0(year_start,"-01-01"), paste0(year_end,"-12-31"))) %>%
    filter(t >= min_date) %>%
    detect_event()

  # Occasionally a large added trend will cause the focus event to become too large for the code to recognise
  # So we must invert the search pattern and run it again
  res_check <- res$event %>%
    filter(date_peak >= focus_dates$date_start,
           date_peak <= focus_dates$date_end)
  if(nrow(res_check) == 0){
    res_check <- res$event %>%
      filter(date_start <= focus_dates$date_start,
             date_end >= focus_dates$date_end)
  }

  # The metrics for the event(s) occurring during the largest control event
  res_focus <- res_check %>%
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
#


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


# Summary stats -----------------------------------------------------------

# Summary stats for broad results
# tester...
# df <- sst_clim_metric %>%
# ungroup() %>%
# filter(test == "window_20") %>%
# select(-test)
summary_stats <- function(df){

  # Determine the control group
  if(30 %in% df$index_vals){
    control_val <- 30
  } else{
    control_val <- 0
  }

  # Create gapless index of test values
  # This is necessary for the missing data as some steps have no MHWs
  index_step <- unique(df$index_vals)[2]-unique(df$index_vals)[1]
  if(control_val == 0){
    if(max(df$index_vals) > 0.3){
      gapless_vector <- seq(0, 0.5, by = index_step)
    } else{
      gapless_vector <- seq(0, 0.3, by = index_step)
    }
  } else if(control_val == 30){
    gapless_vector <- seq(10, max(unique(df$index_vals)), by = index_step)
  }
  index_gapless <- data.frame(index_vals = gapless_vector)

  # Count of control values
  control_count <- df %>%
    filter(index_vals == control_val,
                var == "duration") %>%
    nrow()

  # Count of values
  res_count <- df %>%
    group_by(index_vals) %>%
    count(var) %>%
    filter(var == "duration") %>%
    select(-var) %>%
    unique() %>%
    right_join(index_gapless, by = "index_vals") %>%
    mutate(n = replace_na(n, 0)) %>%
    mutate(var = "count",
           n_diff = n-control_count,
           n_perc = (n-control_count)/abs(control_count)*100) %>%
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

  # Find proportion of change
  res_perc <- left_join(res_base, res_control, by = c("var", "id")) %>%
    mutate(perc = (val-cont_val)/abs(cont_val)*100,
           perc = replace_na(perc, 0),
           perc = if_else(is.infinite(perc), -1, perc), # Caused by non-existent focus MHW
           perc_name = paste0(id,"_perc")) %>%
    select(index_vals, var, perc_name, perc) %>%
    dplyr::rename(id = perc_name, val = perc)

  # Combine, round, and exit
  res <- rbind(res_base, res_perc, res_count) %>%
    mutate(val = round(val, 3)) %>%
    arrange(index_vals, var, id)
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


# Full analyses -----------------------------------------------------------

# This single function runs through and outputs all of the desired tests as a list
# testers...
# df <- sst_Med
# df <- sst[sst$lat == sst$lat[1],]
# full_seq = T
# clim_metric = F
# count_miss = T
# windows = T
# Bad pixels for testing - "140.375 0.625", "-73.625_-77.125", ", "-112.125 -28.875"(window width),
# "-149.375 10.625"(focus event dissapears with large decadal trend value)
# which(c(seq(0.125, 179.875, by = 0.25), seq(-179.875, -0.125, by = 0.25)) == 34.625)
# df <- load_noice_OISST(OISST_files[139]) %>%
# filter(lat == -33.875)
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
  if(nrow(focus_event) == 0) return() # Some nearly prozen pixels manage to slip through and cause issues
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
  # NB: The warnings that pop up here are fine and are caused by missing data completely
  # removing all of the MHWs from a time series. This is accounted for in the summary stats
  # system.time(
  sst_base_res <- rbind(sst_length, sst_missing, sst_trend) %>%
    group_by(test, index_vals) %>%
    group_modify(~clim_metric_focus_calc(.x, focus_dates = focus_event)) %>%
    ungroup()
  # ) # 21 seconds, extensive testing of data.table was not faster

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

  # Create summary statistics of MHW results
  sst_summary <- sst_clim_metric %>%
    group_by(test) %>%
    group_modify(~summary_stats(.x)) %>%
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


# This function runs the full analysis on a randomly selected pixel
  # NB: empty_integer is here just to satisfy plyr::ldply()
random_analysis <- function(empty_integer, base_period = F){

  # Load a random lon slice
  sst <- load_noice_OISST(sample(OISST_files, 1))

  # Randomly pick a lat slice
  pixel <- filter(sst, lat == sample(sst$lat, 1))
  rm(sst); gc()

  # Make the function more chatty for error trapping
  # print(paste0(pixel$lon[1],"_",pixel$lat[1]))

  # Run a full analysis and exit
  if(base_period){
    res <- base_period_analysis(pixel, clim_metric = F) %>%
      mutate(lon = pixel$lon[1],
             lat = pixel$lat[1]) %>%
      select(lon, lat, everything())
  } else{
    res <- single_analysis(pixel, full_seq = T, clim_metric = F, count_miss = T, windows = T) %>%
      mutate(lon = pixel$lon[1],
             lat = pixel$lat[1]) %>%
      select(lon, lat, everything())
  }
  return(res)
}


# This function runs the analysis on a full lon slice
# nc_file <- OISST_files[0040]
# par_op = T
global_analysis <- function(nc_file, par_op = F){

  # Load and prep data
  sst <- load_noice_OISST(nc_file)

  # Run the analysis on each lon/lat pixel
  # system.time(
  res <- plyr::ddply(sst, c("lon", "lat"), single_analysis, .parallel = par_op, .progress = "text")
  # ) # 2 minutes in parallel, 30 minutes not in parallel
  # ~5 seconds for one pixel, ~xxx seconds not in parallel
  # Times may vary by 50% due to change in pixel count per longitude step
  return(res)
}


# Trend calculations ------------------------------------------------------

# Fills in gaps in data
# This is expected to be fed a data.frame with only and `index_vals` and `val` column
# This reduces the need to create larger more complex data.frames when plugging holes
gapless_data <- function(df){
  # Determine the control group
  if(30 %in% df$index_vals){
    control_val <- 30
  } else{
    control_val <- 0
  }

  # Create gapless index of test values
  # This is necessary for the missing data as some steps have no MHWs
  index_step <- unique(df$index_vals)[2]-unique(df$index_vals)[1]
  if(control_val == 0){
    if(max(df$index_vals) > 0.3){
      gapless_vector <- seq(0, 0.5, by = index_step)
    } else{
      gapless_vector <- seq(0, 0.3, by = index_step)
    }
  } else if(control_val == 30){
    gapless_vector <- seq(10, max(unique(df$index_vals)), by = index_step)
  }
  index_gapless <- data.frame(index_vals = gapless_vector)

  # Gapless data.frame
  df_full <- df %>%
    right_join(index_gapless, by = "index_vals") %>%
    mutate(val = replace_na(val, 0))
  return(df_full)
}

# Custom function for trend calculation
trend_stats <- function(df){
  # This function must first ensure that any gaps in the data are filled in explicitly
    res_model <- lm(val ~ index_vals, data = df)
    res <- data.frame(intercept = round(broom::tidy(res_model)$estimate[1], 4),
                      slope = round(broom::tidy(res_model)$estimate[2], 4),
                      R2 = round(broom::glance(res_model)$adj.r.squared, 2),
                      p = round(broom::tidy(res_model)$p.value[2], 2)) %>%
      mutate(R2 = ifelse(R2 < 0, 0, R2),
             R2 = ifelse(slope == 0 & intercept == 0, 1, R2),
             p = ifelse(slope == 0, 1, p))
  return(res)
}

# The function that will calculate the slopes correctly based on the name of the test
# Custom function for trend calculation
# Note: The slopes are in units of 1 by R default
# So this makes sense for the length slope being in units of 1 year
# But not for missing data, where 1 is 100%,
# or for decadal trend where 1 is 1C/dec.
# Therefore the missing/interp and decadal trend slopes need to be devided by 100
# This produces slopes that represent 1% and 0.01C/dec changes
# testers...
# full_seq = T
# full_seq = F
# df <- random_quant %>%
#   filter(test == "missing",
#          var == "intensity_max",
#          id == "mean_perc") %>%
#   gather(key = "stat", value = "val", -c(test:id)) %>%
#   mutate(test2 = test) %>%
#   filter(stat == "q50") %>%
#   select(-test2, -var, -id, -stat)
# df <- lon_correct %>%
#   filter(lat == 40.375,
#          test == "missing",
#          var == "count",
#          id == "n_perc")
trend_correct <- function(df, full_seq = T){
  if(df$test[1] == "length"){
    df_A <- df %>%
      filter(index_vals %in% 10:30) %>%
      mutate(index_vals = abs(index_vals-30))
    res_A <- trend_stats(df_A) %>%
      mutate(range = "30 - 10")
    if(full_seq){
      df_B <- df %>%
        filter(index_vals %in% 30:37) %>%
        mutate(index_vals = abs(index_vals-30))
      res_B <- trend_stats(df_B) %>%
        mutate(range = "30 - 37")
      res <- rbind(res_A, res_B)
    } else{
     res <- res_A
    }
  } else if(df$test[1] == "missing" & df$var[1] %in% c("count", "focus_count")){
    if(full_seq){
      miss_seq_A <- seq(0, 0.25, by = 0.01)
      miss_seq_B <- seq(0.26, 0.50, by = 0.01)
    } else{
      miss_seq_A <- seq(0, 0.3, by = 0.1)
      miss_seq_B <- seq(0.3, 0.5, by = 0.1)
    }
    df_A <- df %>%
      filter(index_vals <= max(miss_seq_A)+0.001) # There is some weird float rounding hapening off screen...
    df_B <- df %>%
      filter(index_vals  >= max(miss_seq_A)-0.001)
    res_A <- trend_stats(df_A) %>%
      mutate(range = "0.00 - 0.25")
    res_B <- trend_stats(df_B) %>%
      mutate(range = "0.26 - 0.50")
    res <- rbind(res_A, res_B)
  } else{
    res <- trend_stats(df) %>%
      mutate(range = paste0(min(df$index_vals),".00 - ",max(df$index_vals),"0"))
  }
  if(df$test[1] != "length"){
    res$slope <- res$slope/100 # linear models default slope output is in steps of 1
  }
  res$slope <- round(res$slope, 2)
  return(res)
}

# A convenience function to load a lon slice of global results
# Then filter out the variables not used in the results and run linear models
# file_name <- "data/global/test_1130.Rda"
# file_name <- "data/global/slice_0050.Rda"
var_trend <- function(file_name){
  # The subsetting index
  var_choice <- data.frame(var = c("count", "duration", "intensity_max",
                                   "focus_count", "focus_duration", "focus_intensity_max"),
                           id = c("n", "sum_perc", "mean_perc",
                                  "mean_perc", "sum_perc", "mean_perc"),
                           stringsAsFactors = F)
  # Still using an older run of the global output
  # So need to calculate n_perc manually here first
  lon_slice <- readRDS(file_name) %>%
    right_join(var_choice, by = c("var", "id")) %>%
    group_by(lon, lat, test, var, id) %>%
    group_modify(~gapless_data(.x))

  # Find the control count values
  control_count <- lon_slice %>%
    filter(var == "count",
           index_vals %in% c(0, 30)) %>%
    dplyr::rename(control_val = val) %>%
    select(-index_vals)

  # Calculate n_perc
  n_perc_df <- lon_slice %>%
    filter(var == "count") %>%
    left_join(control_count, by = c("lon", "lat", "test", "var", "id")) %>%
    ungroup() %>%
    mutate(n_perc = round((val-control_val)/abs(control_val)*100, 3),
           id = "n_perc") %>%
    dplyr::select(-val, -control_val) %>%
    dplyr::rename(val = n_perc)

  # Stick it back together
  lon_correct <- rbind(data.frame(lon_slice), n_perc_df) %>%
    filter(id != "n") %>%
    mutate(val = ifelse(is.infinite(val), -100, val))

  # system.time(
    res <- lon_correct %>%
      mutate(test2 = test, var2 = var) %>%
      group_by(lon, lat, test2, var2, id) %>%
      group_modify(~trend_correct(.x, full_seq = F)) %>%
      dplyr::rename(test = test2, var = var2)
    # ) # 60 seconds
  return(res)
}


# Figure functions --------------------------------------------------------

# Function for rounding decimal places of plot labels
fmt_dcimals <- function(decimals = 0){
  function(x) as.character(round(x, decimals))
}


# Bootstrap CI calculations for Figure 2 & 3
# NB: The problem with this is that we DONT want the CI for the mean
# We want the CI for the range of values
custom_boot <- function(df){
  res <- data.frame(boot.ci(boot(data = df$val, function(x,i) mean(x[i]), R = 1000),
                            type = "perc")$percent) %>%
    select(V4, V5) %>%
    dplyr::rename(lower = V4, upper = V5)
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


# The code that creates most other figures
# This shows the percentage change in sea/thresh and dur/max.int
# from control for any of the tests
# testers...
# df = full_results
# tests  = "base"
# tests  = "miss_comp"
# result_choice = "average"
# result_choice = "focus"

fig_box_plot <- function(df = full_results, tests, result_choice){

  # Choose if the figure will show the base tests or the interp. comp.
  if(tests == "base"){
    test_choice <- c("length", "missing", "trend")
    test_levels <- c("Time series length (years)", "Missing data (%)", "Trend (°C/dec)")
  } else if(tests == "miss_comp"){
    test_choice <- c("missing", "interp")
    test_levels <- c("Missing data (%)", "Interpolated data (%)")
  } else if(tests == "windows"){
    test_choice <- c("length", "window_10", "window_20", "window_30")
    test_levels <- c("Time series length (years)", "Window width 10", "Window width 20", "Window width 30")
  } else if(tests == "base_period"){
    test_choice <- "base_period"
    test_levels <- "Difference from WMO base period (years)"
  } else{
    stop("Provide a valid 'tests' argument")
  }

  # Prep data for focus or mean MHW results
  if(result_choice == "focus"){
    var_choice <- data.frame(var = c("focus_count", "focus_duration", "focus_intensity_max"),
                             id = c("mean_perc", "sum_perc", "mean_perc"),
                             stringsAsFactors = F)
    var_levels <- c("Count (% n)", "Duration (% sum of days)", "Max. intensity (% of mean °C)")
    y_axis_title <- "Change in focus MHW"
  } else if(result_choice == "average"){
    var_choice <- data.frame(var = c("count", "duration", "intensity_max"),
                             id = c("n_perc", "sum_perc", "mean_perc"),
                             stringsAsFactors = F)
    var_levels <- c("Count (% n)", "Duration (% sum of days)", "Max. intensity (% of mean °C)")
    y_axis_title <- "Change in average MHWs"
  } else if(result_choice == "clims"){
    var_choice <- data.frame(var = c("seas", "thresh"),
                             id = c("sd_perc", "mean_perc"),
                             stringsAsFactors = F)
    var_levels <- c("Seasonal clim. (SD; °C)", "Threshold clim. (mean; °C)")
    y_axis_title <- "Change in thresholds"
  } else{
    stop("Provide a valid 'result_choice' argument")
  }

  # Manually remove a couple of rediculous pixels
  # NB: These pixels were determined in the first run of the code to be anomalous
  # due to the structure of the time series and should have been filtered earlier
  # bad_pixels <- c("140.375_0.625", "-73.625_-77.125", "-44.375_1.125", "-149.375_10.625")
  # pixel_hunt <- filter(df,
  #                test == "trend",
  #                index_vals == 0.20,
  #                var == "focus_duration",
  #                id == "mean_perc",
  #                val < -90)

  # Prep reference results for pretty plotting
  df_prep <- df %>%
    filter(test %in% test_choice) %>% #,
           # is.finite(val),
           #!(site %in% bad_pixels)) %>%
    right_join(var_choice, by = c("var", "id")) %>%
    mutate(index_vals = ifelse(test %in% c("missing", "interp"), index_vals*100, index_vals),
           var_label = case_when(var %in% c("duration", "focus_duration") ~ "Duration (% sum of days)",
                                 var %in% c("intensity_max", "focus_intensity_max" ) ~ "Max. intensity (% of mean °C)",
                                 var %in% c("count", "focus_count") ~ "Count (% n)",
                                 var == "seas" ~ "Seasonal clim. (SD; °C)",
                                 var == "thresh" ~ "Threshold clim. (mean; °C)"),
           var_label = factor(var_label, levels = var_levels),
           test_label = case_when(test == "length" & "length" %in% test_choice ~ "Time series length (years)",
                                  test == "missing" ~ "Missing data (%)",
                                  test == "trend" ~ "Trend (°C/dec)",
                                  test == "interp" ~ "Interpolated data (%)",
                                  test == "base_period" ~ "Difference from WMO base period (years)",
                                  test == "window_10" ~ "Window width 10",
                                  test == "window_20" ~ "Window width 20",
                                  test == "window_30" ~ "Window width 30"),
           test_label = factor(test_label, levels = test_levels))

  # Set the panel labels
  if(tests == "base"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var %in% c("count", "focus_count") & test == "length" ~ "A",
                                     var %in% c("duration", "focus_duration") & test == "length" ~ "B",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "length" ~ "C",
                                     var %in% c("count", "focus_count") & test == "missing" ~ "D",
                                     var %in% c("duration", "focus_duration") & test == "missing" ~ "E",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "missing" ~ "F",
                                     var %in% c("count", "focus_count") & test == "trend" ~ "G",
                                     var %in% c("duration", "focus_duration") & test == "trend" ~ "H",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "trend" ~ "I"))
  } else if(tests == "miss_comp" & result_choice == "average"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var  == "count" & test == "missing" ~ "A",
                                     var == "duration"& test == "missing" ~ "B",
                                     var == "intensity_max" & test == "missing" ~ "C",
                                     var =="count" & test == "interp" ~ "D",
                                     var == "duration" & test == "interp" ~ "E",
                                     var == "intensity_max" & test == "interp" ~ "F"))
  } else if(tests == "miss_comp" & result_choice == "focus"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var == "focus_count" & test == "missing" ~ "H",
                                     var == "focus_duration" & test == "missing" ~ "I",
                                     var == "focus_intensity_max" & test == "missing" ~ "J",
                                     var == "focus_count" & test == "interp" ~ "K",
                                     var == "focus_duration" & test == "interp" ~ "L",
                                     var == "focus_intensity_max" & test == "interp" ~ "M"))
  } else if(tests == "windows"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var %in% c("count", "focus_count") & test == "length" ~ "A",
                                     var %in% c("duration", "focus_duration") & test == "length" ~ "B",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "length" ~ "C",
                                     var %in% c("count", "focus_count") & test == "window_10" ~ "D",
                                     var %in% c("duration", "focus_duration") & test == "window_10" ~ "E",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "window_10" ~ "F",
                                     var %in% c("count", "focus_count") & test == "window_20" ~ "G",
                                     var %in% c("duration", "focus_duration") & test == "window_20" ~ "H",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "window_20" ~ "I",
                                     var %in% c("count", "focus_count") & test == "window_30" ~ "J",
                                     var %in% c("duration", "focus_duration") & test == "window_30" ~ "K",
                                     var %in% c("intensity_max", "focus_intensity_max") & test == "window_30" ~ "L"))
  } else if(tests == "base_period" & result_choice == "average"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var %in% c("count", "focus_count") ~ "A",
                                     var %in% c("duration", "focus_duration") ~ "B",
                                     var %in% c("intensity_max", "focus_intensity_max") ~ "C"))
  } else if(tests == "base_period" & result_choice == "focus"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var %in% c("count", "focus_count") ~ "D",
                                     var %in% c("duration", "focus_duration") ~ "E",
                                     var %in% c("intensity_max", "focus_intensity_max") ~ "F"))
  }
  if(result_choice == "clims"){
    df_prep <- df_prep %>%
      mutate(panel_label = case_when(var %in% c("seas") & test == "length" ~ "A",
                                     var %in% c("thresh") & test == "length" ~ "B",
                                     var %in% c("seas") & test == "missing" ~ "C",
                                     var %in% c("thresh") & test == "missing" ~ "D",
                                     var %in% c("seas") & test == "trend" ~ "E",
                                     var %in% c("thresh") & test == "trend" ~ "F"))
  }

  # Mean values
  reference_df <- df_prep %>%
    filter(site %in% c("WA", "NWA", "Med"),
           index_vals != 50) # For better plotting without a lot of fuss...
  random_df <- df_prep %>%
    filter(!(site %in% c("WA", "NWA", "Med")),
           index_vals != 50)

  # Upper and lower quantile values
  quant_df <- df_prep %>%
    group_by(test, index_vals, var, id, test_label, var_label, panel_label) %>%
    summarise(q05 = quantile(val, 0.05),
              q25 = quantile(val, 0.25),
              q50 = quantile(val, 0.50),
              q75 = quantile(val, 0.75),
              q95 = quantile(val, 0.95),
              iqr50 = q75-q25,
              iqr90 = q95-q05) %>%
    ungroup() %>%
    filter(index_vals != 50)

  # Dummy data.frame to ensure that plot labels are drawn at the correct height
  dummy_height <- reference_df %>%
    group_by(var) %>%
    summarise(max_val = max(val, na.rm = T))

  # Plot labels
  labels_df <- quant_df %>%
    group_by(test) %>%
    mutate(panel_label_x = quantile(index_vals, 0.05)) %>%
    ungroup() %>%
    group_by(var) %>%
    mutate(panel_label_y = max(q95)*0.99) %>%
    # mutate(panel_label_y = ifelse(max(upper)*0.99 > ) %>%
    ungroup() %>%
    select(test, var, panel_label_x, panel_label_y, panel_label, test_label, var_label) %>%
    unique() %>%
    left_join(dummy_height, by = c("var")) %>%
    mutate(panel_label_y = ifelse(max_val > panel_label_y, max_val*0.99, panel_label_y),
           panel_label_y = round(panel_label_y, 2)) %>%
    select(-max_val)

  # Coefficients of determination
  # coef_det <- df_prep %>%
    # group_by(test, index_vals, var, id, test_label, var_label)

  # geom_crosbar is not very clever, so we need to help it along with a logic gate...
  if("trend" %in% quant_df$test){
    # Create figure
    fig_plot_base <- ggplot(quant_df, aes(x = index_vals, y = val)) +
      # 90 CI crossbars
      # Need different lines for tests due to the different x-axis interval sizes
      geom_crossbar(data = filter(quant_df, test != "trend"),
                    aes(x = index_vals, y = 0, ymin = q05, ymax = q95),
                    fatten = 0, fill = "grey70", colour = NA, width = 1) +
      geom_crossbar(data = filter(quant_df, test == "trend"),
                    aes(x = index_vals, y = 0, ymin = q05, ymax = q95),
                    fatten = 0, fill = "grey70", colour = NA, width = 0.01) +
      # IQR Crossbars
      geom_crossbar(data = filter(quant_df, test != "trend"),
                    aes(x = index_vals, y = 0, ymin = q25, ymax = q75),
                    fatten = 0, fill = "grey50", width = 1) +
      geom_crossbar(data = filter(quant_df, test == "trend"),
                    aes(x = index_vals, y = 0, ymin = q25, ymax = q75),
                    fatten = 0, fill = "grey50", width = 0.01) +
      # Median segments
      geom_crossbar(data = filter(quant_df, test != "trend"),
                    aes(x = index_vals, y = 0, ymin = q50, ymax = q50),
                    fatten = 0, fill = NA, colour = "black", width = 1) +
      geom_crossbar(data = filter(quant_df, test == "trend"),
                    aes(x = index_vals, y = 0, ymin = q50, ymax = q50),
                    fatten = 0, fill = NA, colour = "black", width = 0.01)
  } else{
    # Create figure
    fig_plot_base <- ggplot(quant_df, aes(x = index_vals, y = val)) +
      # 90 CI crossbars
      geom_crossbar(data = quant_df, aes(x = index_vals, y = 0, ymin = q05, ymax = q95),
                    fatten = 0, fill = "grey70", colour = NA, width = 1) +
      # IQR Crossbars
      geom_crossbar(data = quant_df, aes(x = index_vals, y = 0, ymin = q25, ymax = q75),
                    fatten = 0, fill = "grey50", width = 1) +
      # Median segments
      geom_crossbar(data = quant_df, aes(x = index_vals, y = 0, ymin = q50, ymax = q50),
                    fatten = 0, fill = NA, colour = "black", width = 1)
  }


  # Finish off the figure
  fig_plot <- fig_plot_base +
    # Random results
    # geom_line(data = random_df, aes(group = site), size = 0.5, alpha = 0.05) +
    # geom_boxplot(data = random_df, aes(group = index_vals), outlier.colour = NA) +

    # Reference results
    geom_line(data = reference_df, aes(colour = site), size = 1.0, alpha = 0.5) +

    # Horizontal itercept at 0
    geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +

    # Labels, scales, facets, and themes
    geom_label(data = labels_df,
               aes(label = panel_label, y = panel_label_y, x = panel_label_x)) +
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(limits = c(-100, 100)) +
    # coord_cartesian(ylim = c(-100, 100)) +
    facet_grid(var_label~test_label, scales = "free", switch = "both") +
    labs(y = y_axis_title, x = NULL, colour = "Site") +
    # coord_cartesian(ylim = c(-3, 3)) + # Correct for extreme line overalys messing up labels by coord-ing to the 90CI range
    theme(legend.position = "top")
  # fig_plot
  return(fig_plot)
}


# Function for easily plotting subsets from the global slope results
# testers...
# test_sub <- "length"
# var_sub <- "intensity_max"
# var_sub <- "duration"
# var_sub <- "count"
# var_sub <- "focus_count"
trend_plot <- function(test_sub, var_sub,
                       df = global_var_trend) {

  if(!exists("global_var_trend")) global_var_trend <- readRDS("data/global_var_trend.Rda")

  # Filter base data
  base_sub <- df %>%
    filter(test == test_sub, var == var_sub) %>%
    # correct focus_count trend back to count from percentage
    mutate(trend = ifelse(var == "focus_count", trend/100, trend))

  # Prepare legend title bits
  sen_change <- "Percent change "

  # if(prop){
    # sen_change <- "Percent change "
    # if(!(var_sub %in% c("count", "focus_count"))){
      # base_sub$trend <- base_sub$trend * 100
    # }
  # } else{
    # sen_change <- "Change "
  # }
  if(test_sub == "length"){
    sen_test <- "per year"
  } else if(test_sub == "missing"){
    sen_test <- "per 1% missing data"
  } else if(test_sub == "trend"){
    sen_test <- "per 0.1°C/dec trend"
  }

  # Prepare colour palette
  if(var_sub %in% c("duration", "focus_duration")) {
    vir_op <- "C"
    col_split <- c("purple", "forestgreen")
    sen_var <- "\nin duration (days)"
  } else if(var_sub %in% c("count", "focus_count")) {
    vir_op <- "B"
    col_split <- c("turquoise4", "chocolate")
    sen_var <- "\nin count (n)"
    sen_change <- "Change " # Intentional overwrite of above value
  } else if(var_sub %in% c("intensity_max", "focus_intensity_max")){
    vir_op <- "A"
    col_split <- c("blue", "red")
    sen_var <- "\nin max. intensity (°C)"
  }

  # Find quantiles
  trend_quantiles <- quantile(base_sub$trend, na.rm = T,
                              probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0))

  # Create legend break labels
  # break_labels <- as.numeric(trend_quantiles[2:6])

  # Correct base data to quantiles as the tails are very long
  base_quantile <- base_sub %>%
    mutate(trend = case_when(trend > trend_quantiles[6] ~ trend_quantiles[6],
                             trend < trend_quantiles[2] ~ trend_quantiles[2],
                             trend <= trend_quantiles[6] | trend >= trend_quantiles[2] ~ trend))

  # The map
  trend_map <- ggplot(base_quantile, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = trend)) +
    geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
    # Test latitudes
    # geom_hline(aes(yintercept = 80)) +
    # geom_hline(aes(yintercept = -70)) +
    # scale_fill_viridis_c(option = vir_op) +
    # scale_fill_gradientn(colors = scales::viridis_pal(option = vir_op)(9),
    # limits = c(as.numeric(trend_quantiles[2]),
    # as.numeric(trend_quantiles[4])),
    # breaks = c(as.numeric(trend_quantiles[2:6]))) +
    scale_fill_gradient2(low = col_split[1], high = col_split[2],
                         breaks = c(round(as.numeric(trend_quantiles[2:6]), 2))) +
    coord_equal(expand = F) +
    theme_void() +
    labs(fill = paste0(sen_change, sen_test, sen_var)) +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "cm"),
          panel.background = element_rect(fill = "grey80"))
  # trend_map

  # ggsave(trend_map,
  #        filename = paste0("output/",test_sub,"_",var_sub,"_",type_sub,"_plot.png"), height = 6, width = 10)
  return(trend_map)
}


# Visualise the consecutive count of missing days in a time series


# Supplementary 2 ---------------------------------------------------------

# A convenience wrapper to use a vector to ply through a range of base periods
control_base <- function(year_spread, df, focus_event){
  res <- clim_metric_focus_calc(df, focus_dates = focus_event,
                                year_start = 1982 + year_spread, year_end = 2011 + year_spread) %>%
    mutate(index_vals = year_spread)
  return(res)
}

# This function compares results for different 30 year base periods
base_period_analysis <- function(df, clim_metric = F){

  # Detrend the ts
  sst_flat <- detrend(df)

  # Calculate MHWs from detrended ts
  # NB: Note that the proper WMO base period is used here
  sst_flat_MHW <- detect_event(ts2clm(sst_flat, climatologyPeriod = c("1982-01-01", "2011-12-31")))

  # Pull out the largest event in the ts
  focus_event <- sst_flat_MHW$event %>%
    filter(date_start >= "2009-01-01") %>%
    filter(intensity_cumulative == max(intensity_cumulative)) %>%
    select(event_no, date_start:date_end, duration, intensity_cumulative, intensity_max) %>%
    mutate(intensity_cumulative = round(intensity_cumulative, 2),
           intensity_max = round(intensity_max, 2))
  if(nrow(focus_event) == 0) return() # Some nearly prozen pixels manage to slip through and cause issues
  if(nrow(focus_event) > 1){
    focus_event <- slice(focus_event, nrow(focus_event))
  }

  # Calculate MHWs in most recent 10 years of data and return the desired clims and metrics
  sst_clim_metric <- plyr::ldply(0:7, control_base, df = sst_flat, focus_event = focus_event) %>%
    mutate(index_vals = index_vals + 30) # Correct this to work with the project wide standard

  # Create summary statistics of MHW results
  sst_summary <- sst_clim_metric %>%
    # group_by(test) %>%
    group_modify(~summary_stats(.x)) %>%
    data.frame() %>%
    filter(index_vals >= 30) %>% # Convert back to the index_val exception used for these data
    mutate(index_vals = index_vals - 30)

  # Include clim/metric data if requested
  if(clim_metric){
    sst_summary <- rbind(sst_summary, sst_clim_metric)
  }

  # Exit
  return(sst_summary)
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
