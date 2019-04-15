# code/functions.R
# This script contains all of the functions used throughout 'code/workflow.R'


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(broom)
library(heatwaveR, lib.loc = "~/R-packages/")
# cat(paste0("heatwaveR version = ",packageDescription("heatwaveR")$Version))
library(lubridate) # This is intentionally activated after data.table
# library(fasttime)
library(ggpubr)
library(boot)
library(FNN)
library(mgcv)
library(doMC); doMC::registerDoMC(cores = 50)
library(ggridges)
library(ncdf4)
# library(rcompanion)


# Meta-data ---------------------------------------------------------------

MHW_colours <- c(
  "I Moderate" = "#ffc866",
  "II Strong" = "#ff6900",
  "III Severe" = "#9e0000",
  "IV Extreme" = "#2d0000"
)

OISST_files <- dir(path = "~/data/OISST", full.names = T)

category_files <- as.character(dir(path = "~/data/cat_clim", pattern = "cat.clim",
                                   full.names = TRUE, recursive = TRUE))

sst_ALL_coords <- data.frame(site = c("WA", "NW_Atl", "Med"),
                             lon = c(112.5, -67, 9),
                             lat = c(-29.5, 43, 43.5))


# Re-sampling -------------------------------------------------------------

# Function for creating a re-sample from 37 years of data
sample_37 <- function(rep){
  year_filter <- c(1,2,3)
  while(length(unique(year_filter)) != 37){
    year_filter <- factor(sample(1982:2018, size = 37, replace = F))
  }
  res <- data.frame()
  for(i in 1:length(year_filter)){
    df_sub <- sst_ALL %>%
      filter(year(t) == year_filter[i],
             # Remove leap-days here
             paste0(month(t),"-",day(t)) != "2-29") %>%
      mutate(year = seq(1982,2018)[i],
             month = month(t),
             day = day(t),
             year_orig = year(t)) %>%
      mutate(t = as.Date(fastPOSIXct(paste(year, month, day, sep = "-")))) %>%
      select(-year, -month, -day)
    res <- rbind(res, df_sub)
  }
  res$rep <- as.character(rep)
  return(res)
}


# De-trending -------------------------------------------------------------

detrend <- function(df){
  resids <- broom::augment(lm(temp ~ t, df))
  res <- df %>%
    mutate(temp = temp - resids$.fitted) %>%
    select(-year_orig)
  return(res)
}


# Random knockout ---------------------------------------------------------

# Function for knocking out data but maintaing the time series consistency
random_knockout <- function(prop, df = sst_ALL_flat){
  # NB: Don't allow sampling of first and last value to ensure
  # all time series are the same length
  ts_length <- nrow(filter(df, site == "WA"))
  miss_index <- sample(seq(2, ts_length-1, 1), ts_length*prop, replace = F)
  # df$x1[sample(nrow(df),250)] <- NA
  res <- df %>%
    group_by(site, rep) %>%
    mutate(row_index = 1:n(),
           temp = replace(temp, which(row_index %in% miss_index), NA)) %>%
    mutate(index_vals = prop) %>%
    select(-row_index)
  return(res)
}


# Count consecutive days --------------------------------------------------

# Quantify consecutive missing days
con_miss <- function(df){
  ex1 <- rle(is.na(df$temp))
  ind1 <- rep(seq_along(ex1$lengths), ex1$lengths)
  s1 <- split(1:nrow(df), ind1)
  res <- do.call(rbind,
                 lapply(s1[ex1$values == TRUE], function(x)
                   data.frame(index_start = min(x), index_end = max(x))
                 )) %>%
    mutate(duration = index_end - index_start + 1) %>%
    group_by(duration) %>%
    summarise(count = n())
  return(res)
}


# Add trends --------------------------------------------------------------

# Function for adding trends to data
add_trend <- function(rate, df = sst_ALL_flat){
  daily_step <- rate/3652.5
  res <- df %>%
    group_by(site, rep) %>%
    mutate(row_index = 1:n(),
           temp = temp + (row_index*daily_step)) %>%
    mutate(index_vals = rate) %>%
    select(-row_index)
  return(res)
}


# Calculate clims, events, and cats ---------------------------------------

# The climatologies themselves are much smaller and easier to handle on their own
# So we want to pull out only the 366 day clims per run
clim_only <- function(df){
  res <- df$climatology %>%
    select(doy, seas:thresh) %>%
    unique() %>%
    mutate_if(is.numeric, round, 3) %>%
    arrange(doy)
  return(res)
}

# Likewise we want to grab only the event values we are interested in
event_only <- function(df){
  res <- df$event %>%
    select(date_start, date_peak, date_end, duration, intensity_mean, intensity_max, intensity_cumulative) %>%
    mutate_if(is.numeric, round, 3)
  return(res)
}

# Calculate climatologies, events, and categories on full length time series
# testers...
# df <- sst_ALL_knockout %>%
  # filter(site == "WA", rep == "1", index_vals == 0.4)
# fix <- "none"
# fix <- "missing"
clim_event_cat_calc <- function(df, fix = "none"){
  res_base <- df %>%
    nest()
  if(fix == "none"){
    res_clim <- res_base %>%
      mutate(clims = map(data, ts2clm,
                         climatologyPeriod = c("1982-01-01", "2011-12-31"))) #%>%
      # select(-data) %>%
      # unnest()
  } else if(fix == "missing"){
    res_clim <- res_base %>%
      mutate(clims = map(data, ts2clm, maxPadLength = 9999, # dummy length
                         climatologyPeriod = c("1982-01-01", "2011-12-31"))) #%>%
      # select(-data) %>%
      # unnest()
  } else{
    stop("Make sure the argument provided to 'fix' is correct.")
  }
  res <- res_clim %>%
    mutate(events = map(clims, detect_event),
           cat = map(events, category),
           clim = map(events, clim_only),
           event = map(events, event_only)) %>%
    select(-data, -clims, -events)
  return(res)
}


# Calculate clim, events, and cats on short time series -------------------

# NB: The shortened time series require their own function for calculating results
# in order to match the length of any given time series
# testers...
# df <- sst_ALL_flat %>%
  # filter(site == "WA", rep == "1")
# year_begin <- 2000
# set_width <- 10
shrinking_results <- function(year_begin, df = sst_ALL_flat, set_width = 5){
  res <- df %>%
    filter(year(t) >= year_begin) %>%
    mutate(index_vals = 2018-year_begin+1) %>%
    group_by(site, rep, index_vals) %>%
    nest() %>%
    mutate(clims = map(data, ts2clm, maxPadLength = 1, windowHalfWidth = set_width,
                       climatologyPeriod = c(paste0(year_begin,"-01-01"), "2018-12-31")),
           events = map(clims, detect_event),
           cat = map(events, category),
           clim = map(events, clim_only),
           event = map(events, event_only)) %>%
    select(site, rep, index_vals, clim, event, cat)
  return(res)
}


# KS tests ----------------------------------------------------------------

# The problem with running KS tests on category data is that it doesn't appear
# to be bothered by differences in sample sizes
# Meaning that a comparison across the count of events matters more on what
# counts are present, rather than the difference in total counts

# Function for running a pairwise KS test against the index_standard
# Lazily I have it run on chosen columns by name given certain tests
# as getting it to do this programmatically was becoming annoying
# testers...
# df_2 <- filter(df_long, index_vals == 0.2)
# df_1 <- df_standard
KS_sub <- function(df_2, df_1){
  if("seas" %in% names(df_1)){
    suppressWarnings( # Suppress warnings about ties
    res <- data.frame(seas = round(ks.test(df_1$seas, df_2$seas)$p.value, 4),
                      thresh = round(ks.test(df_1$thresh, df_2$thresh)$p.value, 4))
    )
  }
  if("intensity_mean" %in% names(df_1)){
    suppressWarnings( # Suppress warnings about ties
      res <- data.frame(duration = round(ks.test(df_1$duration, df_2$duration)$p.value, 4),
                        intensity_mean = round(ks.test(df_1$intensity_mean, df_2$intensity_mean)$p.value, 4),
                        intensity_max = round(ks.test(df_1$intensity_max, df_2$intensity_max)$p.value, 4),
                        intensity_cumulative = round(ks.test(df_1$intensity_cumulative, df_2$intensity_cumulative)$p.value, 4))
    )
  }
  if("category" %in% names(df_1)){
    suppressWarnings( # Suppress warnings about ties
    res <- data.frame(p_moderate = round(ks.test(df_1$p_moderate, df_2$p_moderate)$p.value, 4),
                      p_strong = round(ks.test(df_1$p_strong, df_2$p_strong)$p.value, 4),
                      p_severe = round(ks.test(df_1$p_severe, df_2$p_severe)$p.value, 4),
                      p_extreme = round(ks.test(df_1$p_extreme, df_2$p_extreme)$p.value, 4))
    )
  }
  return(res)
}

# Wrapper function to run all pairwise KS tests
# df <- sst_ALL_clim_event_cat[ ,c(1:4,7)] %>%
  # filter(test == "trended", site == "WA", rep == "1") %>%
  # select(-test, -rep, -site)
# df <- sst_res[ ,c(1:3)] %>%
#   filter(test == "length") %>%
#   select(-test)
# df <- sst_res[ ,c(1:2,4)] %>%
#   filter(test == "missing") %>%
#   select(-test)
KS_p <- function(df){
  # Unnest the clim data
  suppressWarnings( # Suppress warning about different factor levels for cat data
  df_long <- df %>%
    # select(-event, -cat) %>%
    unnest()
  )
  # This determines if the data are from the length test or not
  # It then sets the benchmark value (30 years length or 0% missing/ 0C-dec)
  if(max(df_long$index_vals) > 1){
    index_filter <- 30

  } else {
    index_filter <- 0
  }
  # We then must screen out events only for the last decade for the length and trend tests
  if(max(df_long$index_vals) != 0.5){
    if("intensity_mean" %in% names(df_long)){
      df_long <- df_long %>%
        filter(date_start >= "2009-01-02", date_end <= "2018-12-31")
    }
    if("category" %in% names(df_long)){
      df_long <- df_long %>%
        filter(peak_date >= "2009-01-01", peak_date <= "2018-12-31")
    }
  }
  # Lastly remove events from time series calculated with fewer than 10 years of data
  # to ensire equitable sample sizes
  if(max(df_long$index_vals) > 1){
    df_long <- df_long %>%
      filter(index_vals >= 10)
  }
  # Set the standard results for comparison
  df_standard <- df_long %>%
    filter(index_vals == index_filter)

  # Exit if no events occurred in the last decade of data
    # This is possible in a few icey areas
  if(nrow(df_standard) == 0) return()

  # Run the KS tests
  res <- df_long %>%
    filter(index_vals != index_filter) %>%
    group_by(index_vals) %>%
    nest() %>%
    mutate(KS = map(data, KS_sub, df_1 = df_standard)) %>%
    select(-data) %>%
    unnest()
  return(res)
}


# AOV and Tukey tests -----------------------------------------------------

# Run an ANOVA on each metric of the combined event results and get the p-value
# df <- res
aov_p <- function(df){
  aov_models <- df[ , -grep("index_vals", names(df))] %>%
    map(~ aov(.x ~ df$index_vals)) %>%
    map_dfr(~ broom::tidy(.), .id = 'metric') %>%
    mutate(p.value = round(p.value, 4)) %>%
    filter(term != "Residuals") %>%
    select(metric, p.value)
  return(aov_models)
}

# Run an ANOVA on each metric and then a Tukey test
tukey_p <- function(df){
  tukey_models <- df[ , -grep("index_vals", names(df))] %>%
    map(~ TukeyHSD(aov(.x ~ df$index_vals))) %>%
    map_dfr(~ broom::tidy(.), .id = 'metric') %>%
    mutate(adj.p.value = round(adj.p.value, 4)) %>%
    select(metric, comparison, adj.p.value) %>%
    arrange(metric, adj.p.value)
  return(tukey_models)
}

# Quick wrapper for getting results for ANOVA and Tukey on clims
# df <- sst_ALL_clim_event_cat %>%
  # filter(test == "length", site == "WA", rep == "1") %>%
  # select(-test, -rep, -site)
aov_tukey <- function(df){
  # Extract MHW event data
  prep <- df %>%
    select(index_vals, event) %>%
    unnest()
  # Filter out results not from the most recent decade for equal comparison
  # when looking at the changing lengtths test
  if(max(prep$index_vals) > 1){
    prep <- prep %>%
      filter(index_vals >= 10) %>%
      # 2009-01-02 is intentional
      filter(date_start >= "2009-01-02", date_end <= "2018-12-31")
  }
  # Calculate the ANOVA and Tukey p values
  res <- prep %>%
    select(-date_start, -date_peak, -date_end) %>%
    mutate(index_vals = as.factor(as.character(index_vals))) %>%
    nest() %>%
    mutate(aov = map(data, aov_p),
           tukey = map(data, tukey_p)) %>%
    select(-data)
  return(res)
}



# Fisher test for category counting ---------------------------------------


# Fisher pairwise function
# Takes one rep of missing data and compares it against the complete data
# testers...
# df <- unnest(slice(select(res, data), 1))
# df_comp <- df_standard
fisher_pair <- function(df, df_comp){
  # Prep data
  df_joint <- table(rbind(df_comp, df))
  if(ncol(df_joint) == 1){
    df_joint <- as.table(cbind(df_joint, 'III & IV' = c(0,0)))
  }
  # Run tests
  # res <- round(fisher.test(df_joint)$p.value, 4)
  res <- fisher.test(df_joint)
  # res_broom <- broom::augment(res) #%>%
    # mutate_if(is.numeric, round, 4)
  return(res)
}

# Wrapper to get data ready for pair-wise chi-squared tests
# testers...
# df <- sst_ALL_clim_event_cat %>%
  # filter(test == "length", site == "WA", rep == "1") %>%
  # select(-test, -rep, -site)
fisher_test <- function(df){
  suppressWarnings( # Suppress warning about category levels not all being present
    df_long <- df %>%
      select(index_vals, cat) %>%
      unnest() %>%
      # select(index_vals, category) %>%
      mutate(category = ifelse(category %in% c("I Moderate","II Strong"), "I & II", category),
             category = ifelse(category %in% c("III Severe","IV Extreme"), "III & IV", category),
             category = as.character(category))
  )
  # This determines if the data are from the length test or not
  # It then sets the benchmark value (30 years length or 0% missing/ 0C-dec)
  if(max(df_long$index_vals) > 1){
    index_filter <- 30
    df_long <- df_long %>%
      filter(index_vals >= 10,
             peak_date >= "2009-01-01", peak_date <= "2018-12-31")
  } else {
    index_filter <- 0
  }
  res_table <- table(select(df_long, index_vals, category))
  if(ncol(res_table) == 1){
    res_table <- as.table(cbind(res_table, 'III & IV' = rep(0, nrow(res_table))))
  }
  res <- rcompanion::pairwiseNominalIndependence(res_table, fisher = TRUE, gtest = FALSE,
                                     chisq = FALSE, method = "bonferroni",
                                     digits = 3)
  res_tidy <- separate(res, col = "Comparison", into = c("group_1", "group_2"), sep = " : ") %>%
    filter(group_1 == index_filter | group_2 == index_filter) %>%
    mutate(index_vals = ifelse(group_1 == index_filter, group_2, group_1),
           p.adj.Fisher = p.Fisher* n(),
           p.adj.Fisher = ifelse(p.adj.Fisher > 1, 1, p.adj.Fisher)) %>%
    select(index_vals, p.Fisher, p.adj.Fisher)
  return(res_tidy)
}


# Model linear relationships ----------------------------------------------

# Wrapper function for nesting
# testers...
# df <- slice(sst_ALL_R2_long, 1) %>%
  # unnest()
lm_p_R2 <- function(df){
  res <- lm(p.value.mean~index_vals, data = df)
  res_broom <- broom::glance(res) %>%
    mutate_if(is.numeric, round, 4)
  return(res_broom)
}



# Convenience functions ---------------------------------------------------

# Wrapper function to melt KS results
KS_long <- function(df){
  res <- df %>%
    gather(key = "metric", value = "p.value", -test, -site, -rep, -index_vals) %>%
    group_by(test, site, metric, index_vals) %>%
    summarise(p.value.mean = mean(p.value))
  return(res)
}


# Event effect functions --------------------------------------------------

# The effect that the three tests have the climatologies
effect_clim_func <- function(df, choice_rep = "1"){
  res <- df %>%
    filter(rep == choice_rep) %>%
    select(-event, -cat) %>%
    unnest(clim) %>%
    select(-doy, -rep) %>%
    gather(key = "metric", value = "val", -site, -test, -index_vals) %>%
    group_by(site, test, index_vals, metric) %>%
    summarise_if(is.numeric,
                 .funs = c("min", "median", "mean", "max")) %>%
    # group_by(site, test) %>%
    mutate_if(is.numeric, round, 3) %>%
    ungroup()
  return(res)
}

# The effect that the three tests have on the focus event
# df <- sst_res
# date_guide <- focus_event
# choice_rep <- "1"
effect_event_func <- function(df, date_guide, choice_rep = "1"){
  res <- df %>%
    filter(rep == choice_rep) %>%
    select(-clim, -cat) %>%
    unnest(event) %>%
    left_join(date_guide, by = "site") %>%
    filter(date_peak >= date_start_control,
           date_peak <= date_end_control) %>%
    group_by(site, test, index_vals) %>%
    #
    # Not sure what to do with this information
    # mutate(date_start_change = date_start_control - date_start,
    #        date_peak_change = date_peak_control - date_peak,
    #        date_end_change = date_end_control - date_end) %>%
    #
    summarise(count = n(),
              duration = sum(duration),
              intensity_mean = mean(intensity_mean),
              intensity_max = max(intensity_max),
              intensity_cumulative = sum(intensity_cumulative)) %>%
    mutate_if(is.numeric, round, 2) %>%
    gather(key = "metric", value = "val", -site, -test, -index_vals) %>%
    ungroup()
  return(res)
}

# The effect that the three tests have on the categories
# df <- sst_res
# date_guide <- focus_event
# choice_rep <- "1"
effect_cat_func <- function(df, date_guide, choice_rep = "1"){
  res <- df %>%
    filter(rep == choice_rep) %>%
    select(-clim, -event) %>%
    unnest(cat) %>%
    left_join(date_guide, by = "site") %>%
    filter(peak_date >= date_start_control,
           peak_date <= date_end_control) %>%
    group_by(site, test, index_vals) %>%
    summarise(count = n(),
              duration = sum(duration),
              i_max = max(i_max),
              p_moderate = mean(p_moderate),
              p_strong = mean(p_strong),
              p_severe = mean(p_severe),
              p_extreme = mean(p_extreme)) %>%
    mutate_if(is.numeric, round, 2) %>%
    gather(key = "metric", value = "val", -site, -test, -index_vals) %>%
    ungroup()
  return(res)
}


# Global functions --------------------------------------------------------

# The following wrapper functions use the above functions,
# but allow them to be used on the NOAA OISST NetCDF file structure

# It's not reasonable to run a decadal trend test as well because
# this value isbeing added in a completely controlled way.
# Meaning that controlling for it is the same as simply not adding it
# in the first place.
# To that end we want to provide a post-hoc correction.
# Because we have seen that max. intensity is a function of the
# decadal trend and the peak date of the event we should be able to
# correct for intensities by subtracting the decadal trend multiplied
# by the peak date as it relates to the length of the time series.
# For example, if the decadal trend is 0.3C, and the peak date of an event
# is in the 25th year of a 30 year time series (0.8333 of the length),
# then the the impact of the decadal trend on the max. intensity should be
# 0.3*0.83 = 0.25C
# So one would subtract this value in order to correct the results to
# match a de-trended time series.
# The problem then becomes, why would anyone actually want to do this?
# No. Rather we must provide a post-hoc fix for the potential impact of
# a decadal trend on MHWs in conjunction with shot time series.
# And then on top of that add in the correction for missing data.
# Mercifully the correction for missing data is very simple and should
# play nice with the other issues.
# This is because the linear interpolation of NA gaps will provide
# data that are matching the temperatures around them so the problem of
# having more missing data on one end of a time series over the other
# should not be pronounced.

# Function for knocking out data but maintaing the time series consistency
# df <- sst
# prop <- 0.1
random_knockout_global <- function(prop, df){
  # NB: Don't allow samppling of first and last value to ensure
  # all time series are the same length
  ts_length <- nrow(df)
  miss_index <- sample(seq(2, ts_length-1, 1), ts_length*prop, replace = F)
  # df$x1[sample(nrow(df),250)] <- NA
  res <- df %>%
    # group_by(site, rep) %>%
    mutate(row_index = 1:n(),
           temp = replace(temp, which(row_index %in% miss_index), NA)) %>%
    mutate(index_vals = prop) %>%
    select(-row_index)
  return(res)
}

# Function for calculating shrinking time series from the global data
shrinking_results_global <- function(year_begin, df){
  res <- df %>%
    filter(year(t) >= year_begin) %>%
    mutate(index_vals = 2018-year_begin+1) %>%
    group_by(index_vals) %>%
    nest() %>%
    mutate(clims = map(data, ts2clm,
                       climatologyPeriod = c(paste0(year_begin,"-01-01"), "2018-12-31")),
           events = map(clims, detect_event),
           cat = map(events, category),
           clim = map(events, clim_only),
           event = map(events, event_only)) %>%
    select(index_vals, clim, event, cat)
  return(res)
}

# The function that runs all of the tests on a single pixel
# lat_step <- 78
global_analysis_sub <- function(lat_step, nc_file){

  nc <- nc_open(nc_file)
  lat <- nc$dim$lat$vals[lat_step]
  lon <- nc$dim$lon$vals[1]

  # Extract raw SST
  sst_raw <- ncvar_get(nc, "sst", start = c(lat_step,1,1), count = c(1,-1,-1))
  dimnames(sst_raw) <- list(t = nc$dim$time$vals)
  nc_close(nc)

  # Prep SST for further use
  sst <- as.data.frame(reshape2::melt(sst_raw, value.name = "temp"), row.names = NULL) %>%
    mutate(t = as.Date(t, origin = "1970-01-01")) %>%
    na.omit() %>%
    filter(t <= "2018-12-31") %>%
      mutate(site = as.character(lat), rep = "1")
  if(nrow(sst) == 0) return()

  # Some icey bits don't start on 1982-01-01, or go until 2018-12-31, which is problematic
    # So for now we are simply removing them
  if(min(sst$t) > as.Date("1982-01-01"))  return()
  if(max(sst$t) < as.Date("2018-12-31"))  return()

  # Calculate the secular trend
  dec_trend <- round(as.numeric(broom::tidy(lm(temp ~ t, sst))[2,2]*3652.5), 3)

  # Remove random data
  set.seed(666)
  sst_knockout <- plyr::ldply(.data = c(0.10, 0.25, 0.50, 0.75, 0.90),
                              .fun = random_knockout_global, df = sst)

  # Manually force certain key dates to be present
  sst_knockout$temp[sst_knockout$t == "1982-01-01"] <- sst$temp[sst$t == "1982-01-01"]
  sst_knockout$temp[sst_knockout$t == "1989-01-01"] <- sst$temp[sst$t == "1989-01-01"]
  sst_knockout$temp[sst_knockout$t == "1999-01-01"] <- sst$temp[sst$t == "1999-01-01"]
  sst_knockout$temp[sst_knockout$t == "2009-01-01"] <- sst$temp[sst$t == "2009-01-01"]
  sst_knockout$temp[sst_knockout$t == "2018-12-31"] <- sst$temp[sst$t == "2018-12-31"]

  # # Make base calculations
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

  # Run KS tests
  sst_KS_clim <- plyr::ddply(.data = sst_res[ ,c(1:3)], .variables = c("test"), .fun = KS_p)
  sst_KS_event <- plyr::ddply(.data = sst_res[ ,c(1:2,4)], .variables = c("test"), .fun = KS_p)
  sst_KS_cat <- plyr::ddply(.data = sst_res[ ,c(1:2,5)], .variables = c("test"), .fun = KS_p)

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
              KS_clim = sst_KS_clim,
                 KS_event = sst_KS_event,
                 KS_cat = sst_KS_cat,
                 effect_clim = effect_clim_res,
                 effect_event = effect_event_res,
                 effect_cat = effect_cat_res)
  return(res)
}

# One function to rule them all
# nc_file <- OISST_slice
global_analysis <- function(nc_file){
  # Function that performs all of the analyses in one llply call
    # Run through the latitudes as individual time series
  # test_run <- global_analysis_sub(100, nc_file)
  # system.time(
    suppressWarnings( # Ignore the coercing factor to character messages
      res_pixel <- plyr::llply(.data = seq(1, 720),
                               .fun = global_analysis_sub,
                               nc_file = nc_file, .parallel = T)
    )
  # )  # ~7 seconds for one, ~78 seconds for 720
  return(res_pixel)
}


# Figure convenience functions --------------------------------------------

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

# The code two create figure 2, which shows the p-values for the KS tests
fig_2_plot <- function(df){
  fig_2 <- ggplot(df, aes(x = index_vals, y = p.value.mean, colour = metric)) +
    geom_line() +
    geom_point() +
    geom_point(data = filter(df, p.value.mean <= 0.05), shape = 15, colour = "red", size = 2) +
    geom_hline(yintercept = 0.05, colour = "red", linetype = "dashed") +
    scale_colour_manual(name = "Metric",
                        values = c("skyblue", "navy",
                                   "springgreen", "forestgreen",
                                   "#ffc866", "#ff6900", "#9e0000", "#2d0000"),
                        labels = c("Seasonal climatology", "Threshold climatology",
                                   "Duration (days)", "Max. intensity (°C)",
                                   "Prop. moderate", "Prop. strong", "Prop. severe", "Prop. extreme")) +
    # geom_errorbarh(aes(xmin = year_long, xmax = year_short)) +
    guides(colour = guide_legend(override.aes = list(shape = 15, linetype = NA, size = 3))) +
    labs(y = "Mean p-value (n = 100)", x = NULL) +
    facet_grid(site~test, scales = "free_x", switch = "x") +
    theme(legend.position = "bottom")
}
