# code/functions.R
# This script contains all of the functions used throughout 'code/workflow.R'


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(broom)
library(heatwaveR, lib.loc = "~/R-packages/")
# cat(paste0("heatwaveR version = ",packageDescription("heatwaveR")$Version))
library(lubridate) # This is intentionally activated after data.table
library(fasttime)
library(ggpubr)
library(boot)
library(FNN)
library(mgcv)
library(doMC); doMC::registerDoMC(cores = 50)
library(ggridges)


# Meta-data ---------------------------------------------------------------

MHW_colours <- c(
  "I Moderate" = "#ffc866",
  "II Strong" = "#ff6900",
  "III Severe" = "#9e0000",
  "IV Extreme" = "#2d0000"
)

OISST_files <- dir(path = "../data/OISST", full.names = T)


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
  # NB: Don't allow samppling of first and last value to ensure
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
    arrange(doy)
  return(res)
}

# Likewise we want to grab only the event values we are interested in
event_only <- function(df){
  res <- df$event %>%
    select(duration, intensity_mean, intensity_max, intensity_cumulative)
  return(res)
}

# Calculate climatologies, events, and categories on full length time series
# Set `maxPadLength = 1` so that ts2clm() will fill in the missing leap year days
# caused by the re-sampling step above
# testers...
# df <- sst_ALL_knockout %>%
  # filter(site == "WA", rep == "1", index_vals == 0)
clim_event_cat_calc <- function(df){
  res <- df %>%
    nest() %>%
    mutate(clims = map(data, ts2clm, maxPadLength = 1,
                       climatologyPeriod = c("1982-01-01", "2011-12-31")),
           events = map(clims, detect_event),
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
shrinking_results <- function(year_begin, df = sst_ALL_flat){
  res <- df %>%
    filter(year(t) >= year_begin) %>%
    mutate(index_vals = 2018-year_begin+1) %>%
    group_by(site, rep, index_vals) %>%
    nest() %>%
    mutate(clims = map(data, ts2clm, maxPadLength = 1,
                       climatologyPeriod = c(paste0(year_begin,"-01-01"), "2018-12-31")),
           events = map(clims, detect_event),
           cat = map(events, category),
           clim = map(events, clim_only),
           event = map(events, event_only)) %>%
    select(-data, -clims, -events)
  return(res)
}


# KS tests ----------------------------------------------------------------

# Function for running a pairwise KS test against the index_standard
KS_sub <- function(df_2, df_1){
  suppressWarnings( # Suppress warnings about perfect matches
    res <- data.frame(seas = round(ks.test(df_1$seas, df_2$seas)$p.value, 4),
                      thresh = round(ks.test(df_1$thresh, df_2$thresh)$p.value, 4))
  )
  return(res)
}

# # Wrapper function to run all pairwise KS tests
# # testers...
# df <- sst_ALL_flat %>%
#   filter(site == "WA", rep == "1")
# df <- sst_ALL_clim_only %>%
#   filter(site == "WA", test == "trended", rep == "1")
KS_p <- function(df){
  # Unnest the clim data
  df_long <- df %>%
    select(-event, -cat) %>%
    unnest()
  # This determines if the data are from the length test or not
  # It then sets the benchmark value (30 years length or 0% missing/ 0C-dec)
  if(max(df_long$index_vals) > 1){
    index_filter <- 30
  } else {
    index_filter <- 0
  }
  df_standard <- df_long %>%
    filter(index_vals == index_filter)
  res <- df %>%
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
  aov_tukey <- df[ , -grep("index_vals", names(df))] %>%
    map(~ TukeyHSD(aov(.x ~ df$index_vals))) %>%
    map_dfr(~ broom::tidy(.), .id = 'metric') %>%
    mutate(p.value = round(adj.p.value, 4)) %>%
    select(metric, comparison, adj.p.value) %>%
    arrange(metric, adj.p.value)
  return(aov_tukey)
}

# Quick wrapper for getting results for ANOVA and Tukey on clims
# df <- sst_ALL_clim_only %>%
#   filter(test == "length" rep == "1", site == "WA") #%>%
  # select(-rep, -site)
aov_tukey <- function(df){
  # Extract MHW event data
  prep <- df %>%
    select(-cats) %>%
    unnest(events) %>%
    filter(row_number() %% 2 == 0) %>%
    unnest(events) %>%
    select(index_vals, duration, intensity_mean, intensity_max, intensity_cumulative) %>%
    mutate(index_vals = as.factor(index_vals))
  # Filter out results not from the most recent decade for equal comparison
  # when looking at the changing lengtths test
  if(df$test[1] == "length")
    prep <- prep %>%
      filter(index_vals >= 10) %>%
      filter(date_start >= "2009-01-02", date_end <= "2018-12-31")
  # Calculate the ANVA and Tukey p values
  res <- prep %>%
    nest() %>%
    mutate(aov = map(data, aov_p),
           tukey = map(data, tukey_p)) %>%
    select(-data)
  return(res)
}


# Chi-squared category counting -------------------------------------------

# Prep the data for better chi-squared use
# suppressWarnings( # Suppress warning about category levels not all being present
# sst_ALL_cat_prep <- sst_ALL_clim_event_cat %>%
#   select(-events) %>%
#   unnest() %>%
#   select(test:index_vals, category) %>%
#   mutate(category = ifelse(category == "III Severe", "III & IV", category),
#          category = ifelse(category == "IV Extreme", "III & IV", category))
# )

# Extract the true climatologies
# sst_ALL_miss_cat_only_0 <- sst_ALL_miss_cat_prep %>%
#   filter(miss == 0, rep == "1") %>%
#   select(site, miss, category)

# chi-squared pairwise function
# Takes one rep of missing data and compares it against the complete data
# df <- sst_ALL_miss_cat_prep %>%
# filter(site == "WA", rep == "20", miss == 0.1)
# df <- unnest(slice(sst_ALL_miss_cat_chi,1)) %>%
# select(-site2, -rep, -miss2)
# chi_pair <- function(df){
#   # Prep data
#   df_comp <- sst_ALL_miss_cat_only_0 %>%
#     filter(site == df$site[1])
#   df_joint <- rbind(df, df_comp)
#   # Run tests
#   res <- chisq.test(table(df_joint$index_vals, df_joint$category), correct = FALSE)
#   res_broom <- broom::augment(res) %>%
#     mutate(p.value = broom::tidy(res)$p.value)
#   return(res_broom)
# }

# The pairwise chai-squared results
# suppressWarnings( # Suppress poor match warnings
#   sst_ALL_miss_cat_chi <- sst_ALL_miss_cat_prep %>%
#     filter(miss != 0) %>%
#     # mutate(miss = as.factor(miss)) %>%
#     mutate(site2 = site,
#            miss2 = miss) %>%
#     group_by(site2, miss2, rep) %>%
#     nest() %>%
#     mutate(chi = map(data, chi_pair)) %>%
#     select(-data) %>%
#     unnest() %>%
#     dplyr::rename(site = site2, miss = miss2,
#                   miss_comp = Var1, category = Var2)
# )


# Global functions --------------------------------------------------------

# The following wrapper functions use the above functions,
# but allow them to be used on the NOAA OISST NetCDF file structure
