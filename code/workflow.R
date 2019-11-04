# code/workflow.R
# This script shows the high-level overview of the full analysis for the manuscript
# Much of this code invokes functions that may be found in 'code/functions.R'
# Note that line 10 will likely not run correctly, see 'code/functions.R'
# for an explanation why and a quick fix for one's local machine


# Libraries ---------------------------------------------------------------

# Load all meta-data and functions
source("code/functions.R")

# Remove scientific notation from results
options(scipen=999)


# Reference analysis ------------------------------------------------------

# Combine the three reference time series, run analysis, and save
sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
                 mutate(sst_NW_Atl, site = "NW_Atl"),
                 mutate(sst_Med, site = "Med"))
system.time(
  sst_ALL_results <- plyr::ddply(sst_ALL, c("site"), single_analysis, .parallel = T,
                                 full_seq = T, clim_metric = F, count_miss = T, windows = T)
) # 65 seconds
saveRDS(sst_ALL_results, "data/sst_ALL_results.Rda")


# Random analysis ---------------------------------------------------------

# Calculate the full analysis on 1000 random pixels
doMC::registerDoMC(cores = 25) # 50 appears to use too much RAM
set.seed(666)
system.time(
  random_results <- plyr::ldply(1:1000, random_analysis, .parallel = T)
) # ~54 minutes
saveRDS(random_results, "data/random_results_1000.Rda")


# Global analysis ---------------------------------------------------------

# Wrapper function with a chatty output as it chugs along
global_analysis_single <- function(file_sub, par_op = F){
  OISST_slice <- OISST_files[file_sub]
  lon_row_pad <- str_pad(file_sub, width = 4, pad = "0", side = "left")
  print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))
  slice_res <- global_analysis(OISST_slice, par_op = par_op)
  saveRDS(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rda"))
  print(paste0("Finished run on step ",lon_row_pad," at ",Sys.time()))
  rm(slice_res); gc()
}

# Run ALL
doParallel::registerDoParallel(cores = 50)
plyr::l_ply(1:1440, global_analysis_single, .parallel = T) # This took 37.5 hours to run


# Global trends -----------------------------------------------------------

# Calculate the simple linear trends for the different tests at each pixel
doParallel::registerDoParallel(cores = 50)
system.time(
  global_var_trend <- plyr::ldply(dir("data/global", full.names = T, pattern = "slice"),
                                  .fun = var_trend, .parallel = T)
) # 102 seconds for one lon slice, ~90 minutes for all
saveRDS(global_var_trend, "data/global_var_trend.Rda")


# Figures -----------------------------------------------------------------

# All of the figures that use the above results are made in 'code/figures.R'


# Results -----------------------------------------------------------------

# Code that generates numeric results referred to in the text outside of figures

# Load random (non-global) results
system.time(
  random_results <- readRDS("data/random_results_1000.Rda") %>%
    unite("site", c(lon, lat))
) # 68 seconds, 15 seconds without the "site" column

# The choice variables for focussing on
var_choice <- data.frame(var = c("count", "duration", "intensity_max",
                                 "focus_count", "focus_duration", "focus_intensity_max"),
                         id = c("n_perc", "sum_perc", "mean_perc",
                                "mean_perc", "sum_perc", "mean_perc"),
                         stringsAsFactors = F)

# Upper and lower quantile ranges for sub-optimal tests
quant_subopt <- random_results %>%
  right_join(var_choice, by = c("var", "id")) %>%
  group_by(test, index_vals, var, id) %>%
  summarise(lower = round(quantile(val, 0.05), 2),
            upper = round(quantile(val, 0.95), 2)) %>%
  ungroup()

# The quantiles by test
quant_length <- filter(quant_subopt, test == "length")
quant_miss <- filter(quant_subopt, test == "missing")
quant_interp <- filter(quant_subopt, test == "interp")
quant_trend <- filter(quant_subopt, test == "trend")


# Best practices ----------------------------------------------------------

# Code that produces the tables and supports the statements made in the Best Practices section
# Create table showing the rates of change int the results
# Find the 5th, 25th, 50th,75th, and 95th quantile at each step
# Fit linear models to those and provide those trends + R2 values
# Also show where the inflection points may be where the trends change

# Load the random 1000 data
system.time(
  random_results <- readRDS("data/random_results_1000.Rda") %>%
    unite("site", c(lon, lat))
) # 68 seconds, 15 seconds without the "site" column

# The choice variables for focussing on
var_choice <- data.frame(var = c("count", "duration", "intensity_max", "focus_count", "focus_duration", "focus_intensity_max"),
                         id = c("n_perc", "sum_perc", "mean_perc", "mean_perc", "sum_perc", "mean_perc"),
                         stringsAsFactors = F)

# Calculate the full range of quantiles
random_quant <- random_results %>%
  right_join(var_choice, by = c("var", "id")) %>%
  mutate(test = as.character(test)) %>%
  filter(test %in% c("length", "missing", "interp", "trend")) %>%
  group_by(test, index_vals, var, id) %>%
  summarise(q05 = quantile(val, 0.05),
            q25 = quantile(val, 0.25),
            q50 = quantile(val, 0.50),
            q75 = quantile(val, 0.75),
            q95 = quantile(val, 0.95),
            iqr50 = q75-q25,
            iqr90 = q95-q05) %>%
  ungroup()

# Calculate the table for the Best Practices section
slope_final <- random_quant %>%
  gather(key = "stat", value = "val", -c(test:id)) %>%
  mutate(test2 = test,
         var2 = var) %>%
  group_by(test2, var2, id, stat) %>%
  group_modify(~trend_correct(.x)) %>%
  dplyr::rename(Test = test2,
                Variable = var2,
                Range = range) %>%
  ungroup() %>%
  select(-id) %>%
  select(-intercept, -p) %>%
  filter(!(stat %in% c("iqr50", "iqr90"))) %>%
  group_by(Test, Variable, Range, stat) %>%
  mutate(slope = sprintf(slope, fmt = "%0.2f", how = "replace"),
         R2 = sprintf(R2, fmt = "%0.2f", how = "replace")) %>%
  mutate(slope = paste0(slope,"%"),
         R2 = paste0("(",R2,")")) %>%
  unite(slope, R2, col = "slope_R2", sep = " ") %>%
  spread(stat, slope_R2) %>%
  select(Test, Variable, Range, q05, q25, q50, q75, q95) %>%
  ungroup() %>%
  mutate(Test = factor(Test, levels = c("length", "missing", "interp", "trend"))) %>%
  arrange(Test)

# The table for the average MHW results
best_table_average <- slope_final %>%
  filter(!grepl("focus", Variable)) %>%
  mutate(Variable = str_replace(Variable, "intensity_max", "max. intensity"))
saveRDS(best_table_average, "data/best_table_average.Rda")
write_csv(best_table_average, "data/table_1.csv")

# The table for the focal MHW results
best_table_focus <- slope_final %>%
  filter(grepl("focus", Variable)) %>%
  mutate(Variable = str_remove(Variable, "focus_"),
         Variable = str_replace(Variable, "intensity_max", "max. intensity"))
saveRDS(best_table_focus, "data/best_table_focus.Rda")
write_csv(best_table_focus, "data/table_2.csv")


# Discussion --------------------------------------------------------------

# No one here but us chickens...


# Supplementary 1 ---------------------------------------------------------

# The effect of the sub-optimal tests on seas/thresh
# This figure is created in 'code/figures.R'


# Supplementary 2 ---------------------------------------------------------

# The difference between the proper 30 year base period and all other 30 year base periods

# Combine the three reference time series, run analysis, and save
sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
                 mutate(sst_NW_Atl, site = "NW_Atl"),
                 mutate(sst_Med, site = "Med"))
system.time(
  sst_ALL_results <- plyr::ddply(sst_ALL, c("site"), base_period_analysis, .parallel = T)
) # 3 seconds
saveRDS(sst_ALL_results, "data/sst_ALL_bp_results.Rda")

# Calculate the base period analysis on 1000 random pixels
doParallel::registerDoParallel(cores = 50)
set.seed(666)
system.time(
  random_results <- plyr::ldply(1:1000, random_analysis, .parallel = T, base_period = T)
) # 312 seconds
saveRDS(random_results, "data/random_bp_results_1000.Rda")

# This figure is created in 'code/figures.R'


# Supplementary 3 ---------------------------------------------------------

# The global patterns in missing data are unremarkable and generally consistent across the oceans
# This figure is created in 'code/figures.R'


# Supplementary 4 ---------------------------------------------------------

# The global patterns in added decadal trends generally show that MHW metrics increase
# This figure is created in 'code/figures.R'


# Supplementary 5 ---------------------------------------------------------

# The effect of widening window half-widths on the focal MHWs
# This figure is created in 'code/figures.R'

