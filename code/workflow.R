# workflow.R


# Libraries ---------------------------------------------------------------

source("code/functions.R")


# Base data ---------------------------------------------------------------

# Combine the three reference time series into one file
sst_ALL <- rbind(sst_Med, sst_NW_Atl, sst_WA) %>%
  mutate(site = as.factor(rep(c("Med", "NW_Atl", "WA"), each = nrow(sst_WA))))

# Save and clear
save(sst_ALL, file = "data/sst_ALL.Rdata")
# rm(sst_ALL); gc()


# Re-sample the data ------------------------------------------------------

# Re-sample the data 100 times
doMC::registerDoMC(cores = 50)
set.seed(666)
sst_ALL_repl <- plyr::ldply(2:100, sample_37, .parallel = T)

# Add the original data as rep = "1"
sst_ALL_1 <- sst_ALL %>%
  mutate(year_orig = year(t),
         rep = "1")
sst_ALL_repl <- rbind(sst_ALL_1, sst_ALL_repl) %>%
  mutate(rep = factor(rep))

# Save and clear
save(sst_ALL_repl, file = "data/sst_ALL_repl.Rdata")
# rm(sst_ALL_1, sst_ALL_repl); gc()


# De-trend the data -------------------------------------------------------

# Create de-trended anomaly time series from all re-samples
doMC::registerDoMC(cores = 50)
sst_ALL_flat <- plyr::ddply(sst_ALL_repl, c("site", "rep"), detrend)

# Save and clear
save(sst_ALL_flat, file = "data/sst_ALL_flat.Rdata")
# rm(sst_ALL_flat); gc()


# Knockout random days ----------------------------------------------------

# Randomly knockout 0 - 50% of each of the 100 re-samples
doMC::registerDoMC(cores = 26)
set.seed(666)
sst_ALL_knockout <- plyr::ldply(seq(0.00, 0.50, 0.01), random_knockout, .parallel = T)

# Save and clear
save(sst_ALL_knockout, file = "data/sst_ALL_knockout.Rdata")
# rm(sst_ALL_knockout); gc()


# Count consecutive missing days ------------------------------------------

# Calculate the consecutive missing days
doMC::registerDoMC(cores = 50)
sst_ALL_consec <- plyr::ddply(filter(sst_ALL_knockout, index_vals != 0),
                           c("site", "rep", "index_vals"), con_miss, .parallel = T)

# Save and clear
save(sst_ALL_consec, file = "data/sst_ALL_consec.Rdata")
# rm(sst_ALL_consec); gc()


# Add trends --------------------------------------------------------------

# Add the trends
doMC::registerDoMC(cores = 50)
sst_ALL_add_trend <- plyr::ldply(seq(0.00, 0.30, 0.01), add_trend, .parallel = T)

# Save and clear
save(sst_ALL_add_trend, file = "data/sst_ALL_add_trend.Rdata")
# rm(sst_ALL_add_trend); gc()


# Calculate clims, events, and categories ---------------------------------

# NB: These are to large to be run as one command
  # They must be manually selected and run
  # Wait for each to finish before running the next

# Begin fresh here if desired
# load("data/sst_ALL_flat.Rdata")
# unique(sst_ALL_flat$rep)
# load("data/sst_ALL_knockout.Rdata")
# unique(sst_ALL_knockout$index_vals)
# load("data/sst_ALL_add_trend.Rdata")
# unique(sst_ALL_add_trend$index_vals)

## Duration
# NB: The time series duration results need to be calculated with a different
# function because the climatology period needs to adjust to the length
# of the shortened time series
# Run this preferably with 35 cores for speed, RAM allowing
doMC::registerDoMC(cores = 35)
sst_ALL_length <- plyr::ldply(1982:2016, shrinking_results, .parallel = T) %>%
  mutate(test = as.factor("length"))
# test that it ran correctly
# names(sst_ALL_length)
# unnest(slice(sst_ALL_length, 1), clim)
# unnest(slice(sst_ALL_length, 1), event)
# unnest(slice(sst_ALL_length, 1), cat)
# save(sst_ALL_length, file = "data/sst_ALL_length.Rdata")

## Missing data
doMC::registerDoMC(cores = 50)
sst_ALL_missing <- plyr::ddply(sst_ALL_knockout, c("site", "rep", "index_vals"),
                               clim_event_cat_calc, .parallel = T) %>%
  mutate(test = as.factor("missing"))
# save(sst_ALL_missing, file = "data/sst_ALL_length.Rdata")

## Trended data
doMC::registerDoMC(cores = 50)
sst_ALL_trended <- plyr::ddply(sst_ALL_add_trend, c("site", "rep", "index_vals"),
                             clim_event_cat_calc, .parallel = T) %>%
  mutate(test = as.factor("trended"))
# save(sst_ALL_trended, file = "data/sst_ALL_trended.Rdata")

## Combine all results
sst_ALL_clim_event_cat <- rbind(sst_ALL_length, sst_ALL_missing, sst_ALL_trended) %>%
  select(test, site, rep, index_vals, clim, event, cat)

## Save and clear
save(sst_ALL_clim_event_cat, file = "data/sst_ALL_clim_event_cat.Rdata")
# rm(sst_ALL_length, sst_ALL_miss, sst_ALL_trend, sst_ALL_clim_event_cat); gc()


# Climatology results -----------------------------------------------------

# Kolmogorov-Smirnov tests for similarity of climatologies
doMC::registerDoMC(cores = 50)
sst_ALL_KS_clim <- plyr::ddply(sst_ALL_clim_event_cat[ ,c(1:5)],
                               c("test", "site", "rep"), KS_p, .parallel = T)

# Save and clear
save(sst_ALL_KS_clim, file = "data/sst_ALL_KS_clim.Rdata")
# rm(sst_ALL_KS_clim); gc()


# Event metric results ----------------------------------------------------

# Kolmogorov-Smirnov tests for similarity of climatologies
doMC::registerDoMC(cores = 50)
sst_ALL_KS_event <- plyr::ddply(sst_ALL_clim_event_cat[ ,c(1:4,6)],
                               c("test", "site", "rep"), KS_p, .parallel = T)

# Save and clear
save(sst_ALL_KS_event, file = "data/sst_ALL_KS_event.Rdata")
# rm(sst_ALL_KS_event); gc()

# AOV and Tukey
doMC::registerDoMC(cores = 50)
sst_ALL_aov_tukey <- plyr::ddply(sst_ALL_clim_event_cat, c("test", "site", "rep"),
                                 aov_tukey, .parallel = T)

# test that it ran correctly
# names(sst_ALL_aov_tukey)
# unnest(slice(sst_ALL_aov_tukey, 1), aov)
# unnest(slice(sst_ALL_aov_tukey, 1), tukey)

# Save and clear
save(sst_ALL_aov_tukey, file = "data/sst_ALL_aov_tukey.Rdata")
# rm(sst_ALL_aov_tukey); gc()


# Category count results --------------------------------------------------

# Fisher
doMC::registerDoMC(cores = 50)
sst_ALL_fisher <- plyr::ddply(sst_ALL_clim_event_cat,
                              c("test", "site", "rep"), fisher_test, .parallel = T)

# Save and clear
save(sst_ALL_fisher, file = "data/sst_ALL_fisher.Rdata")
# rm(sst_ALL_fisher); gc()

# Post-hoc


# Linear relationships ----------------------------------------------------




# Global ------------------------------------------------------------------

# And now we run the above workflow on the global data
