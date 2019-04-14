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

## Length data
# NB: The time series length results need to be calculated with a different
# function because the climatology period needs to adjust to the length
# of the shortened time series
# Run this preferably with 35 cores for speed, RAM allowing
doMC::registerDoMC(cores = 35)
system.time(
sst_ALL_length <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results, .parallel = T)
)  # ~92 seconds
# NB: Adding the test column in a pipe seems to make ldply sad...
sst_ALL_length$test <- as.factor("length")
# test that it ran correctly
# names(sst_ALL_length)
# unnest(slice(sst_ALL_length, 1), clim)
# unnest(slice(sst_ALL_length, 1), event)
# unnest(slice(sst_ALL_length, 1), cat)
# save(sst_ALL_length, file = "data/sst_ALL_length.Rdata")
# load("data/sst_ALL_length.Rdata")

## Missing data
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_missing <- plyr::ddply(sst_ALL_knockout, c("site", "rep", "index_vals"),
                               clim_event_cat_calc, .parallel = T) %>%
  mutate(test = as.factor("missing"))
) # ~235 seconds
save(sst_ALL_missing, file = "data/sst_ALL_missing.Rdata")
# load("data/sst_ALL_missing.Rdata")

## Trended data
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_trended <- plyr::ddply(sst_ALL_add_trend, c("site", "rep", "index_vals"),
                               clim_event_cat_calc, .parallel = T) %>%
  mutate(test = as.factor("trended"))
) # 140 seconds
save(sst_ALL_trended, file = "data/sst_ALL_trended.Rdata")
# load("data/sst_ALL_trended.Rdata")

## Combine all results
sst_ALL_clim_event_cat <- rbind(sst_ALL_length, sst_ALL_missing, sst_ALL_trended) %>%
  select(test, site, rep, index_vals, clim, event, cat)

## Save and clear
save(sst_ALL_clim_event_cat, file = "data/sst_ALL_clim_event_cat.Rdata")
# rm(sst_ALL_length, sst_ALL_miss, sst_ALL_trend, sst_ALL_clim_event_cat); gc()


# Calculate clim, event, and cats with fixes ------------------------------

## Length data
doMC::registerDoMC(cores = 35)
system.time(
sst_ALL_length_width_10 <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results,
                                       .parallel = T, set_width = 10) %>%
  mutate(test = as.factor("length_width_10"))
) # ~105 seconds
# save(sst_ALL_length_width_10, file = "data/sst_ALL_length_width_10.Rdata")
# load("data/sst_ALL_length_width_10.Rdata")
doMC::registerDoMC(cores = 35)
system.time(
sst_ALL_length_width_20 <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results,
                                       .parallel = T, set_width = 20) %>%
  mutate(test = as.factor("length_width_20"))
) # ~104 seconds
# save(sst_ALL_length_width_20, file = "data/sst_ALL_length_width_20.Rdata")
# load("data/sst_ALL_length_width_20.Rdata")
doMC::registerDoMC(cores = 35)
system.time(
sst_ALL_length_width_30 <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results,
                                       .parallel = T, set_width = 30) %>%
  mutate(test = as.factor("length_width_30"))
) # ~104 seconds
# save(sst_ALL_length_width_30, file = "data/sst_ALL_length_width_30.Rdata")
# load("data/sst_ALL_length_width_30.Rdata")

## Missing data
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_missing_fix <- plyr::ddply(sst_ALL_knockout, c("site", "rep", "index_vals"),
                                   clim_event_cat_calc, .parallel = T, fix = "missing") %>%
  mutate(test = as.factor("missing_fix"))
) # ~xxx seconds
save(sst_ALL_missing, file = "data/sst_ALL_missing_fix.Rdata")
# load("data/sst_ALL_missing_fix.Rdata")

## Combine all results
sst_ALL_clim_event_cat_fix <- rbind(sst_ALL_length_width_10, sst_ALL_length_width_20,
                                    sst_ALL_length_width_30, sst_ALL_missing_fix) %>%
  select(test, site, rep, index_vals, clim, event, cat)

## Save and clear
save(sst_ALL_clim_event_cat_fix, file = "data/sst_ALL_clim_event_cat_fix.Rdata")
# rm(sst_ALL_length_width_10, sst_ALL_length_width_15,
#    sst_ALL_length_width_20, sst_ALL_missing_fix); gc()


# Climatology results -----------------------------------------------------

## Kolmogorov-Smirnov tests for similarity of climatologies
# sub-optimal data
doMC::registerDoMC(cores = 50)
sst_ALL_KS_clim <- plyr::ddply(sst_ALL_clim_event_cat[ ,c(1:5)],
                               c("test", "site", "rep"), KS_p, .parallel = T)
# Save and clear
save(sst_ALL_KS_clim, file = "data/sst_ALL_KS_clim.Rdata")
# rm(sst_ALL_KS_clim); gc()

# Fixed data
doMC::registerDoMC(cores = 50)
sst_ALL_KS_clim_fix <- plyr::ddply(sst_ALL_clim_event_cat_fix[ ,c(1:5)],
                                   c("test", "site", "rep"), KS_p, .parallel = T)
# Save and clear
save(sst_ALL_KS_clim_fix, file = "data/sst_ALL_KS_clim_fix.Rdata")
# rm(sst_ALL_KS_clim); gc()


# Event metric results ----------------------------------------------------

## Kolmogorov-Smirnov tests for similarity of climatologies
# sub-optimal data
doMC::registerDoMC(cores = 50)
sst_ALL_KS_event <- plyr::ddply(sst_ALL_clim_event_cat[ ,c(1:4,6)],
                               c("test", "site", "rep"), KS_p, .parallel = T)
# Save and clear
save(sst_ALL_KS_event, file = "data/sst_ALL_KS_event.Rdata")
# rm(sst_ALL_KS_event); gc()

# Fixed data
doMC::registerDoMC(cores = 50)
sst_ALL_KS_event_fix <- plyr::ddply(sst_ALL_clim_event_cat_fix[ ,c(1:4,6)],
                                    c("test", "site", "rep"), KS_p, .parallel = T)
# Save and clear
save(sst_ALL_KS_event_fix, file = "data/sst_ALL_KS_event_fix.Rdata")
# rm(sst_ALL_KS_event); gc()


## AOV and Tukey
## NB: Decided against tests for central tendency
# doMC::registerDoMC(cores = 50)
# sst_ALL_aov_tukey <- plyr::ddply(sst_ALL_clim_event_cat, c("test", "site", "rep"),
#                                  aov_tukey, .parallel = T)

# test that it ran correctly
# names(sst_ALL_aov_tukey)
# unnest(slice(sst_ALL_aov_tukey, 1), aov)
# unnest(slice(sst_ALL_aov_tukey, 1), tukey)

# Save and clear
# save(sst_ALL_aov_tukey, file = "data/sst_ALL_aov_tukey.Rdata")
# rm(sst_ALL_aov_tukey); gc()


# Category count results --------------------------------------------------

## KS test for proportion of days within each category
# Sub-optimal data
doMC::registerDoMC(cores = 50)
sst_ALL_KS_cat <- plyr::ddply(sst_ALL_clim_event_cat[ ,c(1:4,7)],
                              c("test", "site", "rep"), KS_p, .parallel = T)
# Save and clear
save(sst_ALL_KS_cat, file = "data/sst_ALL_KS_cat.Rdata")
# rm(sst_ALL_KS_event); gc()

# Sub-optimal data
doMC::registerDoMC(cores = 50)
sst_ALL_KS_cat_fix <- plyr::ddply(sst_ALL_clim_event_cat_fix[ ,c(1:4,7)],
                                  c("test", "site", "rep"), KS_p, .parallel = T)
# Save and clear
save(sst_ALL_KS_cat_fix, file = "data/sst_ALL_KS_cat_fix.Rdata")
# rm(sst_ALL_KS_event); gc()


## Fisher
## NB: Decided against testing for probability due to sample size constraints
# doMC::registerDoMC(cores = 50)
# sst_ALL_fisher <- plyr::ddply(sst_ALL_clim_event_cat,
#                               c("test", "site", "rep"), fisher_test, .parallel = T)
# Save and clear
# save(sst_ALL_fisher, file = "data/sst_ALL_fisher.Rdata")
# rm(sst_ALL_fisher); gc()

# Post-hoc
## Nothing at the moment


# Individual event effect -------------------------------------------------

# Load the results from above
load("data/sst_ALL_clim_event_cat.Rdata")
load("data/sst_ALL_clim_event_cat_fix.Rdata")

# Filter out the re-sampled data
sst_ALL_clim_event_cat_rep_1 <- sst_ALL_clim_event_cat %>%
  filter(rep == "1")
rm(sst_ALL_clim_event_cat); gc()

## Climatologies
sst_ALL_clim_rep_1 <- sst_ALL_clim_event_cat_rep_1 %>%
  select(-event, -cat) %>%
  unnest(clim)

## Event metrics
sst_ALL_event_rep_1 <- sst_ALL_clim_event_cat_rep_1 %>%
  select(-clim, -cat) %>%
  unnest(event)

## Categories
sst_ALL_cat_rep_1 <- sst_ALL_clim_event_cat_rep_1 %>%
  select(-clim, -event) %>%
  unnest(cat)

# Extract the control data
  # The 0 trend data are the best choice here
## Climatologies
sst_ALL_clim_control <- sst_ALL_clim_rep_1 %>%
  filter(test == "trended", index_vals == 0)

## Event metrics
sst_ALL_event_control <- sst_ALL_event_rep_1 %>%
  filter(test == "trended", index_vals == 0)

## Categories
sst_ALL_cat_control <- sst_ALL_cat_rep_1 %>%
  filter(test == "trended", index_vals == 30)

# Specify the infamous event
focus_Med <- sst_ALL_event_control %>%
  filter(site == "Med", date_end <= "2005-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative))
focus_WA <- sst_ALL_event_control %>%
  filter(site == "WA", date_end >= "2010-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative))
focus_NW_Atl <- sst_ALL_event_control %>%
  filter(site == "NW_Atl",
         date_start >= "2010-01-01", date_start <= "2014-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative))

# Create infamous event index
focus_ALL <- rbind(focus_Med, focus_NW_Atl, focus_WA) %>%
  select(site, date_start:date_end) %>%
  dplyr::rename(date_start_control = date_start,
                date_peak_control = date_peak,
                date_end_control = date_end)

# Quantify changes caused by the three tests
## Climatologies
effect_clim <- sst_ALL_clim_rep_1 %>%
  select(-doy, -rep) %>%
  gather(key = "metric", value = "val", -site, -test, -index_vals) %>%
  group_by(site, test, index_vals, metric) %>%
  summarise_if(is.numeric,
               .funs = c("min", "median", "mean", "max")) %>%
  # group_by(site, test) %>%
  mutate_if(is.numeric, round, 3)

## Event metrics
effect_event <- sst_ALL_event_rep_1 %>%
  left_join(focus_ALL, by = "site") %>%
  filter(date_peak >= date_start_control,
         date_peak <= date_end_control) %>%
  group_by(site, test, index_vals) %>%
  #
  # Not sure what to do with this information
  mutate(date_start_change = date_start_control - date_start,
         date_peak_change = date_peak_control - date_peak,
         date_end_change = date_end_control - date_end) %>%
  #
  summarise(count = n(),
            duration = sum(duration),
            intensity_mean = mean(intensity_mean),
            intensity_max = max(intensity_max),
            intensity_cumulative = sum(intensity_cumulative)) %>%
  mutate_if(is.numeric, round, 2) %>%
  gather(key = "metric", value = "val", -site, -test, -index_vals) %>%
  ungroup() %>%
  filter(metric %in% c("count", "duration", "intensity_max"),
         !index_vals %in% seq(1, 9)) %>%
  mutate(metric = case_when(metric == "intensity_max" ~ "max. intensity (°C)",
                            metric == "duration" ~ "duration (days)",
                            metric == "count" ~ "count (event)"),
         test = case_when(test == "length" ~ "length (years)",
                          test == "missing" ~ "missing data (proportion)" ,
                          test == "trended" ~ "added trend (°C/dec)"),
         test = as.factor(test),
         test = factor(test, levels = levels(test)[c(2,3,1)]))

## Categories
effect_cat <- sst_ALL_cat_rep_1 %>%
  left_join(focus_ALL, by = "site") %>%
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
  gather(key = "metric", value = "val", -site, -test, -index_vals)


# Visualise
## Climatologies
ggplot(effect_clim, aes(x = index_vals)) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  geom_line(aes(y = mean, colour = site)) +
  # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
  facet_grid(metric~test, scales = "free")


## Event metrics
plot_event_effect <- ggplot(effect_event, aes(x = index_vals)) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
  stat_smooth(aes(y = val, colour = site), geom = "line",
              method = "lm", alpha = 0.5, size = 1) +
  geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1.2) +
  # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
  facet_grid(metric~test, scales = "free", switch = "both") +
  labs(x = NULL, y = NULL, colour = "Site") +
  theme(legend.position = "bottom")
plot_event_effect
ggsave(plot_event_effect, filename = "output/effect_event.pdf", height = 5, width = 10)


## Categories
ggplot(effect_cat, aes(x = index_vals)) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
  stat_smooth(aes(y = val, colour = site), geom = "line",
              method = "lm", alpha = 0.5, size = 1) +
  geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1.2) +
  # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
  facet_grid(metric~test, scales = "free")

# Calculate the slope of the change caused by the test with a linear model
# Also calculate the slope of the upper and lower SE of the linear model


# Global ------------------------------------------------------------------

# And now we run the above workflow on the global data
