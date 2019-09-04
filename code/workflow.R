# workflow.R


# Libraries ---------------------------------------------------------------

source("code/functions.R")



# Length experiment -------------------------------------------------------

# Create anomaly and detrended time series
sst_anom <- sst_WA %>%
  mutate(temp = round(temp-mean(temp), 2))
sst_flat <- detrend(sst_WA)

# Visualise
# ggplot(sst_anom, aes(x = t, y = temp)) +
#   geom_line() +
#   geom_smooth(method = "lm")

# Calculate shortened time series results
sst_anom_length <- plyr::ldply(1982:2016, shrinking_results,
                               df = sst_anom, .parallel = T)
sst_flat_length <- plyr::ldply(1982:2016, shrinking_results,
                               df = sst_flat, .parallel = T)

# Summary statistics
sst_anom_length_clim <- sst_anom_length %>%
  select(index_vals, clim) %>%
  unnest() %>%
  filter(index_vals == 37) %>%
  select(-index_vals, - doy) %>%
  summarise_all(c("min", "median", "mean", "max", "sd"), na.rm = T) %>%
  gather(key = "var", value = "val") %>%
  separate(var, c("stat", "test")) %>%
  mutate(val = round(val, 3))



# Base data ---------------------------------------------------------------

# Combine the three reference time series into one file
sst_ALL <- rbind(sst_Med, sst_NW_Atl, sst_WA) %>%
  mutate(site = as.factor(rep(c("Med", "NW_Atl", "WA"), each = nrow(sst_WA))))

# Save and clear
save(sst_ALL, file = "data/sst_ALL.Rdata")
# rm(sst_ALL); gc()


# De-trend the data -------------------------------------------------------

# Create de-trended anomaly time series from all re-samples
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_flat <- plyr::ddply(sst_ALL_repl, c("site", "rep"), detrend)
) # ~30 seconds
# Save and clear
save(sst_ALL_flat, file = "data/sst_ALL_flat.Rdata")
# rm(sst_ALL_flat); gc()


# Shorten the data --------------------------------------------------------



# Knockout random days ----------------------------------------------------

# Randomly knockout 0 - 99% of each of the 100 re-samples
doMC::registerDoMC(cores = 25)
set.seed(666)
system.time(
sst_ALL_knockout <- plyr::ldply(seq(0.00, 0.99, 0.01), random_knockout, .parallel = T)
) # ~194 seconds
# Save and clear
save(sst_ALL_knockout, file = "data/sst_ALL_knockout.Rdata")
# rm(sst_ALL_knockout); gc()


# Count consecutive missing days ------------------------------------------

# Calculate the consecutive missing days
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_consec <- plyr::ddply(filter(sst_ALL_knockout, index_vals != 0),
                              c("site", "rep", "index_vals"), con_miss, .parallel = T)
) # ~614 seconds
# Save and clear
save(sst_ALL_consec, file = "data/sst_ALL_consec.Rdata")
# rm(sst_ALL_consec); gc()


# Add trends --------------------------------------------------------------

# Add the trends
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_add_trend <- plyr::ldply(seq(0.00, 0.50, 0.01), add_trend, .parallel = T)
) # ~74 seconds
# Save and clear
save(sst_ALL_add_trend, file = "data/sst_ALL_add_trend.Rdata")
# rm(sst_ALL_add_trend); gc()


# Calculate clims, events, and categories ---------------------------------

# NB: These are too large to be run as one command
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
sst_ALL_length <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results, .parallel = T) #%>%
  # mutate(test = as.factor("length"))
)  # ~101 seconds
# NB: Adding the test column in a pipe seems to make ldply sad...
sst_ALL_length$test <- as.factor("length")
# test that it ran correctly
# names(sst_ALL_length)
# unnest(slice(sst_ALL_length, 1), clim)
# unnest(slice(sst_ALL_length, 1), event)
# unnest(slice(sst_ALL_length, 1), cat)
save(sst_ALL_length, file = "data/sst_ALL_length.Rdata")
# load("data/sst_ALL_length.Rdata")

## Missing data
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_missing <- plyr::ddply(sst_ALL_knockout, c("site", "rep", "index_vals"),
                               clim_event_cat_calc, .parallel = T) #%>%
  # mutate(test = as.factor("missing"))
) # ~411 seconds
sst_ALL_missing$test <- as.factor("missing")
save(sst_ALL_missing, file = "data/sst_ALL_missing.Rdata")
# load("data/sst_ALL_missing.Rdata")

## Trended data
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_trended <- plyr::ddply(sst_ALL_add_trend, c("site", "rep", "index_vals"),
                               clim_event_cat_calc, .parallel = T) #%>%
  # mutate(test = as.factor("trended"))
) # 196 seconds
sst_ALL_trended$test <- as.factor("trended")
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
  dplyr::mutate(test = as.factor("length_width_10"))
) # ~108 seconds
save(sst_ALL_length_width_10, file = "data/sst_ALL_length_width_10.Rdata")
# load("data/sst_ALL_length_width_10.Rdata")

doMC::registerDoMC(cores = 35)
system.time(
sst_ALL_length_width_20 <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results,
                                       .parallel = T, set_width = 20) %>%
  dplyr::mutate(test = as.factor("length_width_20"))
) # ~104 seconds
save(sst_ALL_length_width_20, file = "data/sst_ALL_length_width_20.Rdata")
# load("data/sst_ALL_length_width_20.Rdata")

doMC::registerDoMC(cores = 35)
system.time(
sst_ALL_length_width_30 <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results,
                                       .parallel = T, set_width = 30) %>%
  dplyr::mutate(test = as.factor("length_width_30"))
) # ~104 seconds
save(sst_ALL_length_width_30, file = "data/sst_ALL_length_width_30.Rdata")
# load("data/sst_ALL_length_width_30.Rdata")

## NB: A window of 40 or greater breaks down for some reason
# doMC::registerDoMC(cores = 35)
# system.time(
  # sst_ALL_length_width_40 <- plyr::ldply(.data = seq(1982, 2016), .fun = shrinking_results,
                                         # .parallel = T, set_width = 40) %>%
    # dplyr::mutate(test = as.factor("length_width_40"))
# ) # ~104 seconds
# save(sst_ALL_length_width_40, file = "data/sst_ALL_length_width_40.Rdata")
# load("data/sst_ALL_length_width_40.Rdata")

## Missing data
doMC::registerDoMC(cores = 50)
system.time(
sst_ALL_missing_fix <- plyr::ddply(sst_ALL_knockout, c("site", "rep", "index_vals"),
                                   clim_event_cat_calc, .parallel = T, fix = "missing") %>%
  dplyr::mutate(test = as.factor("missing_fix"))
) # ~405 seconds
save(sst_ALL_missing, file = "data/sst_ALL_missing_fix.Rdata")
# load("data/sst_ALL_missing_fix.Rdata")

## Combine all results
sst_ALL_clim_event_cat_fix <- rbind(sst_ALL_length_width_10, sst_ALL_length_width_20,
                                    sst_ALL_length_width_30, #sst_ALL_length_width_40,
                                    sst_ALL_missing_fix) %>%
  select(test, site, rep, index_vals, clim, event, cat)

## Save and clear
save(sst_ALL_clim_event_cat_fix, file = "data/sst_ALL_clim_event_cat_fix.Rdata")
# rm(sst_ALL_length_width_10, sst_ALL_length_width_20,
#    sst_ALL_length_width_30, sst_ALL_length_width_50, sst_ALL_missing_fix); gc()


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

## Kolmogorov-Smirnov tests for similarity of event metrics
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


# Individual event effect -------------------------------------------------

# Load the results from above
load("data/sst_ALL_clim_event_cat.Rdata")
load("data/sst_ALL_clim_event_cat_fix.Rdata")

# Specify the infamous event
focus_Med <- sst_ALL_clim_event_cat %>%
  filter(rep == "1", test == "trended", index_vals == 0) %>%
  select(-clim, -cat) %>%
  unnest(event) %>%
  filter(site == "Med",
         date_end <= "2005-01-01", date_start >= "2000-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative))
focus_WA <-  sst_ALL_clim_event_cat %>%
  filter(rep == "1", test == "trended", index_vals == 0) %>%
  select(-clim, -cat) %>%
  unnest(event) %>%
  filter(site == "WA", date_end >= "2010-01-01") %>%
  filter(intensity_max == max(intensity_max))
focus_NW_Atl <-  sst_ALL_clim_event_cat %>%
  filter(rep == "1", test == "trended", index_vals == 0) %>%
  select(-clim, -cat) %>%
  unnest(event) %>%
  filter(site == "NW_Atl",
         date_start >= "2010-01-01", date_start <= "2014-01-01") %>%
  filter(intensity_max == max(intensity_max))

# Create infamous event index
focus_ALL <- rbind(focus_Med, focus_NW_Atl, focus_WA) %>%
  select(site, date_start:date_end) %>%
  dplyr::rename(date_start_control = date_start,
                date_peak_control = date_peak,
                date_end_control = date_end)

# Quantify changes caused by the three tests
## Climatologies
effect_clim <- effect_clim_func(sst_ALL_clim_event_cat)
effect_clim_fix <- effect_clim_func(sst_ALL_clim_event_cat_fix)

## Event metrics
effect_event <- effect_event_func(sst_ALL_clim_event_cat, focus_ALL)
effect_event_fix <- effect_event_func(sst_ALL_clim_event_cat_fix, focus_ALL)

## Categories
effect_cat <- effect_cat_func(sst_ALL_clim_event_cat, focus_ALL)
effect_cat_fix <- effect_cat_func(sst_ALL_clim_event_cat_fix, focus_ALL)


# Prep event data for pretty plotting
effect_event_pretty <- effect_event %>%
  filter(metric %in% c("count", "duration", "intensity_max"),
         !index_vals %in% seq(1, 9),
         index_vals <= 0.5 | index_vals >= 10) %>%
  mutate(panel_label = case_when(metric == "count" & test == "length" ~ "A",
                                 metric == "count" & test == "missing" ~ "B",
                                 metric == "count" & test == "trended" ~ "C",
                                 metric == "duration" & test == "length" ~ "D",
                                 metric == "duration" & test == "missing" ~ "E",
                                 metric == "duration" & test == "trended" ~ "F",
                                 metric == "intensity_max" & test == "length" ~ "H",
                                 metric == "intensity_max" & test == "missing" ~ "I",
                                 metric == "intensity_max" & test == "trended" ~ "J"),
         metric = case_when(metric == "intensity_max" ~ "max. intensity (°C)",
                            metric == "duration" ~ "duration (days)",
                            metric == "count" ~ "count (event)"),
         test = case_when(test == "length" ~ "length (years)",
                          test == "missing" ~ "missing data (proportion)" ,
                          test == "trended" ~ "added trend (°C/dec)"),
         test = as.factor(test),
         test = factor(test, levels = levels(test)[c(2,3,1)]),
         site = as.character(site)) %>%
  group_by(test) %>%
  mutate(panel_label_x = min(index_vals)) %>%
  group_by(metric) %>%
  mutate(panel_label_y = max(val)) %>%
  ungroup()
effect_event_pretty$site[effect_event_pretty$site == "NW_Atl"] <- "NWA"
effect_event_pretty$site <- factor(effect_event_pretty$site, levels = c("WA", "NWA", "Med"))
effect_event_pretty$index_vals[effect_event_pretty$test == "missing"] <- effect_event_pretty$index_vals[effect_event_pretty$test == "missing"]*100

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
plot_event_effect <- ggplot(effect_event_pretty, aes(x = index_vals)) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  # geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
  # stat_smooth(aes(y = val, colour = site), geom = "line",
              # method = "lm", alpha = 0.5, size = 1) +
  geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1) +
  geom_text(aes(label = panel_label, y = panel_label_y, x = panel_label_x)) +
  # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
  scale_colour_brewer(palette = "Dark2") +
  facet_grid(metric~test, scales = "free", switch = "both") +
  labs(x = NULL, y = NULL, colour = "Site") +
  theme(legend.position = "bottom")
plot_event_effect
ggsave(plot_event_effect, filename = "output/effect_event.pdf", height = 5, width = 10)
ggsave(plot_event_effect, filename = "output/effect_event.png", height = 5, width = 10)
ggsave(plot_event_effect, filename = "LaTeX/fig_3.pdf", height = 5, width = 10)
ggsave(plot_event_effect, filename = "LaTeX/fig_3.png", height = 5, width = 10)
ggsave(plot_event_effect, filename = "LaTeX/fig_3.jpg", height = 5, width = 10)

# Fixed data
ggplot(effect_event_fix, aes(x = index_vals)) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
  stat_smooth(aes(y = val, colour = site), geom = "line",
              method = "lm", alpha = 0.5, size = 1) +
  geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1.2) +
  # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
  facet_grid(metric~test, scales = "free", switch = "both") +
  labs(x = NULL, y = NULL, colour = "Site") +
  theme(legend.position = "bottom")

# Missing data only + fix

## Categories
# Sub-optimal data
ggplot(effect_cat, aes(x = index_vals)) +
  # geom_ribbon(aes(ymin = min, ymax = max, fill = metric), alpha = 0.2) +
  geom_smooth(aes(y = val, colour = site), method = "lm", linetype = 0) +
  stat_smooth(aes(y = val, colour = site), geom = "line",
              method = "lm", alpha = 0.5, size = 1) +
  geom_line(aes(y = val, colour = site), alpha = 0.7, size = 1.2) +
  # geom_line(aes(y = median, colour = metric), linetype = "dashed") +
  facet_grid(metric~test, scales = "free")
# Fixed data
ggplot(effect_cat_fix, aes(x = index_vals)) +
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
  # Sub-optimal and fixes

# This is done by sourcing the global script
  # NB: This script should rather be called from an R terminal and not here
source("code/global.R")
