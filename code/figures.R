# code/figures.R
# This script houses the code used for the final publication figures

source("code/functions.R")


# Figure 1 ----------------------------------------------------------------

# Show the 3 focus MHWs and the category layers

# Combine the reference time series
sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
                 mutate(sst_NW_Atl, site = "NW_Atl"),
                 mutate(sst_Med, site = "Med"))

# Load reference results
sst_ALL_results <- readRDS("data/sst_ALL_results.Rds")

# Calculate base results
sst_ALL_res <- sst_ALL %>%
  group_by(site) %>%
  nest() %>%
  mutate(clims = map(data, ts2clm,
                     climatologyPeriod = c("1982-01-01", "2011-12-31")),
         events = map(clims, detect_event),
         cats = map(events, category)) %>%
  select(-data)

# Climatologies
sst_ALL_clim <- sst_ALL_res %>%
  select(-cats) %>%
  unnest(events) %>%
  filter(row_number() %% 2 == 1) %>%
  unnest(events)

# Climatologies doy
sst_ALL_clim_only <- sst_ALL_res %>%
  select(-cats) %>%
  unnest(events) %>%
  filter(row_number() %% 2 == 1) %>%
  unnest(events) %>%
  select(site, doy, seas, thresh) %>%
  unique() %>%
  arrange(site, doy)

# Events
sst_ALL_event <- sst_ALL_res %>%
  select(-cats) %>%
  unnest(events) %>%
  filter(row_number() %% 2 == 0) %>%
  unnest(events)

# Categories
sst_ALL_cat <- sst_ALL_res %>%
  select(site, cats) %>%
  unnest()

# Manually grab each event and combine
focus_WA <- sst_ALL_event %>%
  filter(site == "WA", date_end <= "2013-01-01") %>%
  filter(intensity_max == max(intensity_max)) %>%
  left_join(sst_ALL_coords, by = "site")
focus_NW_Atl <- sst_ALL_event %>%
  filter(site == "NW_Atl", date_end <= "2013-01-01") %>%
  filter(intensity_max == max(intensity_max)) %>%
  left_join(sst_ALL_coords, by = "site")
focus_Med <- sst_ALL_event %>%
  filter(site == "Med", date_end <= "2004-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative)) %>%
  left_join(sst_ALL_coords, by = "site")
focus_ALL <- rbind(focus_WA, focus_NW_Atl, focus_Med)

# Merge with results for better plotting
sst_ALL_clim <- left_join(sst_ALL_clim, focus_ALL[,c("site", "date_start", "date_peak", "date_end")], by = "site") %>%
  mutate(site_label = case_when(site == "WA" ~ "A",
                                site == "NW_Atl" ~ "B",
                                site == "Med" ~ "C"))

# Create the three panelled events
fig_1 <- fig_1_plot(sst_ALL_clim, 183)
ggsave(fig_1, filename = "LaTeX/fig_1.pdf", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.png", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.jpg", width = 8, height = 8)

# Create flattened figure for poster
focus_WA_flat <- fig_1_plot(filter(sst_ALL_clim, site == "WA"), 183, y_label = NULL)
focus_NW_Atl_flat <- fig_1_plot(filter(sst_ALL_clim, site == "NW_Atl"), 183, y_label = "Temp. (°C)")
focus_Med_flat <- fig_1_plot(filter(sst_ALL_clim, site == "Med"), 183, y_label = NULL)
fig_1_flat <- ggarrange(focus_WA_flat, focus_NW_Atl_flat, focus_Med_flat, align = "v",
                        ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(fig_1_flat, filename = "LaTeX/fig_1_flat.png", width = 8, height = 5)


# Figure 2 ----------------------------------------------------------------

# The results from the length tests


# Combine for plotting
sst_ALL_plot_long <- rbind(sst_ALL_KS_clim_long, sst_ALL_KS_event_long,
                           sst_ALL_KS_cat_long, sst_ALL_KS_clim_fix_long,
                           sst_ALL_KS_event_fix_long, sst_ALL_KS_cat_fix_long) %>%
  filter(!metric %in% c("intensity_mean", "intensity_cumulative"),
         index_vals <= 0.5 | index_vals >= 10) %>%
  mutate(metric = factor(metric, levels = c("seas", "thresh",
                                            "duration", "intensity_max")))

# Create plugs for control time series
control_plug_WA <- control_plug("length", "WA")
control_plug_NW_Atl <- control_plug("length", "NW_Atl")
control_plug_Med <- control_plug("length", "Med")

# Plug-em
sst_ALL_plot_long <- rbind(sst_ALL_plot_long, control_plug_WA, control_plug_NW_Atl, control_plug_Med)

# Add label for panels
sst_ALL_plot_long <- sst_ALL_plot_long %>%
  mutate(site_label = case_when(site == "WA" ~ "A",
                                site == "NW_Atl" ~ "B",
                                site == "Med" ~ "C"))

# Create figure and save
fig_2 <- fig_line_plot("length")
ggsave(plot = fig_2, filename = "LaTeX/fig_2.pdf", height = 4, width = 8)
ggsave(plot = fig_2, filename = "LaTeX/fig_2.png", height = 4, width = 8)
ggsave(plot = fig_2, filename = "LaTeX/fig_2.jpg", height = 4, width = 8)


# Figure 3 ----------------------------------------------------------------

# The effects of the three sub-optimal tests on the focus MHW metrics

# This is currently made in the workflow.R script


# Figure 4 ----------------------------------------------------------------

# Create a map at which the year after which a threshold is exceeded in the change in the statistic in question


# The global effect of length on max. intensity

# This figure is currently made in the global.R script


# Figure 5 ----------------------------------------------------------------

# The global effect of length on duration

# This figure is currently made in the global.R script


# Figure 6 ----------------------------------------------------------------

# The results from the 100 re-sampled missing data tests

# Create and save
fig_6 <- fig_line_plot("missing")
ggsave(plot = fig_6, filename = "LaTeX/fig_6.pdf", height = 4, width = 8)
ggsave(plot = fig_6, filename = "LaTeX/fig_6.png", height = 4, width = 8)
ggsave(plot = fig_6, filename = "LaTeX/fig_6.jpg", height = 4, width = 8)


# Figure 7 ----------------------------------------------------------------

# The fix for missing data

# Create and save
fig_7 <- fig_line_plot("missing_fix")
ggsave(plot = fig_7, filename = "LaTeX/fig_7.pdf", height = 4, width = 8)
ggsave(plot = fig_7, filename = "LaTeX/fig_7.png", height = 4, width = 8)
ggsave(plot = fig_7, filename = "LaTeX/fig_7.jpg", height = 4, width = 8)


# Figure 8 ----------------------------------------------------------------

# The results from the 100 re-sampled added trend tests

# Create and save
fig_8 <- fig_line_plot("trended")
ggsave(plot = fig_8, filename = "LaTeX/fig_8.pdf", height = 4, width = 8)
ggsave(plot = fig_8, filename = "LaTeX/fig_8.png", height = 4, width = 8)
ggsave(plot = fig_8, filename = "LaTeX/fig_8.jpg", height = 4, width = 8)


# Figure 9 ----------------------------------------------------------------



# Figure 10 ---------------------------------------------------------------

# The effect of added trends on the focus MHWs

# Include the maps showing the results with the trends left in as supplementary material









# Figure illustrating the change caused in the 90th perc. thresh.
# against the seas. clim. in the reference time series

load("data/sst_ALL_clim_event_cat.Rdata")

clim_only <- sst_ALL_clim_event_cat %>%
  # filter(rep == "1") %>%
  select(test:clim) %>%
  unnest() %>%
  filter(index_vals <= 30,
         test == "length") %>%
  gather(key = "metric", value = "temp", seas, thresh) %>%
  droplevels()

clim_change_plot <- ggplot(data = clim_only, aes(x = doy, y = temp)) +
  geom_line(aes(colour = index_vals, group = index_vals)) +
  scale_colour_viridis_c(direction = -1) +
  facet_grid(metric~site)
clim_change_plot


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

# Percent change away from control values for base three tests
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

# Percent change away from control values with widened windows
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

