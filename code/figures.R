# code/figures.R
# This script houses the code used for the final publication figures

source("code/functions.R")


# Test visuals ------------------------------------------------------------

# The results from all tests
sst_ALL_results <- readRDS("data/sst_ALL_results.Rds")

# Overall mean values
sst_ALL_results %>%
  filter(id == "mean") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  # geom_point(data = filter(sst_ALL_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free")# +
# coord_cartesian(ylim = c(-2, 2))

# Percent change away from control values for base three tests
sst_ALL_results %>%
  filter(test %in% c("length", "missing", "trend"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  # geom_point(data = filter(sst_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Percent change away from control values with interpolation
sst_ALL_results %>%
  filter(test %in% c("missing", "interp"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  # geom_point(data = filter(sst_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Change in count of MHWs from interpolation
sst_ALL_results %>%
  filter(test %in% c("missing", "interp"),
         id == "n_diff") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour= "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  facet_grid(var~test, scales = "free")

# Percent change away from control values with widened windows
sst_ALL_results %>%
  filter(test %in% c("length", "window_10", "window_20", "window_30"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  # geom_point(data = filter(sst_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Change in count of MHWs from window size
sst_ALL_results %>%
  filter(test %in% c("length", "window_10", "window_20", "window_30"),
         id == "n_diff") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour= "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  facet_grid(var~test, scales = "free")

# Change in percent of MHW days
sst_ALL_results %>%
  filter(test %in% c("length", "window_10", "window_20", "window_30"),
         id == "sum_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  facet_grid(var~test, scales = "free")

# Change in focus event
sst_ALL_results %>%
  filter(var == "focus_duration", id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  facet_grid(var~test, scales = "free")


# Figure 1 ----------------------------------------------------------------

# Show the 3 focus MHWs

# It appears as though it would be best to show the largest MHW that occurred in the last ten years
# This will ensure more consistency throughout the results

# Combine the reference time series
sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
                 mutate(sst_NW_Atl, site = "NW_Atl"),
                 mutate(sst_Med, site = "Med"))

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

# Add label for panels
# sst_ALL_plot_long <- sst_ALL_plot_long %>%
#   mutate(site_label = case_when(site == "WA" ~ "A",
#                                 site == "NW_Atl" ~ "B",
#                                 site == "Med" ~ "C"))

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

