# code/figures.R
# This script houses the code used for the final publication figures


# Setup -------------------------------------------------------------------

source("code/functions.R")


# Reference test visuals --------------------------------------------------

# The results from all tests
sst_ALL_results <- readRDS("data/sst_ALL_results.Rda")

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


# Global test visuals -----------------------------------------------------

slice_res <- readRDS("data/global/slice_0001.Rda")

# Overall mean values
slice_res %>%
  filter(id == "mean") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  # geom_point(data = filter(sst_ALL_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free")# +
# coord_cartesian(ylim = c(-2, 2))

# Percent change away from control values for base three tests
slice_res %>%
  filter(test %in% c("length", "missing", "trend"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  # geom_point(data = filter(sst_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Percent change away from control values with interpolation
slice_res %>%
  filter(test %in% c("missing", "interp"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 2, alpha = 0.7) +
  # geom_point(data = filter(sst_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Change in count of MHWs from interpolation
slice_res %>%
  filter(test %in% c("missing", "interp"),
         id == "n_diff") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour= "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  facet_grid(var~test, scales = "free")

# Percent change away from control values with widened windows
slice_res %>%
  filter(test %in% c("length", "window_10", "window_20", "window_30"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  # geom_point(data = filter(sst_summary, difference == TRUE),
  # colour = "red", size = 1, alpha = 0.4) +
  # scale_colour_manual(values = c("black", "red")) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Change in focus event
slice_res %>%
  filter(var == "focus_duration", id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  facet_grid(var~test, scales = "free")


# Figure 1 ----------------------------------------------------------------
# Show the 3 focus MHWs, which also provides a reference for what MHWs look like

# A choice was made to show the largest MHW in the past 10 years,
# rather than the focus MHWs from the literature
# The WA MHW is still the same but the NWA and Med events are different

# Combine the reference time series
sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
                 mutate(sst_NW_Atl, site = "NW_Atl"),
                 mutate(sst_Med, site = "Med"))

# Calculate clims
sst_ALL_clim <- sst_ALL %>%
  group_by(site) %>%
  group_modify(~ts2clm(.x, climatologyPeriod = c("1982-01-01", "2011-12-31")))

# Climatologies doy
sst_ALL_clim_only <- sst_ALL_clim %>%
  select(-t, -temp) %>%
  unique()

# Calculate events
sst_ALL_event <- sst_ALL_clim %>%
  group_by(site) %>%
  group_modify(~detect_event(.x)$event)

# Find largest event in most recent ten years of data
focus_ALL <- sst_ALL_event %>%
  group_by(site) %>%
  filter(intensity_cumulative == max(intensity_cumulative)) %>%
  ungroup()

# The MHWs from the literature
# focus_WA <- sst_ALL_event %>%
#   filter(site == "WA", date_end <= "2013-01-01") %>%
#   filter(intensity_max == max(intensity_max)) %>%
#   left_join(sst_ALL_coords, by = "site")
# focus_NW_Atl <- sst_ALL_event %>%
#   filter(site == "NW_Atl", date_end <= "2013-01-01") %>%
#   filter(intensity_max == max(intensity_max)) %>%
#   left_join(sst_ALL_coords, by = "site")
# focus_Med <- sst_ALL_event %>%
#   filter(site == "Med", date_end <= "2004-01-01") %>%
#   filter(intensity_cumulative == max(intensity_cumulative)) %>%
#   left_join(sst_ALL_coords, by = "site")
# focus_ALL <- rbind(focus_WA, focus_NW_Atl, focus_Med)

# Merge with results for better plotting
sst_ALL_focus <- left_join(sst_ALL_clim,
                           focus_ALL[,c("site", "date_start", "date_peak", "date_end")], by = "site") %>%
  mutate(site_label = case_when(site == "WA" ~ "A",
                                site == "NW_Atl" ~ "B",
                                site == "Med" ~ "C"))

# Create the three panelled events
fig_1 <- fig_1_plot(sst_ALL_focus, 250)
ggsave(fig_1, filename = "LaTeX/fig_1.pdf", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.png", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.jpg", width = 8, height = 8)

# Create flattened figure for poster
focus_WA_flat <- fig_1_plot(filter(sst_ALL_clim, site == "WA"), 183, y_label = NULL)
focus_NW_Atl_flat <- fig_1_plot(filter(sst_ALL_clim, site == "NW_Atl"), 250, y_label = "Temp. (°C)")
focus_Med_flat <- fig_1_plot(filter(sst_ALL_clim, site == "Med"), 183, y_label = NULL)
fig_1_flat <- ggarrange(focus_WA_flat, focus_NW_Atl_flat, focus_Med_flat, align = "v",
                        ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(fig_1_flat, filename = "LaTeX/fig_1_flat.png", width = 8, height = 5)


# Figure 2 ----------------------------------------------------------------
# The effect of the three base tests on the mean of the last 10 years of MHWs

# The results from all tests
sst_ALL_results <- readRDS("data/sst_ALL_results.Rda") %>%
  mutate(site = ifelse(site == "NW_Atl", "NWA", site))
random_results <- readRDS("data/random_results_100.Rda") %>%
  unite("site", c(lon, lat))
full_results <- rbind(sst_ALL_results, random_results)
rm(sst_ALL_results, random_results); gc()

# Create figure and save
fig_2 <- fig_line_plot(full_results, tests = "base", result_choice = "10_years")
ggsave("LaTeX/fig_2.pdf", fig_2, height = 8, width = 10)
ggsave("LaTeX/fig_2.png", fig_2, height = 8, width = 10)
ggsave("LaTeX/fig_2.jpg", fig_2, height = 8, width = 10)


# Figure 3 ----------------------------------------------------------------
# The effects of the three sub-optimal tests on the focus MHW metrics

fig_3 <- fig_line_plot(full_results, tests = "base", result_choice = "focus")
ggsave("LaTeX/fig_3.pdf", fig_3, height = 8, width = 10)
ggsave("LaTeX/fig_3.png", fig_3, height = 8, width = 10)
ggsave("LaTeX/fig_3.jpg", fig_3, height = 8, width = 10)


# Figure 4 ----------------------------------------------------------------
# The global effect of length on count, max. intensity, and duration

fig_4A <- trend_plot(test_sub = "length", var_sub = "count")
fig_4B <- trend_plot(test_sub = "length", var_sub = "duration")
fig_4C <- trend_plot(test_sub = "length", var_sub = "intensity_max")
fig_4D <- trend_plot(test_sub = "length", var_sub = "focus_count")
fig_4E <- trend_plot(test_sub = "length", var_sub = "focus_duration")
fig_4F <- trend_plot(test_sub = "length", var_sub = "focus_intensity_max")
fig_4 <- ggarrange(fig_4A, fig_4D, fig_4B, fig_4E, fig_4C, fig_4F, ncol = 2, nrow = 3,
                   labels = c("A", "D", "B", "E", "C", "F"))
ggsave("LaTeX/fig_4.pdf", fig_4, height = 12, width = 16)
ggsave("LaTeX/fig_4.png", fig_4, height = 12, width = 16)
ggsave("LaTeX/fig_4.jpg", fig_4, height = 12, width = 16)


# Figure 5 ----------------------------------------------------------------
# The effect of interpolating missing data

fig_5 <- fig_line_plot(tests = "miss_comp", result_choice = "10_years")
ggsave("LaTeX/fig_5.pdf", fig_5, height = 8, width = 7)
ggsave("LaTeX/fig_5.png", fig_5, height = 8, width = 7)
ggsave("LaTeX/fig_5.jpg", fig_5, height = 8, width = 7)


# Figure 6 ----------------------------------------------------------------


# Appendix 1 --------------------------------------------------------------
# The difference between the proper 30 year clim and all other 30 year clim periods
# the difference between the focus results when using a 30 year base period vs. a 37 year base period


# Appendix 2 --------------------------------------------------------------
# The effect of the base tests on seas/thresh

# NB: This still needs some work

app_2 <- fig_line_plot(tests = "base", result_choice = "clims")
ggsave("LaTeX/app_2.pdf", app_2, height = 8, width = 7)
ggsave("LaTeX/app_2.png", app_2, height = 8, width = 7)
ggsave("LaTeX/app_2.jpg", app_2, height = 8, width = 7)


# Appendix 3 --------------------------------------------------------------

# A: The global effect of missing data on count, max. intensity, and duration
app_3A_A <- trend_plot(test_sub = "missing", var_sub = "count")
app_3A_B <- trend_plot(test_sub = "missing", var_sub = "duration")
app_3A_C <- trend_plot(test_sub = "missing", var_sub = "intensity_max")
app_3A_D <- trend_plot(test_sub = "missing", var_sub = "focus_count")
app_3A_E <- trend_plot(test_sub = "missing", var_sub = "focus_duration")
app_3A_F <- trend_plot(test_sub = "missing", var_sub = "focus_intensity_max")
app_3A <- ggarrange(app_3A_A, app_3A_D, app_3A_B, app_3A_E, app_3A_C, app_3A_F, ncol = 2, nrow = 3,
                   labels = c("A", "D", "B", "E", "C", "F"))
ggsave("LaTeX/app_3A.pdf", app_3A, height = 12, width = 16)
ggsave("LaTeX/app_3A.png", app_3A, height = 12, width = 16)
ggsave("LaTeX/app_3A.jpg", app_3A, height = 12, width = 16)

# B: The global effect of decadal trend on count, max. intensity, and duration
app_3B_A <- trend_plot(test_sub = "trend", var_sub = "count")
app_3B_B <- trend_plot(test_sub = "trend", var_sub = "duration")
app_3B_C <- trend_plot(test_sub = "trend", var_sub = "intensity_max")
app_3B_D <- trend_plot(test_sub = "trend", var_sub = "focus_count")
app_3B_E <- trend_plot(test_sub = "trend", var_sub = "focus_duration")
app_3B_F <- trend_plot(test_sub = "trend", var_sub = "focus_intensity_max")
app_3B <- ggarrange(app_3B_A, app_3B_D, app_3B_B, app_3B_E, app_3B_C, app_3B_F, ncol = 2, nrow = 3,
                   labels = c("A", "D", "B", "E", "C", "F"))
ggsave("LaTeX/app_3B.pdf", app_3B, height = 12, width = 16)
ggsave("LaTeX/app_3B.png", app_3B, height = 12, width = 16)
ggsave("LaTeX/app_3B.jpg", app_3B, height = 12, width = 16)


# Appendix 4 --------------------------------------------------------------
# The effect of widening the clim window


