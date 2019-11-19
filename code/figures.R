# code/figures.R
# This script houses the code used for the final publication figures
# This includes the supplementary figures
# There is also code at the end used for test visuals of results
# Note that line 11 will likely not run correctly, see 'code/functions.R'
# for an explanation why and a quick fix for one's local machine


# Setup -------------------------------------------------------------------

source("code/functions.R")


# Figure 1 ----------------------------------------------------------------

# Show the 3 focal MHWs, which also provides a reference for what MHWs look like

# A choice was made to show the largest MHW in the past 10 years,
# rather than the focal MHWs from the literature
# The WA MHW is still the same but the NWA and Med events are different

# Combine the reference time series
sst_ALL <- rbind(mutate(sst_WA, site = "WA"),
                 mutate(sst_NW_Atl, site = "NWA"),
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

# The MHWs from the literature, for interest sake
literature_WA <- sst_ALL_event %>%
  filter(site == "WA", date_end <= "2013-01-01") %>%
  filter(intensity_max == max(intensity_max)) %>%
  left_join(sst_ALL_coords, by = "site")
literature_NWA <- sst_ALL_event %>%
  filter(site == "NWA", date_end <= "2013-01-01") %>%
  filter(intensity_max == max(intensity_max)) %>%
  left_join(sst_ALL_coords, by = "site")
literature_Med <- sst_ALL_event %>%
  filter(site == "Med", date_end <= "2004-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative)) %>%
  left_join(sst_ALL_coords, by = "site")

# Merge with results for better plotting
sst_ALL_focus <- left_join(sst_ALL_clim,
                           focus_ALL[,c("site", "date_start", "date_peak", "date_end")], by = "site") %>%
  mutate(site_label = case_when(site == "WA" ~ "A",
                                site == "NWA" ~ "B",
                                site == "Med" ~ "C"))

# Create the three panelled events
fig_1 <- fig_1_plot(sst_ALL_focus, 250)
ggsave(fig_1, filename = "LaTeX/fig_1.pdf", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.png", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.jpg", width = 8, height = 8)

# Create flattened figure for poster
focus_WA_flat <- fig_1_plot(filter(sst_ALL_focus, site == "WA"), 183, y_label = NULL) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
focus_NWA_flat <- fig_1_plot(filter(sst_ALL_focus, site == "NWA"), 250, y_label = "Temp. (°C)") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
focus_Med_flat <- fig_1_plot(filter(sst_ALL_focus, site == "Med"), 183, y_label = NULL) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
fig_1_flat <- ggpubr::ggarrange(focus_WA_flat, focus_NWA_flat, focus_Med_flat, align = "v",
                        ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(fig_1_flat, filename = "LaTeX/fig_1_flat.pdf", width = 8, height = 5)
ggsave(fig_1_flat, filename = "LaTeX/fig_1_flat.png", width = 8, height = 5)
ggsave(fig_1_flat, filename = "LaTeX/fig_1_flat.jpg", width = 8, height = 5)


# Figure 2 ----------------------------------------------------------------

# The effect of the three base tests on the mean of the last 10 years of MHWs
sst_ALL_results <- readRDS("data/sst_ALL_results.Rda") %>%
  mutate(site = ifelse(site == "NW_Atl", "NWA", site))
system.time(
  random_results <- readRDS("data/random_results_1000.Rda") %>%
    unite("site", c(lon, lat))
) # 68 seconds
full_results <- rbind(sst_ALL_results, random_results)
rm(sst_ALL_results, random_results); gc()

# Create figure and save
fig_2 <- fig_box_plot(full_results, tests = "base", result_choice = "average")
ggsave("LaTeX/fig_2.pdf", fig_2, height = 8, width = 12)
ggsave("LaTeX/fig_2.png", fig_2, height = 8, width = 12)
ggsave("LaTeX/fig_2.jpg", fig_2, height = 8, width = 12)


# Figure 3 ----------------------------------------------------------------

# The effects of the three sub-optimal tests on the focal MHW metrics
fig_3 <- fig_box_plot(full_results, tests = "base", result_choice = "focus")
ggsave("LaTeX/fig_3.pdf", fig_3, height = 8, width = 12)
ggsave("LaTeX/fig_3.png", fig_3, height = 8, width = 12)
ggsave("LaTeX/fig_3.jpg", fig_3, height = 8, width = 12)


# Figure 4 ----------------------------------------------------------------

# The global effect of length on count, max. intensity, and duration

# The global trend files
global_var_trend <- readRDS("data/global_var_trend.Rda")

# The panels
fig_4A <- trend_plot(test_sub = "length", var_sub = "count")
fig_4B <- trend_plot(test_sub = "length", var_sub = "duration")
fig_4C <- trend_plot(test_sub = "length", var_sub = "intensity_max")
fig_4D <- trend_plot(test_sub = "length", var_sub = "focus_count")
fig_4E <- trend_plot(test_sub = "length", var_sub = "focus_duration")
fig_4F <- trend_plot(test_sub = "length", var_sub = "focus_intensity_max")

# The full figure
fig_4 <- ggpubr::ggarrange(fig_4A, fig_4D, fig_4B, fig_4E, fig_4C, fig_4F, ncol = 2, nrow = 3,
                           labels = c("A", "D", "B", "E", "C", "F"))
ggsave("LaTeX/fig_4.pdf", fig_4, height = 12, width = 16)
ggsave("LaTeX/fig_4.png", fig_4, height = 12, width = 16)
ggsave("LaTeX/fig_4.jpg", fig_4, height = 12, width = 16)


# Figure 5 ----------------------------------------------------------------

# The effect of widening the clim window

# The effect on the average MHWs
fig_5 <- fig_box_plot(tests = "windows", result_choice = "average")
ggsave("LaTeX/fig_5.pdf", fig_5, height = 8, width = 12)
ggsave("LaTeX/fig_5.png", fig_5, height = 8, width = 12)
ggsave("LaTeX/fig_5.jpg", fig_5, height = 8, width = 12)


# Figure 6 ----------------------------------------------------------------

# The effect of interpolating missing data on the average MHWs
fig_6A <- fig_box_plot(tests = "miss_comp", result_choice = "average")

# The effect of interpolating on the focal MHW
fig_6B <- fig_box_plot(tests = "miss_comp", result_choice = "focus")

# Stick them together
fig_6 <- ggpubr::ggarrange(fig_6A, fig_6B, common.legend = T)
ggsave("LaTeX/fig_6.pdf", fig_6, height = 8, width = 12)
ggsave("LaTeX/fig_6.png", fig_6, height = 8, width = 12)
ggsave("LaTeX/fig_6.jpg", fig_6, height = 8, width = 12)


# Supplementary 1 ---------------------------------------------------------

# The effect of the base tests on seas/thresh
fig_S1 <- fig_box_plot(tests = "base", result_choice = "clims")

# Create the caption and save
fig_S1_cap <-  grid::textGrob(paste0(strwrap(fig_cap_S1, 160), sep = "", collapse = "\n"),
                              x = 0.01, just = "left", gp = grid::gpar(fontsize = 10))
fig_S1_cap <- ggpubr::ggarrange(fig_S1, fig_S1_cap,
                                heights = c(1, 0.2), nrow = 2)
ggsave("LaTeX/fig_S1.pdf", fig_S1_cap, height = 6.5, width = 10)
ggsave("LaTeX/fig_S1.png", fig_S1_cap, height = 6.5, width = 10)
ggsave("LaTeX/fig_S1.jpg", fig_S1_cap, height = 6.5, width = 10)


# Supplementary 2 ---------------------------------------------------------

# The difference between the proper 30 year base period and all other 30 year base periods
sst_ALL_results <- readRDS("data/sst_ALL_bp_results.Rda") %>%
  mutate(site = ifelse(site == "NW_Atl", "NWA", site))
random_results <- readRDS("data/random_bp_results_1000.Rda") %>%
  unite("site", c(lon, lat))
full_results <- rbind(sst_ALL_results, random_results) %>%
  mutate(test = "base_period")
rm(sst_ALL_results, random_results); gc()

# The effect on the 10 years of MHWs
fig_S2A <- fig_box_plot(tests = "base_period", result_choice = "average")

# The effect on the focal MHW
fig_S2B <- fig_box_plot(tests = "base_period", result_choice = "focus")

# Stick them together
fig_S2 <- ggpubr::ggarrange(fig_S2A, fig_S2B, common.legend = T)

# Create the caption and save
fig_S2_cap <-  grid::textGrob(paste0(strwrap(fig_cap_S2, 140), sep = "", collapse = "\n"),
                              x = 0.01, just = "left", gp = grid::gpar(fontsize = 10))
fig_S2_cap <- ggpubr::ggarrange(fig_S2, fig_S2_cap,
                                heights = c(1, 0.1), nrow = 2)
ggsave("LaTeX/fig_S2.pdf", fig_S2_cap, height = 8, width = 9)
ggsave("LaTeX/fig_S2.png", fig_S2_cap, height = 8, width = 9)
ggsave("LaTeX/fig_S2.jpg", fig_S2_cap, height = 8, width = 9)


# Supplementary 3 ---------------------------------------------------------

# The global effect of missing data on count, max. intensity, and duration

# The panels
fig_S3_A <- trend_plot(test_sub = "missing", var_sub = "count")
fig_S3_B <- trend_plot(test_sub = "missing", var_sub = "duration")
fig_S3_C <- trend_plot(test_sub = "missing", var_sub = "intensity_max")
fig_S3_D <- trend_plot(test_sub = "missing", var_sub = "focus_count")
fig_S3_E <- trend_plot(test_sub = "missing", var_sub = "focus_duration")
fig_S3_F <- trend_plot(test_sub = "missing", var_sub = "focus_intensity_max")
fig_S3 <- ggpubr::ggarrange(fig_S3_A, fig_S3_D, fig_S3_B, fig_S3_E, fig_S3_C, fig_S3_F, ncol = 2, nrow = 3,
                            labels = c("A", "D", "B", "E", "C", "F"))

# Create the caption and save
fig_S3_cap <-  grid::textGrob(paste0(strwrap(fig_cap_S3, 200), sep = "", collapse = "\n"),
                              x = 0.01, just = "left", gp = grid::gpar(fontsize = 12))
fig_S3_cap <- ggpubr::ggarrange(fig_S3, fig_S3_cap,
                                heights = c(1, 0.05), nrow = 2)
ggsave("LaTeX/fig_S3.pdf", fig_S3_cap, height = 13, width = 17)
ggsave("LaTeX/fig_S3.png", fig_S3_cap, height = 13, width = 17)
ggsave("LaTeX/fig_S3.jpg", fig_S3_cap, height = 13, width = 17)


# Supplementary 4 ---------------------------------------------------------

# The global effect of long-term trend on count, max. intensity, and duration

# The panels
fig_S4_A <- trend_plot(test_sub = "trend", var_sub = "count")
fig_S4_B <- trend_plot(test_sub = "trend", var_sub = "duration")
fig_S4_C <- trend_plot(test_sub = "trend", var_sub = "intensity_max")
fig_S4_D <- trend_plot(test_sub = "trend", var_sub = "focus_count")
fig_S4_E <- trend_plot(test_sub = "trend", var_sub = "focus_duration")
fig_S4_F <- trend_plot(test_sub = "trend", var_sub = "focus_intensity_max")
fig_S4 <- ggpubr::ggarrange(fig_S4_A, fig_S4_D, fig_S4_B, fig_S4_E, fig_S4_C, fig_S4_F, ncol = 2, nrow = 3,
                            labels = c("A", "D", "B", "E", "C", "F"))

# Create the caption and save
fig_S4_cap <-  grid::textGrob(paste0(strwrap(fig_cap_S4, 200), sep = "", collapse = "\n"),
                              x = 0.01, just = "left", gp = grid::gpar(fontsize = 12))
fig_S4_cap <- ggpubr::ggarrange(fig_S4, fig_S4_cap,
                                heights = c(1, 0.05), nrow = 2)
ggsave("LaTeX/fig_S4.pdf", fig_S4_cap, height = 13, width = 17)
ggsave("LaTeX/fig_S4.png", fig_S4_cap, height = 13, width = 17)
ggsave("LaTeX/fig_S4.jpg", fig_S4_cap, height = 13, width = 17)


# Supplementary 5 ---------------------------------------------------------

# The results from all tests
sst_ALL_results <- readRDS("data/sst_ALL_results.Rda") %>%
  mutate(site = ifelse(site == "NW_Atl", "NWA", site))
system.time(
  random_results <- readRDS("data/random_results_1000.Rda") %>%
    unite("site", c(lon, lat))
) # 68 seconds
full_results <- rbind(sst_ALL_results, random_results)
rm(sst_ALL_results, random_results); gc()

# The effect on the focal MHWs
# The partial changes (e.g. -33%) for the count of events from the 30 year control is possible
# because the changing of the window is already splitting up the focal event into multiples in the control
fig_S5 <- fig_box_plot(tests = "windows", result_choice = "focus")

# Create the caption and save
fig_S5_cap <-  grid::textGrob(paste0(strwrap(fig_cap_S5, 180), sep = "", collapse = "\n"),
                              x = 0.01, just = "left", gp = grid::gpar(fontsize = 10))
fig_S5_cap <- ggpubr::ggarrange(fig_S5, fig_S5_cap,
                                heights = c(1, 0.1), nrow = 2)
ggsave("LaTeX/fig_S5.pdf", fig_S5_cap, height = 9, width = 12)
ggsave("LaTeX/fig_S5.png", fig_S5_cap, height = 9, width = 12)
ggsave("LaTeX/fig_S5.jpg", fig_S5_cap, height = 9, width = 12)


# Reference test visuals --------------------------------------------------

# The results from all tests
sst_ALL_results <- readRDS("data/sst_ALL_results.Rda")

# Overall mean values
sst_ALL_results %>%
  filter(id == "mean") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  facet_grid(var~test, scales = "free")

# Percent change away from control values for base three tests
sst_ALL_results %>%
  filter(test %in% c("length", "missing", "trend"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
  facet_grid(var ~ test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Percent change away from control values with interpolation
sst_ALL_results %>%
  filter(test %in% c("missing", "interp"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(colour = site), size = 2, alpha = 0.7) +
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

# Change in focal event
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
  facet_grid(var~test, scales = "free")

# Percent change away from control values for base three tests
slice_res %>%
  filter(test %in% c("length", "missing", "trend"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Percent change away from control values with interpolation
slice_res %>%
  filter(test %in% c("missing", "interp"),
         id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 2, alpha = 0.7) +
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
  facet_grid(var~test, scales = "free") +
  coord_cartesian(ylim = c(-3, 3))

# Change in focal event
slice_res %>%
  filter(var == "focus_duration", id == "mean_perc") %>%
  ggplot(aes(x = index_vals, y = val)) +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_line(aes(group = lat), size = 1, alpha = 0.1) +
  facet_grid(var~test, scales = "free")

