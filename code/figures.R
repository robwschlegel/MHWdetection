# code/figures.R
# This script houses the code used for the final publication figures

source("code/functions.R")


# Figure 1 ----------------------------------------------------------------

# Show the 3 focus MHWs and the category layers

# Load base data
load("data/sst_ALL.Rdata")

sst_WA_event <- detect_event(ts2clm(filter(sst_ALL, site == "WA"),
                                    climatologyPeriod = c("1982-01-01", "2011-12-31")))
sst_NW_Atl_event <- detect_event(ts2clm(filter(sst_ALL, site == "NW_Atl"),
                                    climatologyPeriod = c("1982-01-01", "2011-12-31")))
sst_Med_event <- detect_event(ts2clm(filter(sst_ALL, site == "Med"),
                                    climatologyPeriod = c("1982-01-01", "2011-12-31")))

# Calculate base results
sst_ALL_res <- sst_ALL %>%
  group_by(site) %>%
  nest() %>%
  mutate(clims = map(data, ts2clm, maxPadLength = 1,
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

# Manually grab each event
focus_Med <- sst_ALL_event %>%
  filter(site == "Med", date_end <= "2004-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative)) %>%
  left_join(sst_ALL_coords, by = "site")
focus_WA <- sst_ALL_event %>%
  filter(site == "WA", date_end <= "2013-01-01") %>%
  filter(intensity_max == max(intensity_max)) %>%
  left_join(sst_ALL_coords, by = "site")
focus_NW_Atl <- sst_ALL_event %>%
  filter(site == "NW_Atl", date_end <= "2013-01-01") %>%
  filter(intensity_max == max(intensity_max)) %>%
  left_join(sst_ALL_coords, by = "site")

# Create the three panelled events
focus_WA_plot <- fig_1_plot(sst_WA_event, focus_WA$date_peak, 183)
focus_NW_Atl_plot <- fig_1_plot(sst_NW_Atl_event, focus_NW_Atl$date_peak, 183)
focus_Med_plot <- fig_1_plot(sst_Med_event, focus_Med$date_peak, 183)
fig_1 <- ggarrange(focus_WA_plot, focus_NW_Atl_plot, focus_Med_plot,
                   ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(fig_1, filename = "LaTeX/fig_1.pdf", width = 8, height = 8)
ggsave(fig_1, filename = "LaTeX/fig_1.png", width = 8, height = 8)


# Create flattened figure for poster
focus_WA_flat <- fig_1_plot(sst_WA_event, focus_WA$date_peak, 183, y_label = NULL)
focus_NW_Atl_flat <- fig_1_plot(sst_NW_Atl_event, focus_NW_Atl$date_peak, 183, y_label = "Temp. (°C)")
focus_Med_flat <- fig_1_plot(sst_Med_event, focus_Med$date_peak, 183, y_label = NULL)
fig_1_flat <- ggarrange(focus_WA_flat, focus_NW_Atl_flat, focus_Med_flat, align = "v",
                        ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(fig_1_flat, filename = "LaTeX/fig_1_flat.png", width = 8, height = 5)


# Figure 2 ----------------------------------------------------------------

# The results from the 100 re-sampled length tests

# Sub-optimal data
load("data/sst_ALL_KS_clim.Rdata")
load("data/sst_ALL_KS_event.Rdata")
load("data/sst_ALL_KS_cat.Rdata")
# Fixed data
load("data/sst_ALL_KS_clim_fix.Rdata")
load("data/sst_ALL_KS_event_fix.Rdata")
load("data/sst_ALL_KS_cat_fix.Rdata")

# Climatology results
sst_ALL_KS_clim_long <- KS_long(sst_ALL_KS_clim)
sst_ALL_KS_clim_fix_long <- KS_long(sst_ALL_KS_clim_fix)

# Event results
sst_ALL_KS_event_long <- KS_long(sst_ALL_KS_event)
sst_ALL_KS_event_fix_long <- KS_long(sst_ALL_KS_event_fix)

# Cagtegory results
sst_ALL_KS_cat_long <- KS_long(sst_ALL_KS_cat)
sst_ALL_KS_cat_fix_long <- KS_long(sst_ALL_KS_cat_fix)

# Combine for plotting
sst_ALL_plot_long <- rbind(sst_ALL_KS_clim_long, sst_ALL_KS_event_long,
                           sst_ALL_KS_cat_long, sst_ALL_KS_clim_fix_long,
                           sst_ALL_KS_event_fix_long, sst_ALL_KS_cat_fix_long) %>%
  filter(!metric %in% c("intensity_mean", "intensity_cumulative"),
         index_vals <= 0.5 | index_vals >= 10) %>%
  mutate(metric = factor(metric, levels = c("seas", "thresh",
                                            "duration", "intensity_max",
                                            "p_moderate", "p_strong", "p_severe", "p_extreme")))

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

# Create individual panels
# fig_2_WA <- fig_line_plot("WA", "length")
# fig_2_WA
# fig_2_NW_Atl <- fig_line_plot("NW_Atl", "length")
# fig_2_NW_Atl
# fig_2_Med <- fig_line_plot("Med", "length")
# fig_2_Med

# Combine and save
# fig_2 <- ggarrange(fig_2_WA, fig_2_NW_Atl, fig_2_Med,
                   # ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
fig_2 <- fig_line_plot("length")
ggsave(plot = fig_2, filename = "LaTeX/fig_2.pdf", height = 4, width = 8)
ggsave(plot = fig_2, filename = "LaTeX/fig_2.png", height = 4, width = 8)


# Figure 3 ----------------------------------------------------------------

# The effects of the three sub-optimal tests on the focus MHW metrics

# This is currently made in the workflow.R script


# Figure 4 ----------------------------------------------------------------

# The global effect of length on max. intensity

# This figure is currently made in the global.R script


# Figure 5 ----------------------------------------------------------------

# The global effect of length on duration

# This figure is currently made in the global.R script


# Figure 6 ----------------------------------------------------------------

# The results from the 100 re-sampled missing data tests

# Create individual panels
fig_6_WA <- fig_line_plot("WA", "missing")
# fig_6_WA
fig_6_NW_Atl <- fig_line_plot("NW_Atl", "missing")
# fig_6_NW_Atl
fig_6_Med <- fig_line_plot("Med", "missing")
# fig_6_Med

# Combine and save
fig_6 <- ggpubr::ggarrange(fig_6_WA, fig_6_NW_Atl, fig_6_Med,
                   ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(plot = fig_6, filename = "LaTeX/fig_6.pdf", height = 8, width = 8)
ggsave(plot = fig_6, filename = "LaTeX/fig_6.png", height = 8, width = 8)


# Figure 7 ----------------------------------------------------------------

# The fix for missing data

# Create individual panels
fig_7_WA <- fig_line_plot("WA", "missing_fix")
# fig_7_WA
fig_7_NW_Atl <- fig_line_plot("NW_Atl", "missing_fix")
# fig_7_NW_Atl
fig_7_Med <- fig_line_plot("Med", "missing_fix")
# fig_7_Med

# Combine and save
fig_7 <- ggpubr::ggarrange(fig_7_WA, fig_7_NW_Atl, fig_7_Med,
                           ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(plot = fig_7, filename = "LaTeX/fig_7.pdf", height = 8, width = 8)
ggsave(plot = fig_7, filename = "LaTeX/fig_7.png", height = 8, width = 8)


# Figure 8 ----------------------------------------------------------------

# The results from the 100 re-sampled added trend tests

# Create individual panels
fig_8_WA <- fig_line_plot("WA", "trended")
# fig_8_WA
fig_8_NW_Atl <- fig_line_plot("NW_Atl", "trended")
# fig_8_NW_Atl
fig_8_Med <- fig_line_plot("Med", "trended")
# fig_8_Med

# Combine and save
fig_8 <- ggpubr::ggarrange(fig_8_WA, fig_8_NW_Atl, fig_8_Med,
                   ncol = 1, nrow = 3, labels = "AUTO", common.legend = T, legend = "top")
ggsave(plot = fig_8, filename = "LaTeX/fig_8.pdf", height = 8, width = 8)
ggsave(plot = fig_8, filename = "LaTeX/fig_8.png", height = 8, width = 8)


# Figure 10 ---------------------------------------------------------------

# The effect of added trends on the focus MHWs











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
