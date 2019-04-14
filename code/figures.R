# code/figures.R
# This script houses the code used for the final publication figures

source("code/functions.R")


# Figure 1 ----------------------------------------------------------------

# A multi-paneled figure showing a summary of information for the three reference time series
# These could include:
  # Time series with rug plot showing events
  # Lolliplots
  # The clims only
  # The main event
  # Map/location point
  # Summary/stats table

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

# Focus events
MHW_focus <- sst_ALL_event %>%
  group_by(site) %>%
  filter(intensity_cumulative == max(intensity_cumulative))
# Interestingly this grabs a different MHW for the Med
# so we need to more manually grab that one
focus_Med <- sst_ALL_event %>%
  filter(site == "Med", date_end <= "2004-01-01") %>%
  filter(intensity_cumulative == max(intensity_cumulative)) %>%
  left_join(sst_ALL_coords, by = "site")
focus_WA <- MHW_focus[3,] %>%
  left_join(sst_ALL_coords, by = "site")
focus_NW_Atl <- MHW_focus[2,] %>%
  left_join(sst_ALL_coords, by = "site")


# Time series with rug plot showing events
ts_clim_rug_WA <- ts_clim_rug("WA")
ts_clim_rug_WA
ts_clim_rug_NW_Atl <- ts_clim_rug("NW_Atl")
ts_clim_rug_NW_Atl
ts_clim_rug_Med <- ts_clim_rug("Med")
ts_clim_rug_Med


# Lolliplots
lolli_WA <- lolli_plot(sst_WA_event)
lolli_WA
lolli_NW_Atl <- lolli_plot(sst_NW_Atl_event)
lolli_NW_Atl
lolli_Med <- lolli_plot(sst_Med_event)
lolli_Med


# The main event
event_WA <- event_line(sst_WA_event, category = T)
event_WA
event_NW_Atl <- event_line(sst_NW_Atl_event, category = T)
event_NW_Atl
event_Med <- event_line(sst_Med_event, category = T)
event_Med


# The clims only
clim_WA <- clim_line("WA")
clim_WA
clim_NW_Atl <- clim_line("NW_Atl")
clim_NW_Atl
clim_Med <- clim_line("Med")
clim_Med


# Map/location point
map_WA <- map_point(focus_WA)
map_WA
map_NW_Atl <- map_point(focus_NW_Atl)
map_NW_Atl
map_Med <- map_point(focus_Med)
map_Med


# Summary/stats table
# Ideas for table:
  # Count of events
  # Mean/ max duration
  # mean/ max intensity
  # Seasonal range
  # Threshold range
  # Category count
table_summary_WA <- table_summary("WA")
# table_summary_WA
table_summary_NW_Atl <- table_summary("NW_Atl")
# table_summary_NW_Atl
table_summary_Med <- table_summary("Med")
# table_summary_Med


## Stitch all of the figures together
# Time series with rug plot showing events
# Lolliplots
# The clims only
# The main event
# Map/location point
# Summary/stats table
stitch_plot_WA <- ggpubr::ggarrange(ts_clim_rug_WA, clim_WA, lolli_WA,
                                    event_WA, map_WA, table_summary_WA, labels = "AUTO")
stitch_plot_WA
# ggsave(plot = stitch_plot_WA, filename = "output/stitch_plot_WA.pdf", height = 8, width = 14)
stitch_plot_NW_Atl <- ggpubr::ggarrange(ts_clim_rug_NW_Atl, clim_NW_Atl,
                                        lolli_NW_Atl, event_NW_Atl, map_NW_Atl,
                                        table_summary_NW_Atl, labels = "AUTO")
stitch_plot_NW_Atl
stitch_plot_Med <- ggpubr::ggarrange(ts_clim_rug_Med, clim_Med, lolli_Med,
                                     event_Med, map_Med, table_summary_Med, labels = "AUTO")
stitch_plot_Med

## A more paired down figure
stitch_sub_plot_WA <- ggpubr::ggarrange(event_WA, map_WA, labels = "AUTO")
stitch_sub_plot_WA
ggsave(plot = stitch_sub_plot_WA,
       filename = "output/stitch_sub_plot_WA.pdf",
       height = 4, width = 8)


# Figure 2 ----------------------------------------------------------------

# A synthesis of the three tests: length, missing, trended
# The tests are all the same and the three things investigated were
  # climatologies: KS tests
load("data/sst_ALL_KS_clim.Rdata")
  # event metrics: ANOVA/Tukey OR KS tests
load("data/sst_ALL_KS_event.Rdata")
load("data/sst_ALL_aov_tukey.Rdata")
  # category count: fisher test
load("data/sst_ALL_KS_cat.Rdata")
load("data/sst_ALL_fisher.Rdata")
# I envision here some sort of figure that compares the results side by side
# while highlighting how as the degridation of the three tests increases
# how much more rapidly this effects the results than the other tests
# Currently this is being shown with line plots but they aren't popular...

# Climatology results
sst_ALL_KS_clim_long <- sst_ALL_KS_clim %>%
  gather(key = "metric", value = "p.value", -test, -site, -rep, -index_vals) %>%
  group_by(test, site, metric, index_vals) %>%
  summarise(p.value.mean = mean(p.value))# %>%
  # mutate(index_col = paste(site, test, sep = "-"))
# sst_ALL_KS_clim_long_sig <- sst_ALL_KS_clim_long %>%
#   filter(p.value.mean <= 0.05)

# Event results
sst_ALL_KS_event_long <- sst_ALL_KS_event %>%
  gather(key = "metric", value = "p.value", -test, -site, -rep, -index_vals) %>%
  group_by(test, site, metric, index_vals) %>%
  summarise(p.value.mean = mean(p.value)) #%>%
# mutate(index_col = paste(test, site, metric, sep = "-"))
# sst_ALL_KS_event_long_sig <- sst_ALL_KS_event_long %>%
#   filter(p.value.mean <= 0.05)

# Cagtegory results
sst_ALL_KS_cat_long <- sst_ALL_KS_cat %>%
  gather(key = "metric", value = "p.value", -test, -site, -rep, -index_vals) %>%
  group_by(test, site, metric, index_vals) %>%
  summarise(p.value.mean = mean(p.value)) #%>%
# mutate(index_col = paste(test, site, metric, sep = "-"))
# sst_ALL_KS_cat_long_sig <- sst_ALL_KS_cat_long %>%
#   filter(p.value.mean <= 0.05)

# Combine for plotting
sst_ALL_plot_long <- rbind(sst_ALL_KS_clim_long, sst_ALL_KS_event_long, sst_ALL_KS_cat_long) %>%
  filter(!metric %in% c("intensity_mean", "intensity_cumulative")) %>%
  ungroup() %>%
  mutate(metric = factor(metric, levels = c("seas", "thresh",
                                            "duration", "intensity_max",
                                            "p_moderate", "p_strong", "p_severe", "p_extreme")),
         test = case_when(test == "length" ~ "length (years)",
                          test == "missing" ~ "missing data (proportion)" ,
                          test == "trended" ~ "added trend (°C/dec)"),
         test = as.factor(test),
         test = factor(test, levels = levels(test)[c(2,3,1)]))

# Plot them all together
fig_2 <- ggplot(sst_ALL_plot_long, aes(x = index_vals, y = p.value.mean, colour = metric)) +
  geom_line() +
  geom_point() +
  geom_point(data = filter(sst_ALL_plot_long, p.value.mean <= 0.05), shape = 15, colour = "red", size = 2) +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "dashed") +
  scale_colour_manual(name = "Metric",
                      values = c("skyblue", "navy",
                                 "springgreen", "forestgreen",
                                 "#ffc866", "#ff6900", "#9e0000", "#2d0000"),
                      labels = c("Seasonal climatology", "Threshold climatology",
                                 "Duration (days)", "Max. intensity (°C)",
                                 "Prop. moderate", "Prop. strong", "Prop. severe", "Prop. extreme")) +
  # geom_errorbarh(aes(xmin = year_long, xmax = year_short)) +
  guides(colour = guide_legend(override.aes = list(shape = 15, linetype = NA, size = 3))) +
  labs(y = "Mean p-value (n = 100)", x = NULL) +
  facet_grid(site~test, scales = "free_x", switch = "x") +
  theme(legend.position = "bottom")
fig_2
ggsave(plot = fig_2, filename = "LaTeX/fig_2.pdf", height = 8, width = 12)


# Figure 3 ----------------------------------------------------------------

# A table here showing the R2 values from Figure 2 would be good
sst_ALL_R2_long <- sst_ALL_plot_long %>%
  group_by(test, site, metric) %>%
  filter(index_vals <= 30) %>%
  nest() %>%
  mutate(R2 = map(data, lm_p_R2)) %>%
  select(-data) %>%
  unnest()


# Figure 4 ----------------------------------------------------------------

# I imagine this figure being something that supports an argument that follows
# on from one of the above results
# Specifically I think this figure will be the one that shows the average count
# of consecutive missing days depending on the percent of missing data


# Figure 5 ----------------------------------------------------------------


# Figure 6 ----------------------------------------------------------------

# A sixth figure would be a bit much, but it may be useful to show some sort
# of meta-analysis here
# Such as the relationships between increasing tests and aspects of the time series
# To that end though if such a thing is to be investiagated it would be better
# to run the global analysis and show that



###
# Perhaps some basic representations of what is going on would be best
# Text book like diagrams with arrows showing the general trends in how
# the data are affected when the three variables are tweeked

