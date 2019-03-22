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
ts_clim_rug <- ggplot(data = sst_ALL, aes(x = t, y = temp)) +
  geom_line(colour = "grey20") +
  geom_line(data = sst_ALL_clim, aes(y = seas),
            linetype = "dashed", colour = "steelblue3") +
  geom_line(data = sst_ALL_clim, linetype = "dotted", colour = "tomato3",
            aes(x = t, y = thresh)) +
  geom_rug(data = sst_ALL_event, sides = "b", colour = "red3", size = 2,
           aes(x = date_peak, y = min(sst_ALL$temp))) +
  facet_wrap(~site, ncol = 1) +
  labs(x = NULL, y = "Temperature (°C)")
ts_clim_rug


# Lolliplots
lolli_WA <- lolli_plot(sst_WA_event)
lolli_WA
lolli_NW_Atl <- lolli_plot(sst_NW_Atl_event)
lolli_NW_Atl
lolli_Med <- lolli_plot(sst_Med_event)
lolli_Med


# The main event
event_WA <- event_line(sst_WA_event)
event_WA
event_NW_Atl <- event_line(sst_NW_Atl_event)
event_NW_Atl
event_Med <- event_line(sst_Med_event)
event_Med


# The clims only
clim_only <- ggplot(data = sst_ALL_clim_only, aes(x = doy)) +
  geom_line(aes(y = seas), colour = "steelblue3") +
  geom_line(aes(y = thresh), colour = "tomato3") +
  facet_wrap(~site, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Temperature (°C)")
clim_only


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
# table_WA <- tableGrob(summary_WA, rows = NULL, theme = tt)


# Figure 2 ----------------------------------------------------------------

# A synthesis of the three tests: length, missing, trended
# The tests are all the same and the three things investigated were
  # climatologies: KS tests
  # event metrics: ANOVA/Tukey
  # category count: chi-squared/residuals
# I envision here some sort of figure that compares the results side by side
# while highlighting how as the degridation of the three tests increases
# how much more rapidly this effects the results than the other tests
# Currently this is being shown with line plots but they aren't popular...

# This figure would likely be the climatology results only


# Figure 3 ----------------------------------------------------------------

# This figure would likely be the event results only


# Figure 4 ----------------------------------------------------------------

# This figure would likely be the cagtegory results only


# Figure 5 ----------------------------------------------------------------

# I imagine this figure being something that supports an argument that follows
# on from one of the above results
# Specifically I think this figure will be the one that shows the average count
# of consecutive missing days depending on the percent of missing data


# Figure 6 ----------------------------------------------------------------

# A sixth figure would be a bit much, but it may be useful to show some sort
# of meta-analysis here
# Such as the relationships between increasing tests and aspects of the time series
# To that end though if such a thing is to be investiagated it would be better
# to run the glbal analysis and show that


