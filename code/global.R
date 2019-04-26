# global.R

# Libraries ---------------------------------------------------------------

source("code/functions.R")

# Set cores
doMC::registerDoMC(cores = 50)


# Global analysis ---------------------------------------------------------

# Set NOAA OISST pathway
OISST_files <- dir(path = "~/data/OISST", full.names = T)

# Run sequentially so that each lon slice can be saved en route
# i <- 6
# for(i in 1:length(OISST_files)){
#
#   # Determine file
#   OISST_slice <- OISST_files[i]
#   lon_row_pad <- str_pad(i, width = 4, pad = "0", side = "left")
#   print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))
#
#   # Calculate tests etc.
#   # NB: This runs the MKE and Eddy masks
#   # system.time(
#   slice_res <- global_analysis(OISST_slice)
#   # ) # ~76 seconds for one
#   save(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rdata"))
#   print(paste0("Finished run on step ",lon_row_pad," at ",Sys.time()))
#
#   # Clear up some RAM
#   gc()
# } # ~ 80 seconds each

# NB: I have not curremtly set this up to work
# The for loop above is still the method for global calculation
# plyr::l_ply(.data = OISST_files, .fun = global_analysis, .parallel = T)


# Unpack results ----------------------------------------------------------

# Unpack and save all of the different global bits
# global_unpack()


# Process results ---------------------------------------------------------

# Calculate the simple linear slopes for the different tests at each pixel
# Of primary interest here is the effect on individual events

# Set cores
doMC::registerDoMC(cores = 50)

# Load data
# load("data/global_effect_event.Rdata")

# Calculate slopes
# global_effect_event_slope <- plyr::ddply(global_effect_event,
                                         # .variables = c("lat"),
                                         # .fun = global_slope, .parallel = T)
# save(global_effect_event_slope, file = "data/global_effect_event_slope.Rdata")


# Global relationships ----------------------------------------------------

# In this section we look at the relationships between certain global variables

# Set cores
doMC::registerDoMC(cores = 50)

# Global decadal trends
load("data/global_dec_trend.Rdata")

# The effect on single events
load("data/global_effect_event.Rdata")

# The slope results for single events
load("data/global_effect_event_slope.Rdata")

# Merge and filter data
event_slope_dec_trend <- left_join(global_effect_event_slope, global_dec_trend, by = c("lon", "lat")) %>%
  filter(test == "length") %>%
  na.omit() %>%
  droplevels()

# Linear model of relationship
event_slope_dec_trend_model <- event_slope_dec_trend %>%
  group_by(test, metric) %>%
  do(model = broom::glance(lm(slope ~ dec_trend, data = .))) %>%
  unnest()

# A regression scatterplot between decadal trend and duration/max.int. and ts length slope
event_dec_scatterplot <- ggplot(data = filter(event_slope_dec_trend,
                                              metric %in% c("duration", "intensity_max")),
                                aes(x = dec_trend, y = slope)) +
  geom_point(aes(colour = metric)) +
  geom_smooth(method = "lm") +
  facet_wrap(~metric, scales = "free_y", ncol = 1)
event_dec_scatterplot

# The spread of the results per year measured
effect_event_spread <- global_effect_event %>%
  filter(test == "length") %>%
  droplevels() %>%
  spread(key = index_vals, value = val) %>%
  mutate(val_spread = `30`-`10`)

# Plot showing value spread between 30 and 10 years
duration_spread_plot <- ggplot(filter(effect_event_spread, metric == "intensity_max"),
                                      aes(x = lon, y = lat)) +
  geom_raster(aes(fill = val_spread)) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
  scale_fill_gradient2(low = "blue", high = "red") +
  coord_equal(expand = F) +
  # labs(fill = "Linear trend (°C/dec)") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))
duration_spread_plot

# Join in the event metrics detected from full 30 year time series
event_slope_prop <- left_join(event_slope_dec_trend, effect_event_spread,
                              by = c("lon", "lat", "test", "metric")) %>%
  na.omit() %>%
  mutate(slope_prop_10 = round(slope/`10`, 3),
         slope_prop_20 = round(slope/`20`, 3),
         slope_prop_30 = round(slope/`30`, 3))


# Get just slope ten and change it for easy plotting
event_slope_prop_10 <- event_slope_prop %>%
  dplyr::select(lat:metric, slope_prop_10) %>%
  dplyr::rename(slope = slope_prop_10) %>%
  mutate(metric = "slope_prop_10")

event_slope_10_quantiles <- quantile(event_slope_prop_10$slope, na.rm = T,
                                     probs = c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1.0))

# Global map showing patterns of the slope as a proportion of the overall change
int_max_spread_prop_plot <- ggplot(filter(event_slope_prop, metric == "intensity_max"),
                               aes(x = lon, y = lat)) +
  geom_raster(aes(fill = slope_prop_10)) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
  scale_fill_gradient2(low = "blue", high = "red") +
  coord_equal(expand = F) +
  labs(fill = "Proportion change\yearLinear trend (°C/dec)") +
  scale_fill_gradient2(low = col_split[1], high = col_split[2],
                       breaks = c(as.numeric(slope_quantiles[2:6]))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))
int_max_spread_prop_plot

# Boxplot showing spread of proportion values


# Linear model of relationship between full event metric and rate of change from shortening
event_slope_prop_model <- event_slope_prop %>%
  group_by(test, metric) %>%
  filter(metric %in% c("duration", "intensity_max")) %>%
  do(model = broom::glance(lm(slope ~ dec_trend, data = .))) %>%
  unnest()

# A regression scatterplot between decadal trend and duration/max.int. and ts length slope
event_slope_prop_scatterplot <- ggplot(data = filter(event_slope_prop,
                                              metric %in% c("duration", "intensity_max")),
                                aes(x = val, y = slope)) +
  geom_point(aes(colour = metric)) +
  geom_smooth(method = "lm") +
  facet_wrap(~metric, scales = "free", ncol = 1)
event_slope_prop_scatterplot
