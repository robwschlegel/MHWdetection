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
  na.omit()

# Linear model of relationship
event_slope_dec_trend_model <- event_slope_dec_trend %>%
  group_by(test, metric) %>%
  # Rename columns so we can use the global_slope() function from the previous section
  # dplyr::rename(val = slope, index_vals = dec_trend)
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

## The proportion of the size of the event in relationship to the change over time
# First filter out only the full 30 year time series events
global_effect_event_30 <- global_effect_event %>%
  filter(test == "length",
         index_vals == 30)

event_slope_prop <- left_join(event_slope_dec_trend, global_effect_event_30, by = c("lon", "lat", "test", "metric")) %>%
  droplevels()

# A regression scatterplot between decadal trend and duration/max.int. and ts length slope
event_slope_prop_scatterplot <- ggplot(data = filter(event_slope_prop,
                                              metric %in% c("duration", "intensity_max")),
                                aes(x = val, y = slope)) +
  geom_point(aes(colour = metric)) +
  geom_smooth(method = "lm") +
  facet_wrap(~metric, scales = "free_y", ncol = 1)
event_slope_prop_scatterplot
