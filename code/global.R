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
for(i in 1:length(OISST_files)){

  # Determine file
  OISST_slice <- OISST_files[i]
  lon_row_pad <- str_pad(i, width = 4, pad = "0", side = "left")
  print(paste0("Began run on step ",lon_row_pad," at ",Sys.time()))

  # Calculate tests etc.
  # NB: This runs the MKE and Eddy masks
  # system.time(
  slice_res <- global_analysis(OISST_slice)
  # ) # ~76 seconds for one
  save(slice_res, file = paste0("data/global/slice_",lon_row_pad,".Rdata"))
  print(paste0("Finished run on step ",lon_row_pad," at ",Sys.time()))

  # Clear up some RAM
  gc()
} # ~ 80 seconds each
