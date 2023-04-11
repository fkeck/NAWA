
source("R/load_data.R")


project <- "Trend"
rhine_catchment_only <- TRUE
res_ML <- list()

fraction_type <- "Natural"
status <- "status_q50"
res_ML$Natural <- list()
source("R/RFC_Trend_IBCH_OTU97_Natural.R")
source("R/stats_paper.R")
source("R/plots_paper.R")


# Workspace saved in /home/ecoadmin/Documents/temp/NAWA_ML_frozen.RData
# load("~/Documents/temp/NAWA_ML_frozen.RData")
