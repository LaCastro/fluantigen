rm(list=ls())
source('analysis_functions.R'); source('plotting_functions.R')

# Plotting Timeseries Analysis 
tropics.folder = "../data/tropics/"
tropics.folder = "../newdata/"


plot_timeseries_id(tropics.timeseries)
########### Tropics #######################################
tropics.timeseries = read_outputfiles(tropics.trial, "/out.timeseries.txt")

tropics.timeseries %>%
  group_by(.id) %>%
  mutate(ratio = totalI/totalS,
         normalized.ratio = (ratio-min(ratio))/(max(ratio)-min(ratio)),
         normalized.I = (totalI-min(totalI))/(max(totalI)-min(totalI)),
         day = date*365) -> tropics.timeseries

################################## Different metrics
trackAntigen <- read_outputfiles(tropics.folder, "/out.trackAntigenSeries.txt")
viralFitness <- read_outputfiles(tropics.folder, "/out.viralFitnessSeries.txt")

variable.sets = c("timeSeries", "trackAntigen", "viralFitness")

a_ply(.data = trials, .margins = 1, function(trial) {
  for(set in variable.sets) {
    plot.full.output = plot_full_output(timeseries = tropics.timeseries, trial = trial, variable.set = set)
    save_plot(filename = paste0("exploratory.figs/full.output.", trial, ".", set, ".pdf"),
              plot = plot.full.output,
              base_height = 8, base_aspect_ratio = 1.6)
  }
})
  










