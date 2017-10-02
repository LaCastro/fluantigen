rm(list=ls())
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(broom)


source('analysis_functions.R'); source('plotting_functions.R')

# Plotting Timeseries Analysis 
data.folder = "../data/north/"
fig.folder = "exploratory.figs/north_figures/"

timeseries = read_outputfiles(dir = data.folder, type = "/out.timeSeries.txt")
trials = unique(timeseries$.id)
breaks = seq(from = 1, to = length(unique(timeseries$.id)), length.out = 5)

trials1 = trials[1:25]
trials2 = trials[26:50]
trials3 = trials[51:75]
trials4 = trials[76:100]
trials5 = trials[101:125]
trials6 = trials[126:151]
trials7 = trials[152:175]
trials8 = trials[176:200]
  
  
timeseries %>%
  filter(.id %in% trials8) -> timeseries.sub

twenty.TS.2.8 = plot_timeseries_id(timeseries.sub)
save_plot(twenty.TS.2.8, filename = paste0(fig.folder, "twenty.TS.2.8.pdf"), base_height = 8, base_aspect_ratio = 1.8)

########### Tropics #######################################

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
variable.sets = c("trackAntigen")
trials = unique(tropics.timeseries$.id)
incorrect = c("tropics_trial3", "tropics_trial4", "tropics_trial12")
trials = trials[-which(trials %in% incorrect)]

a_ply(.data = trials, .margins = 1, function(trial) {
  for(set in variable.sets) {
    plot.full.output = plot_full_output(timeseries = tropics.timeseries, trial = trial, variable.set = set)
    save_plot(filename = paste0("exploratory.figs/full.output.", trial, ".", set, ".pdf"),
              plot = plot.full.output,
              base_height = 8, base_aspect_ratio = 1.6)
  }
})
  


################### fixing TMRCA

weeks = c(rep(1:80, each = 9), rep(81,10))

tropics.timeseries$weeks = rep(weeks, 101)

tropics.timeseries %>%
  dplyr::select(.id, date, day, tmrca) %>%
  filter(tmrca < 10 & .id %in% trials1) %>%
  ggplot(aes(x=date, y = tmrca)) + geom_line() + facet_wrap(~.id)

-> tmrca.data
tmrca.data%>%
  summarize(n = n()) -> rows

  summarize(mean = mean(tmrca))
  dplyr::summarize(total.quarters = nrow(.))

  filter(.id %in% trials1) %>%
  ggplot(aes(x=day/365, y = tmrca)) + geom_smooth() + facet_wrap(~.id)
  
  
  ddply(.data = tropics.timeseries, .variables = ".id", function(trial) {
   tmrca = trial$tmrca
   smoothed = rollmean(tmrca,3, align= "right")
   return(smoothed)
  })




