# Plotting Timeseries Analysis 



tropics.folder = "../data/tropics/"
north.folder = "../data/north/"

#create combined meta data for each entry of list  
# Code for comparing north trials that are correct and those that arent' 

########### Tropics #######################################
tropics.timeseries = read_outputfiles(tropics.folder, "/out.timeseries.txt")

tropics.timeseries %>%
  group_by(.id) %>%
  mutate(ratio = totalI/totalS,
         normalized.ratio = (ratio-min(ratio))/(max(ratio)-min(ratio)),
         normalized.I = (totalI-min(totalI))/(max(totalI)-min(totalI))) -> tropics.timeseries

tropics.timeseries %>%
  select(.id, date, normalized.ratio, normalized.I) %>%
  gather(key = metric, value = value, 3:4) %>%
  ggplot(aes(x=date, y = value, color = metric)) + 
  geom_line(size = 1.5)+
  facet_wrap(~.id)

tropics.timeseries %>%
  ggplot(aes(x=date, y = totalI*.0025)) +
  geom_line(size = 1.5) +
  facet_wrap(~.id, scales = "free") -> tropics.timeseries.plot

tropics.trials = unique(tropics.timeseries$.id)
tropics.correct.trials = tropics.trials[-which(tropics.trials == "tropics_11")]

save_plot("../analysis/exploratory.figs/tropics.timeseries.plot.pdf", 
          tropics.timeseries.plot,
          base_height = 8,
          base_aspect_ratio = 1.8)

############### North #######################################
north.timeseries = read_outputfiles(north.folder, "/out.timeseries.txt")

north.trials = as.data.frame(unique(north.timeseries$.id)); colnames(north.trials) = "trial"


north.correct.trial = c("north_2", 
                        "north_4", 
                        "north_10", 
                        "north_11", 
                        "north_12",
                        "north_20",
                        "north_22")

north.trials %>%
  filter(!(trial %in% north.correct.trial)) -> north.incorrect.trial

north.track.antigen = read_outputfiles(north.output.folder, "/out.trackAntigenSeries.txt")
north.track.fitness = read_outputfiles(north.output.folder, "/out.viralFitnessSeries.txt")



total.timeseries = rbind(timeseries.all, north.timeseries, idcol = "source")
correct.sample = sample(north.correct.trial, 1)
incorrect.sample = as.character(sample(north.incorrect.trial$trial, 1))

trial <- display(correct.sample, incorrect.sample)
save_plot(trial, filename = "exploratory.figs/trial6.pdf", base_height = 8)



timeseries.all %>%
  group_by(.id) %>%
  mutate(years = (seq(from = 10, to = 3650, by = 10))/365) %>%
  ungroup() %>%
  group_by(years) %>%
  summarize(mean.diversity = mean(diversity),
            mean.antiDiversity = mean(antigenicDiversity), 
            mean.IS = mean(totalI/totalS),
            mean.tmrca = mean(tmrca),
            cases.per100k = mean(totalI)/100000) %>%
  gather(key = metric, value = value, -years) %>%
  mutate(region = "tropics") -> tropics.timeseries.long


north.timeseries %>%
  group_by(.id) %>%
  filter(.id %in% desired.trials) %>%
  mutate(years = (seq(from = 10, to = 3650, by = 10))/365) %>%
  ungroup() %>%
  group_by(years) %>%
  summarize(mean.diversity = mean(diversity),
            mean.antiDiversity = mean(antigenicDiversity), 
            mean.IS = mean(totalI/totalS),
            mean.tmrca = mean(tmrca),
            cases.per100k = mean(totalI)/100000) %>%
  gather(key = metric, value = value, -years) %>%
  mutate(region = "north") -> north.timeseries.long






timeseries.long = rbind(tropics.timeseries.long, north.timeseries.long)

timeseries.long %>%
  ggplot(aes(x = years, y = value, color = region)) + 
  geom_line()  + facet_wrap(~metric, scales = "free")



