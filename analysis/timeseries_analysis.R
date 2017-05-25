# Plotting Timeseries Analysis 


tropics.output.folder = "../data/tropics/"
north.output.folder = "../data/north/"

#create combined meta data for each entry of list  
# Code for comparing north trials that are correct and those that arent' 
north.timeseries = read.outputfiles(north.output.folder, "/out.timeseries.txt")

north.timeseries %>%
  ggplot(aes(x=date, y = totalI*.0025)) +
  geom_line(size = 1.5) +
  facet_wrap(~.id, scales = "free") -> north.timeseries.plot
save_plot("../analysis/exploratory.figs/north.timeseries.plot.pdf", 
          north.timeseries.plot,
          base_height = 8,
          base_aspect_ratio = 1.8)

correct.trial = c("north_2", "north_4", "north_10")

north.track.antigen = read.outputfiles(north.output.folder, "/out.trackAntigenSeries.txt")
north.track.fitness = read.outputfiles(north.output.folder, "/out.viralFitnessSeries.txt")

total.timeseries = rbind(timeseries.all, north.timeseries, idcol = "source")

correct.sample = sample(correct.trials, 1)
incorrect.sample = sample(incorrect.trials, 1)
trial <- display(correct.sample, incorrect.sample)
save_plot(trial, filename = paste0(exploratory.figures, "trial5.pdf"), base_height = 8)



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



