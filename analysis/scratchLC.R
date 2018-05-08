rm(list=ls())
#####
source('loadLibraries.R')
source('analysis_functions.R')


########## Variance Vs. Infected 

exploratory.figures = "../analysis/exploratory.figs/"
data.folder = "../data/tropics/eligible/"


data.list = dir(data.folder)

timeseries = read_outputfiles(dir = data.folder, type = "/out.timeSeries.txt")
viralfitness = read_outputfiles(dir = data.folder, type = "/out.viralFitnessSeries.txt")

timeseries %<>% mutate(propI = totalI/totalN)
trials = unique(timeseries$.id)

plot_varR.v.infections = function(samples, timeseries, viralfitness) {
  sub.time = timeseries %>% filter(.id %in% samples) %>% select(propI, date, .id)
  sub.fitness = viralfitness %>% filter(.id %in% samples) %>% select(day, simDay, varR)
  data.plot = cbind(sub.time, sub.fitness) 
  # checked that time is roughly matched up 
  data.plot %>%
    mutate(scaled.propI = scale(propI),
           scaled.varR = scale(varR)) %>% 
    ggplot(aes(x =scaled.propI, y = scaled.varR)) + geom_point() + 
    facet_wrap(~.id) + labs(x = "Scaled Infected", y = "Scaled Variance in Fitness") -> scatter.plot
  data.plot %>%
    mutate(scaled.propI = scale(propI),
           scaled.varR = scale(varR)) %>%
    select(date, .id, scaled.propI, scaled.varR) %>%
    gather(key = variable, value = value, scaled.propI, scaled.varR, -date, -.id) %>%
    ggplot(aes(x = date, y = value, color = variable)) + 
      geom_line(size = 1.3, alpha = .8) + facet_wrap(~.id) +
      labs(x = "Year", y = "Scaled Value", color = "") + 
      scale_color_brewer(type = "qual", palette = "Dark2", labels = c("Infected", "Variance in Fitness")) +
      theme(legend.position = c(.8,.8)) + 
      theme(strip.text.x = element_blank())  -> dynamics.plot
  
  plot_grid(scatter.plot, dynamics.plot, nrow =2)
}

varR.v.infections1 = plot_varR.v.infections(samples = sample(trials, size = 3), timeseries, viralfitness) 
varR.v.infections2 = plot_varR.v.infections(samples = sample(trials, size = 3), timeseries, viralfitness)
save_plot(plot = varR.v.infections1, filename = "exploratory.figs/varR.v.infections1.pdf", base_height = 8, base_aspect_ratio = 1.8)



############## 
# Variation vs Frequency of Dominant/Transient

tropics.folder = "../data/tropics/eligible/"
trial.dirs = dir(tropics.folder)

# This is a slow step
antigen.frequencies = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.antigenFrequencies.txt"), header = TRUE)) 
names(antigen.frequencies) = trial.dirs

######## Step 2: For each list entry need to calculate max frequency and days above
thres = .2
rel.frequency.thres = .05 #CHANGE 

# Calculating the two criteria for determining success and if successful how long was it there 
antigen.frequencies %>%
  map(function(x) calculate_days_above(x, threshold = .2)) -> days.above

antigen.frequencies %>%
  map(find_max_frequency) %>%
  map(function(x) mutate_at(x, "antigentype", as.character)) -> max.frequencies

# Creating the combined data frame and using the criteria to assign a label 
full.data = map2(max.frequencies, days.above, left_join) %>%
  map(replace_na_zeros) %>% 
  map(determine_success_labels, max.rf =  rel.frequency.thres)
rm(max.frequencies, days.above)

####### Step 3: Subsetting for just transient and yes; and getting rid of zeros 
subset.data = map(full.data, filter_out_loss) 
subset.df = do.call("rbind", subset.data)
subset.df$name = rep(trial.dirs, sapply(subset.data, nrow))

subset.df %>%
  filter(antigentype != 0) -> subset.analyze

subset.analyze %>%
  mutate(key = paste0(antigentype, "_", name)) -> subset.analyze

#antigen.frequencies.df = do.call("rbind", antigen.frequencies)
#antigen.frequencies.df$name = rep(trial.dirs, sapply(antigen.frequencies, nrow))  

antigen.frequencies.df %>%
  mutate(key = paste0(antigentype, "_", name)) %>%
  filter(key %in% subset.analyze$key) %>%
  group_by(antigentype, name) %>%
  arrange(antigentype) -> antigen.frequencies.df.sub

subset.analyze %>%
  select(success, key) %>%
  left_join(antigen.frequencies.df.sub, by = "key") -> antigen.frequencies.df.sub


sample.key = sample(antigen.frequencies.df.sub$key, size = 10)

# Frequency profiles of difference types 
antigen.frequencies.df.sub %>% 
  filter(key %in% sample.key) %>%
  ggplot(aes(x = day/365, y = frequency, color = success)) + geom_line() + facet_wrap(~key, scales = "free_x") +
  scale_x_continuous(breaks = seq(from = 5, to = 25, by=2))
  

antigen.frequencies.df.sub %>%
  group_by(success, antigentype, name) %>%
  summarize(max.freq = mean(max(frequency))) %>%
  ggplot(aes(x=success, y=max.freq)) + geom_boxplot() + 
  labs(x = "Antigenic Fate", y = "Maximum Frequency")  -> max.frequency.plot.bp
  
viralfitness %>%
  select(.id, day, meanR, varR) %>%
  left_join(antigen.frequencies.df.sub, by = c("day", ".id"="name")) -> antigen.frequencies.df.sub
  
viralTrackType = read_outputfiles(dir = data.folder, type = "/out.typeViralFitness.txt")
colnames(viralTrackType)[5:12] = paste0("type_", colnames(viralTrackType)[5:12])

viralTrackType %>%
  mutate(key = paste0(antigenType, "_", .id)) %>%
  filter(key %in% subset.analyze$key) %>%
  select(-.id) %>%
  left_join(antigen.frequencies.df.sub, by = c("day", "key")) %>%
  group_by(key) %>%
  arrange(key) %>%
  mutate(ratio.meanR = exp(type_meanR)/meanR) ->  metrics.by.freq

library(ggridges)

keys = unique(metrics.by.freq$key)
sample.keys = sample(keys, 10)

metrics.by.freq %>%
  group_by(key) %>%
  slice(which.max(frequency)) %>%
  select(day, key) -> peak.day


# This takes a while 
metrics.trial = ddply(metrics.by.freq, .variables = "key", function(x) {
  peak.day %>%
    ungroup() %>%
    filter(key == x$key[1]) %>%
    select(day) -> antigen.peak.day
  x %<>% mutate(phase = ifelse(day <= antigen.peak.day$day, "Growth","Decline"))
  return(x)
})

metrics.trial$phase = factor(metrics.trial$phase, levels = c("Growth", "Decline"))

metrics.trial %>%
  mutate(scale.varR = scale(varR),
      freq.bin = cut(frequency, breaks = seq(from = 0, to = 1,by =.1))) -> metrics.trial
attributes(metrics.trial$scale.varR) = NULL  

metrics.trial %>%
  filter(success == "Est.") %>%
  select(success, phase,simDay, antigenType, varR, ratio.meanR, frequency, freq.bin, scale.varR) %>%
  ggplot(aes(x = ratio.meanR, y = freq.bin, color = phase, fill = phase)) + 
  geom_density_ridges2(scale = 2,,color = "white", rel_min_height = .01, alpha = .8) +
  scale_y_discrete(limits = rev(levels(metrics.trial$freq.bin)),
                   labels = rev(seq(.1,1,.1))) + 
  theme_ridges(center = TRUE) + scale_x_continuous(limits = c(.8, 1.2)) + 
  labs(x = "Relative Effective R", y = "Frequency", fill = "Viral Stage") +
  theme(legend.position = c(.75,.7)) -> ratio.meanR.freq

save_plot(ratio.meanR.freq, filename = "exploratory.figs/ratio.meanR.plot.pdf", base_height = 8)  
  
  
metrics.trial %>%
  filter(phase == "Growth") %>%
  select(success, phase,simDay, antigenType, varR, ratio.meanR, frequency, freq.bin, scale.varR) %>%
  ggplot(aes(x = scale.varR, y = freq.bin, color = success, fill = success)) + 
  geom_density_ridges2(scale = 2,color = "white", rel_min_height = .01, alpha = .8) +
  scale_y_discrete(limits = rev(levels(metrics.trial$freq.bin)),
                   labels = rev(seq(.1,1,.1))) + 
  theme_ridges(center = TRUE) +  scale_x_continuous(limits = c(-2, 4)) +
  labs(x = "Fitness Variance (scaled)", y = "Frequency", fill = "") +
  theme(legend.position = c(.75,.3)) -> varR.v.growth


save_plot(varR.v.growth, filename = "exploratory.figs/varR.v.growth.pdf", base_height = 8)



###### Sample variance
population.1.R = c(1.02, 1.02, 1.02,  1.02, 1.02, 1.1)
mean(population.1.R)
sum((population.1.R - mean(population.1.R))^2)/length(population.1.R)



