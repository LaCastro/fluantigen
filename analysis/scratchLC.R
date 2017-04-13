#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)


source('analysis_functions.R')
source('plotting_functions.R')

### Timeseries
output.folder = "~/Dropbox/Projects/mutantigen/north_runs/"

## Will need this later for reading in multiple files 
variation.list <- lapply(variation.list, function(.file) {
  file.list = list.files(.file)
  out.timeseries <- read.table(paste0(.file,"/out.trackAntigenSeries.txt"), header = TRUE)
  out.timeseries
})

variation.timeseries = rbindlist(all.timeseries, idcol = TRUE)
variation.timeseries$.id = as.factor(variation.timeseries$.id)



##### Population
analysis.dir <- "~/Dropbox/Projects/mutantigen/population_test/"
population.list = list.dirs(analysis.dir, full.names = FALSE, recursive = FALSE)

population.antigenFrequencies <- lapply(population.list, function(.file) {
  out.timeseries <- read.table(paste0(analysis.dir,.file, "/out.antigenFrequencies.txt"), header = TRUE)
  out.timeseries
})

names(population.antigenFrequencies) = population.list
population.antigenFrequencies = rbindlist(population.antigenFrequencies, idcol = TRUE)


################# Antigen Series
timeseries  = read.table(paste0(output.folder, "/out.timeseries.txt"), header = TRUE)
desired.metrics = c("northI", "antigenicDiversity", "northTmrca", "diversity")


timeseries %>%
  gather(key = metric, value = value, - date) %>%
  filter(metric %in% desired.metrics ) %>%
  ggplot(aes(x = date, y = value)) + geom_line() +
  facet_wrap(~metric, scales = "free")

antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)

successful.types = find.successful.types(frequencies = antigen.frequencies, threshold = .1)
myColors = set.my.colors(length(successful.types))


antigen.frequencies %>%
  filter(antigentype %in% successful.types) %>%
  ggplot(aes(x = day, y = frequency, color = antigentype)) + 
  geom_line() + guides(col = FALSE) + 
  scale_color_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") -> frequency.short


frequency.plot <- plot.successful.frequency(antigen.frequencies, successful.types)
infecteds.plot <- plot.successful.infections(antigen.frequencies, successful.types)


################################ Plotting successful versus non-successful
antigen.specific <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")
pop.dynamics <- c("S", "I")
viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")


meta.data %>%
  gather(key = metric, value = value, -postAntigen, -success) -> meta.data.long

meta.data.long$value = as.numeric(meta.data.long$value)

meta.data.long %>%
  filter(metric %in% pop.dynamics) %>%
  ggplot(aes(success, value)) + geom_boxplot() +
  facet_wrap(~metric, scales = "free")

meta.data.long %>%
  filter(metric %in% antigen.specific) %>%
  ggplot(aes(success, value)) + geom_boxplot() + 
  facet_wrap(~metric, scales = "free")

meta.data.long %>%
  filter(metric %in% viral.fitness.metrics) %>%
  ggplot(aes(success, value)) + geom_boxplot() +
  facet_wrap(~metric, scales = "free")



####################### Histograms of each, divide it by success
meta.data.long %>%
  filter(metric %in% antigen.specific) %>%
  ggplot(aes(x=value, color = success, fill  = success)) + 
  geom_histogram(bins=5, position = "dodge") + 
  facet_wrap(~metric, scales = "free")


meta.data.long %>%
  filter(metric %in% viral.fitness.metrics) %>%
  ggplot(aes(x=value, color = success, fill  = success)) + 
  geom_histogram(bins=5, position = "dodge") + 
  facet_wrap(~metric, scales = "free")


meta.data.long %>%
  filter(metric %in% pop.dynamics) %>%
  ggplot(aes(x=value, color = success, fill  = success)) + 
  geom_histogram(bins=5, position = "dodge") + 
  facet_wrap(~metric, scales = "free")



#### Read multiple runs of similar type 
console.data = read.console.files(output.folder)


#create combined meta data for each entry of list 
trial.meta <- create.meta.all(dir = output.folder)
timeseries.all = read.outputfiles(output.folder, "/out.timeseries.txt")


timeseries.all %>% 
  gather(key = metric, value = value, -.id, -date) %>%
  filter(metric == "totalI") %>%
  ggplot(aes(x = date, y = value, group = .id, color = .id)) + 
  geom_line(size = 1.3)

desired.metrics = c("diversity", "tmrca", "antigenicDiversity")
timeseries.all %>%
  gather(key = metric, value = value, -.id, -date) %>%
  filter(metric %in% desired.metrics) %>%
  ggplot(aes(x = date, y = value, group = .id, color = .id)) +
  geom_smooth() +
  facet_wrap(~metric, scales = "free", ncol = 1)
  


########## Antigen Frequencies START HERE 
source('analysis_functions.R')
output.folder = "~/Dropbox/Projects/mutantigen/north_runs/"

#Read and combine files 
north.data = create.meta.data.all(dir = output.folder)
antigen.freq.all <- read.outputfiles(output.folder, "/out.antigenFrequencies.txt")

north.data %>%
  group_by(.id) %>%
  filter(success == "yes") %>%
  select(.id, postAntigen) -> success.types

### For each successtype, go into antigen frequecny all, filter and full out

##### Read in the antigen frequencies and then store those in a list 

antigen.freq.success.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  antigen.freq.all %>%
    filter(.id == sim$.id[1]) %>%
    filter(antigentype %in% successful.types) -> antigen.freq.sim
  return(antigen.freq.sim)
})

#Determine the maximum number of colors going to need 
antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = length(unique(antigentype))) -> num.transitions

max.color = sum(num.transitions$num.transitions)
myColors = set.my.colors(max.color)

## Need to first go in and fill all the missing values and then combine
ant.freq.success.l = ddply(.data = antigen.freq.success.df, .variables = ".id", function(sim) {
  sim %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2, - infected)  -> antigen.freq.long
  return(antigen.freq.long)
})

ant.freq.success.l$antigentype = as.factor(ant.freq.success.l$antigentype)

## Need to first go in and fill all the missing values and then combine
# Plot Frequency First 
ant.freq.success.l %>%
  mutate(year = day/365) %>%
  filter(frequency > 0) %>%
  ggplot(aes(x = year, y = frequency, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id) +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  guides(col = FALSE) + guides(fill = FALSE) -> freq.plot
freq.plot
save_plot(filename = paste0(output.folder, "exploratory_figures/freq.plot.pdf"), freq.plot,
          base_height = 8, base_aspect_ratio = 1.5)
  
ant.freq.success.l %>%
  mutate(year = day/365) %>%
  mutate( prevalence = infected*frequency) %>%
  filter(prevalence > 0) %>%
  ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  guides(col = FALSE) + guides(fill = FALSE) -> prev.plot

prev.plot
save_plot(filename = paste0(output.folder, "exploratory_figures/prev.plot.pdf"), prev.plot,
          base_height = 8, base_aspect_ratio = 1.5)

##########################
# Histograms of the metrics
 
antigen.specific <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")
pop.dynamics <- c("S", "I")
viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")

north.data %>%
  gather(key = metric, value = value, -postAntigen, -success, -.id) -> north.data.long
north.data.long$value = as.numeric(north.data.long$value)

north.data.long %>%
  filter(metric %in% pop.dynamics) %>%
  filter(metric == "I") %>%
  ggplot(aes(success, value)) + geom_boxplot() +
  facet_wrap(~.id, scales = "free")

# Going to group,

north.data.long %>%
  filter(metric %in% antigen.specific) %>%
  ggplot(aes(success,value)) + 
  geom_violin() + geom_jitter(height = 0, width = 0.1) +
  facet_grid(metric~.id, scales = "free")


north.data.long %>%
  filter(metric %in% viral.fitness.metrics) %>%
  ggplot(aes(success,value)) + 
  geom_violin() + geom_jitter(height = 0, width = 0.1) +
  facet_grid(metric~.id, scales = "free")
