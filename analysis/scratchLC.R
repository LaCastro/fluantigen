#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)

### Timeseries
analysis.dir <- "~/Dropbox/Projects/mutantigen/04-04-2017_09-20-07/"
analysis.dir <- "~/Documents/projects/fluantigen/04-04-2017_03-13-02/"

## Will need this later for reading in multiple files 
variation.list <- lapply(variation.list, function(.file) {
  file.list = list.files(.file)
  out.timeseries <- read.table(paste0(.file,"/out.trackAntigenSeries.txt"), header = TRUE)
  out.timeseries
})

variation.timeseries = rbindlist(all.timeseries, idcol = TRUE)
variation.timeseries$.id = as.factor(variation.timeseries$.id)

myColors <- colorRampPalette(brewer.pal(8, "Dark2"))(10)
### Visualzing the frequency of the wildtype


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
timeseries  = read.table(paste0(analysis.dir, "/out.timeseries.txt"), header = TRUE)
desired.metrics = c("northI", "antigenicDiversity", "northTmrca", "diversity")

timeseries %>%
  gather(key = metric, value = value, - date) %>%
  filter(metric %in% desired.metrics ) %>%
  ggplot(aes(x = date, y = value)) + geom_line() +
  facet_wrap(~metric, scales = "free")

antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)



antigen.frequencies %>%
  filter(antigentype %in% successful.types) %>%
  ggplot(aes(x = day, y = frequency, color = antigentype)) + 
  geom_line() + guides(col = FALSE) + 
  scale_color_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") -> frequency.short


frequency.plot <- plot.successful.frequency(antigen.frequencies, successful.types)
infecteds.plot <- plot.successful.infections(antigen.frequencies, successful.types)


####### Antigenic Change
# one line graph, number of times they occur
# histogram, the size
library(readr)
antigenic.mutations <- read_delim("~/Documents/projects/fluantigen/03-27-2017_09-26/out.console.txt", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
antigenic.mutations$distance = as.numeric(antigenic.mutations$distance)

# Plotting timeseries of antigenic mutations 
antigenic.mutations %>% ggplot(aes(x = day, y = distance)) + 
  geom_point() + facet_wrap(~oriAntigenType)

antigenic.mutations %>% ggplot(aes(x = distance)) + geom_histogram() + geom_vline(xintercept = .012)

antigenic.mutations %>% filter(oriAntigenType == "0") %>%
  ggplot(aes(x = distance)) + geom_histogram() + geom_vline(xintercept = .012)

antigenic.mutations %>% group_by(day) %>%
  summarise(num.mutations = length(day)) %>%
  ggplot(aes(x = day, y = num.mutations)) + geom_line()


##### Plotting frequencies of different antigenic 


### Mean Viral Fitness Series 

viral.fitness = read.table(paste0(analysis.dir, "out.viralFitnessSeries.txt"), header = TRUE)
viral.fitness = viral.fitness[!duplicated(viral.fitness),]


metric.mean = c("meanR", "varR")
metric.beta = c("meanBeta", "varBeta")
metric.sigma = c("meanSigma", "varSigma")

viral.fitness %>%
  gather(key = metric, value = value, -date) %>%
  mutate(value.exp = exp(value)) %>%
  filter(metric %in% metric.mean) %>%
  select(date, metric, value.exp) %>%
  spread(key = metric, value = value.exp) %>%
  ggplot(aes(x = date, y = meanR)) + geom_line() 

viral.fitness %>%
  gather(key = metric, value = value, -date) %>%
  mutate(value.exp = exp(value)) %>%
  filter(metric %in% metric.beta) %>%
  select(date, metric, value.exp) %>%
  spread(key = metric, value = value.exp) -> viral.fitness.beta

viral.fitness %>%
  gather(key = metric, value = value, -date) %>%
  mutate(value.exp = exp(value)) %>%
  filter(metric %in% metric.sigma) %>%
  select(date, metric, value.exp) %>%
  spread(key = metric, value = value.exp) -> viral.fitness.sigma


viral.fitness.beta %>%
  ggplot(aes(x = date, y = meanBeta)) + geom_ribbon(aes(ymin = meanBeta-varBeta, ymax = meanBeta + varBeta), fill = "lightblue") + 
  geom_line(aes(y = meanBeta))



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
