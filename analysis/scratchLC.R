rm(list=ls())
#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(broom)

source('analysis_functions.R')
source('plotting_functions.R')


########## Antigen Frequencies START HERE 
tropics.folder = "../data/tropics/"
north.folder = "../data/north/"

exploratory.figures = "../analysis/exploratory.figs/"

## Single Geo Analysis 
#Read and combine files 
antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")


tropics.data = create.meta.data.all(dir = tropics.folder)
tropics.antigen.frequencies <- read.outputfiles(tropics.folder, "/out.antigenFrequencies.txt")


tropics.data %>%
  select(-.id, -day, -oriAntigen, -N, -R, -cases, -simDay) %>%
  gather(key = metric, value = value,
         -final.max, -life.length, -success, -postAntigen) -> data.l
data.l$value = as.numeric(data.l$value)

Tropics.ant.density = plot_metric_density(data.l = data.l, metrics = antigen.specific.metrics)
save_plot(antigenspecific.density, 
          filename = "../analysis/exploratory.figs/Tantigenspecific.density.pdf", 
          base_height = 8, base_aspect_ratio = 2)

pop.dynamics <- c("S", "I")
Tropics.pop.dynamics.density <- plot_metric_density(data.l, pop.dynamics)


viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")

emergence.viral.fitness.l %>%
viral.fitness.density  = plot_metric_density(data.l, viral.fitness.metrics)
save_plot(viral.fitness.density, 
          filename = "../analysis/exploratory.figs/Tviralfitness.density.pdf", 
          base_height = 8, base_aspect_ratio = 1.5)


######### Successful  Dynamics 
tropics.data %>%
  group_by(.id) %>%
  filter(success == "yes") %>%
  select(.id, postAntigen) -> success.types

### For each successtype, go into antigen frequency all, filter and full out
##### Read in the antigen frequencies and then store those in a list 
antigen.freq.success.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  tropics.antigen.frequencies %>%
    filter(.id == sim$.id[1]) %>%
    filter(antigentype %in% successful.types) -> antigen.freq.sim
  return(antigen.freq.sim)
})



#Determine the maximum number of colors going to need 
antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = n_distinct(antigentype)) -> num.transitions
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

# Calculate life spans
tropics.life.spans = calculate.life.spans.success(ant.freq.success.l) 

tropics.life.spans %>%
  ggplot(aes(years)) +
  facet_wrap(~.id)+
  geom_histogram(binwidth = 1)
  

all.lifespans = calculate.total.life(tropics.antigen.frequencies)
# Designate which ones are successful 
all.lifespans.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  all.lifespans %>%
    filter(.id == sim$.id[1]) %>%
    mutate(success = ifelse(antigentype %in% successful.types, "yes", "no")) -> all.lifespans
  return(all.lifespans)
})

all.lifespans.df %>%
  ggplot(aes(life.length/365))+
  facet_wrap(~success, scales = "free") + 
  geom_histogram(bins = 10) -> life.span.histogram
save_plot(filename = "../analysis/exploratory.figs/life.span.hist.pdf",
          life.span.histogram, base_aspect_ratio = 1.5)


ant.freq.success.l %>%
  mutate(year = day/365) %>%
  filter(frequency > 0) %>%
  ggplot(aes(x = year, y = frequency, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id) +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  guides(col = FALSE) + guides(fill = FALSE) -> freq.plot.tropics

save_plot(filename = "../analysis/exploratory.figs/freq.tropics.plot.pdf",
          freq.plot.tropics,
          base_height = 8, base_aspect_ratio = 1.5)
  
ant.freq.success.l %>%
  mutate(year = day/365) %>%
  mutate(prevalence = infected*frequency*.0025) %>%
  filter(prevalence > 0) %>%
  ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  guides(col = FALSE) + guides(fill = FALSE) -> prev.tropics.plot

save_plot(filename = "../analysis/exploratory.figs/prev.plotT.pdf", 
          prev.tropics.plot,
          base_height = 8, base_aspect_ratio = 1.5)


### NORTH 
north.data = create.meta.data.all(dir = north.output.folder)

north.antigen.frequencies <- read.outputfiles(north.output.folder, 
                                              "/out.antigenFrequencies.txt")

north.data %>%
  filter(.id %in% correct.trial) %>%
  select(-.id, -day, -oriAntigen, -N, -R, -cases, -simDay) %>%
  gather(key = metric, value = value,
         -final.max, -life.length, -success, -postAntigen) -> data.l
data.l$value = as.numeric(data.l$value)

Tropics.ant.density = plot_metric_density(data.l = data.l, metrics = antigen.specific.metrics)
save_plot(antigenspecific.density, 
          filename = "../analysis/exploratory.figs/Tantigenspecific.density.pdf", 
          base_height = 8, base_aspect_ratio = 2)

pop.dynamics <- c("S", "I")
Tropics.pop.dynamics.density <- plot_metric_density(data.l, pop.dynamics)


viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")

emergence.viral.fitness.l %>%
  viral.fitness.density  = plot_metric_density(data.l, viral.fitness.metrics)
save_plot(viral.fitness.density, 
          filename = "../analysis/exploratory.figs/Tviralfitness.density.pdf", 
          base_height = 8, base_aspect_ratio = 1.5)


######### Successful  Dynamics 
tropics.data %>%
  group_by(.id) %>%
  filter(success == "yes") %>%
  select(.id, postAntigen) -> success.types

### For each successtype, go into antigen frequency all, filter and full out
##### Read in the antigen frequencies and then store those in a list 
antigen.freq.success.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  tropics.antigen.frequencies %>%
    filter(.id == sim$.id[1]) %>%
    filter(antigentype %in% successful.types) -> antigen.freq.sim
  return(antigen.freq.sim)
})



#Determine the maximum number of colors going to need 
antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = n_distinct(antigentype)) -> num.transitions
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

# Calculate life spans
tropics.life.spans = calculate.life.spans.success(ant.freq.success.l) 

tropics.life.spans %>%
  ggplot(aes(years)) +
  facet_wrap(~.id)+
  geom_histogram(binwidth = 1)


all.lifespans = calculate.total.life(tropics.antigen.frequencies)
# Designate which ones are successful 
all.lifespans.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  all.lifespans %>%
    filter(.id == sim$.id[1]) %>%
    mutate(success = ifelse(antigentype %in% successful.types, "yes", "no")) -> all.lifespans
  return(all.lifespans)
})

all.lifespans.df %>%
  ggplot(aes(life.length/365))+
  facet_wrap(~success, scales = "free") + 
  geom_histogram(bins = 10) -> life.span.histogram
save_plot(filename = "../analysis/exploratory.figs/life.span.hist.pdf",
          life.span.histogram, base_aspect_ratio = 1.5)


ant.freq.success.l %>%
  mutate(year = day/365) %>%
  filter(frequency > 0) %>%
  ggplot(aes(x = year, y = frequency, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id) +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  guides(col = FALSE) + guides(fill = FALSE) -> freq.plot.tropics

save_plot(filename = "../analysis/exploratory.figs/freq.tropics.plot.pdf",
          freq.plot.tropics,
          base_height = 8, base_aspect_ratio = 1.5)

ant.freq.success.l %>%
  mutate(year = day/365) %>%
  mutate(prevalence = infected*frequency*.0025) %>%
  filter(prevalence > 0) %>%
  ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  guides(col = FALSE) + guides(fill = FALSE) -> prev.tropics.plot

save_plot(filename = "../analysis/exploratory.figs/prev.plotT.pdf", 
          prev.tropics.plot,
          base_height = 8, base_aspect_ratio = 1.5)










##########################
# Histograms of the metrics
antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")
pop.dynamics <- c("S", "I")
viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")



################## Creating the quantile metrics ----  Need a function for this 
emergence.data %>%
  gather(key = metric, value = value, -postAntigen, -success, -.id) -> data.long
data.long$value = as.numeric(data.long$value)

data.long %>%
  group_by(success, metric, .id) %>%
  nest(-metric) -> north.nested

north.nested %>%
  mutate(Quantiles = map(data, ~ quantile(.$value, probs = c(.025,.5,.975), na.rm = TRUE))) %>%
  unnest(map(Quantiles, tidy)) %>%
  spread(names, value = x) -> north.data.quantiles

colnames(north.data.quantiles) = c("success", "metric", "trial", "lower", "median", "upper")


# Combine the two meta popualtions 
combined.data = bind_rows(tropics.data.quantiles, north.data.quantiles, .id = "source")
combined.data$source[which(combined.data$source == 2)] = "north"
combined.data$source[which(combined.data$source == 1)] = "tropics"


combined.data %>%
  filter(metric %in% antigen.specific.metrics) %>%
  ggplot(aes(x = success, y = median, group = trial, color = source)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("purple", "orange")) +
  facet_wrap(~metric, scales = "free") +
  labs(x = "Successful Transition", y = "Number of Individuals") +
  ggtitle("Viral Fitness Metrics 10 Trials, 40 mil people") -> antigen.specific.plot

save_plot(filename = paste0(exploratory.figures, "antigen.specific.plotCulled.pdf"), antigen.specific.plot, base_height = 8, base_aspect_ratio = 1.5)

##############################################################
#################### Trying to calculate days between successull ---- work on this 
tropics.data$day = as.numeric(tropics.data$day)

tropics.data %>%
  filter(success == "yes") %>%
  group_by(.id) %>%
  mutate(time.difference = day - lag(day)) -> tropics.data.day

north.data$day = as.numeric(north.data$day)

north.data %>%
  filter(success == "yes") %>%
  group_by(.id) %>%
  mutate(time.difference = day - lag(day)) -> north.data.day
na.index = which(is.na(north.data.day$time.difference))
north.data.day$time.difference[na.index] = north.data.day$day[na.index]


north.data.day %>%
  select(.id, time.difference) -> north.time.difference

north.time.difference$region = "north"
tropics.time.difference$region = "tropics"
tropics.data.day %>%
  select(.id, time.difference) -> tropics.time.difference

time.difference <-rbind(north.time.difference, tropics.time.difference)

time.difference %>%
  
  head(time.difference)
ggplot(time.difference, aes(time.difference, color = region, fill = region)) + geom_histogram(bins = 20) +
  scale_color_manual(values = c("purple", "orange")) +
  scale_fill_manual(values = c("purple",  "orange")) -> time.difference.plot



# Summary for Threshold Plot 
#tropics.antigen.summary = antigen.success.summary(threshold.levels = threshold.levels, antigen.frequencies = tropics.antigen)
#north.antigen.summary = antigen.success.summary(threshold.levels = threshold.levels, antigen.frequencies = north.antigen)

#combined.antigen.summary = bind_rows(north.antigen.summary, tropics.antigen.summary, .id = "source")
#combined.antigen.summary$source[which(combined.antigen.summary$source == 1)] = "north"
#combined.antigen.summary$source[which(combined.antigen.summary$source == 2)] = "tropics"
#colnames(combined.antigen.summary)[3:5] = c("lower", "median", "upper")

#geom_point(position = position_dodge(width = 0.5)) +
#  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = position_dodge(width = 0.5)) +

#combined.antigen.summary %>%
#  gather(metric, value, -source, -threshold) %>%
#  mutate(type = ifelse(value < 1, "percentage", "counts")) %>%
#  filter(type == "counts") %>%
#  spread(key = metric, value = value) %>%
#  ggplot(aes(x = threshold, y = mean, color = source)) + geom_jitter(position=position_dodge(width=0.01)) +
#  geom_errorbar(aes(ymin = min, ymax = max), position=position_dodge(width=0.01)) +
#  scale_color_manual(values = c("purple", "orange")) +
#  scale_x_continuous(breaks = seq(from = .05, to = .25,by = .05)) +
#  labs(x = "Threshold Level of Success", y = "Number of Success") +
#  guides(color = FALSE) -> summary.plot1
  
#combined.antigen.summary %>%
#  gather(metric, value, -source, -threshold) %>%
#  mutate(type = ifelse(value < 1, "percentage", "counts")) %>%
#  filter(type == "percentage") %>%
#  spread(key = metric, value = value) %>%
#  ggplot(aes(x = threshold, y = median, color = source)) + geom_jitter(position=position_dodge(width=0.01)) +
#  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.01)) +
#  scale_color_manual(values = c("purple", "orange")) +
#  scale_x_continuous(breaks = seq(from = .05, to = .25,by = .05)) +
#  labs(x = "Threshold Levels of Success", y = "Percentage of Antigens") -> summary.plot2


#plot_grid(summary.plot1, summary.plot2, nrow = 1) -> summary.grid
#save_plot(filename = paste0(exploratory.figures, "antigen.success.summary.pdf"), 
#          summary.grid, base_height = 4, base_aspect_ratio = 1.8)

######## histograms of the frequencies

#antigen.freq.all %>%
#  group_by(antigentype, .id) %>%
#  summarize(max.freq = max(frequency)) %>%
#  ggplot(aes(max.freq))+ geom_histogram(bins = 50) + 
#  coord_cartesian(ylim = c(0,50))+
#  facet_wrap(~.id) +
#  geom_vline(xintercept = .1, color = "purple") -> max.freq.hist
#save_plot("../analysis/exploratory.figs/max.freq.histT.pdf", plot = max.freq.hist,base_height = 8, base_aspect_ratio = 2)


#antigen.freq.all %>%
#  group_by(.id) %>%
#  summarize(unique.antigens = n_distinct(antigentype)) -> n.antigens


#antigen.freq.all %>%
#  group_by(.id, antigentype) %>%
#  summarize(max.freq = max(frequency)) %>%
#  ungroup(antigentype) %>%
#  group_by(.id) %>%
#  summarize(above.threshold  = sum(max.freq > .25)) -> above.threshold


desired.metrics = c("diversity", "tmrca", "netau", "serialInterval", "antigenicDiversity")


north.timeseries %>%
  mutate(north.type = ifelse(.id %in% desired.trials, "correct.run", "incorrect.run")) %>%
  mutate(i.s = totalI/totalS) %>%
  gather(key = metric, value = value, -date, -.id, -north.type) %>%
  group_by(date, metric, north.type) %>%
  summarize(average = mean(value),
            sd = sd(value)) %>%
  filter(metric %in% desired.metrics | metric == "i.s") %>%
  ggplot(aes(x = date, y = average, color = north.type, group = north.type)) + geom_line(size =  1.3) + 
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = set.my.colors(2))



north.viralfitness %>%
  mutate(north.type = ifelse(.id %in% desired.trials, "correct.run", "incorrect.run")) %>%
  gather(key = metric, value = value, -day, -.id, -north.type, -simDay) %>%
  group_by(day, metric, north.type) %>%
  summarize(average = mean(value),
            sd = sd(value)) %>%
  ggplot(aes(x = day, y = average, color = north.type, group = north.type)) + geom_line(size =  1.3) + 
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = set.my.colors(2))

  
######
north.trackAntigenic = read.outputfiles(north.output.folder, type = "/out.trackAntigenSeries.txt")

desired.metrics = c("northS", "northI", "antigenicDiversity", "diversity")

north.timeseries %>%
  filter(.id == "north3" | .id == "north4" | .id == "north1")%>%
  gather(key = metric, value = value, -.id, -date) %>%
  filter(metric == desired.metrics) %>%
  ggplot(aes(x = date, y = value, color = .id)) + geom_line() + facet_wrap(~metric, scales = "free")

north.trackAntigenic %>%
  mutate(north.type = ifelse(.id %in% desired.trials, "correct.run", "incorrect.run")) %>%
  gather(key = metric, value = value, -day, -.id, -north.type) %>%
  group_by(day, metric, north.type) %>%
  summarize(average = mean(value),
            sd = sd(value)) %>%
  head()
  ggplot(aes(x = day, y = average, color = north.type, group = north.type)) + geom_line(size =  1.3) + 
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = set.my.colors(2))
  
  
  north.lifespans <- calculate.life.spans(ant.freq.success.l)
  north.lifespans %>%
    ggplot(aes(years)) + geom_histogram(bins = 20)
