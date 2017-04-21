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

### Timeseries
output.folder = "~/Dropbox/Projects/mutantigen/north_runs/"



#create combined meta data for each entry of list 
trial.meta <- create.meta.all(dir = output.folder)
timeseries.all = read.outputfiles(output.folder, "/out.timeseries.txt")



timeseries.all %>% 
  mutate(correct.dynamics = ifelse(.id %in% correct, "correct", "incorrect")) %>%
  gather(key = metric, value = value, -.id, -date, -correct.dynamics) %>%
  filter(metric == "totalI") %>%
  ggplot(aes(x = date, y = value, group = .id, color = correct.dynamics)) + 
  geom_line(size = 1.3)

desired.metrics = c("diversity", "tmrca", "antigenicDiversity")
colnames(timeseries.all)

timeseries.all %>%
  gather(key = metric, value = value, -.id, -date) %>%
  filter(metric %in% desired.metrics) %>%
  ggplot(aes(x = date, y = value, group = .id, color = .id)) +
  geom_smooth() +
  facet_grid(~metric, scales = "free", ncol = 1)
  

###### Would also have to read in viral fitness here to see what's going on 


########## Antigen Frequencies START HERE 
source('analysis_functions.R')
output.folder = "~/Dropbox/Projects/mutantigen/north/"
output.folder = "~/Dropbox/Projects/mutantigen/tropics/"

#Read and combine files 
tropics.data = create.meta.data.all(dir = output.folder)
antigen.freq.all <- read.outputfiles(output.folder, "/out.antigenFrequencies.txt")


###### Count the number of days each antigen is above threshold
## Create a combined figure for this 

antigen.freq.all %>%
  group_by(.id) %>%
  filter(frequency > .10) %>%
  group_by(.id, antigentype) %>%
  summarize(n = n(),
            max.freq = max(frequency)) %>%
  filter(n > 180) -> above.threshold


######### Successful  Dynamics 
tropics.data %>%
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
  guides(col = FALSE) + guides(fill = FALSE) -> freq.plot.tropics
freq.plot.tropics
save_plot(filename = paste0(exploratory.figures, "freq.plotT.pdf"), freq.plot.tropics,
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
  guides(col = FALSE) + guides(fill = FALSE) -> prev.plot.tropics

save_plot(filename = paste0(exploratory.figures, "prev.plotT.pdf"), prev.plot.tropics,
          base_height = 8, base_aspect_ratio = 1.5)



##########################
# Histograms of the metrics
antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")
pop.dynamics <- c("S", "I")
viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")



################## Creating the quantile metrics ----  Need a function for this 
tropics.data %>%
  gather(key = metric, value = value, -postAntigen, -success, -.id) -> tropics.data.long
tropics.data.long$value = as.numeric(tropics.data.long$value)

tropics.data.long %>%
  group_by(success, metric, .id) %>%
  nest(-metric) -> tropics.nested

tropics.nested %>%
  mutate(Quantiles = map(data, ~ quantile(.$value, probs = c(.025,.5,.975), na.rm = TRUE))) %>%
  unnest(map(Quantiles, tidy)) %>%
  spread(names, value = x) -> tropics.data.quantiles

colnames(tropics.data.quantiles) = c("success", "metric", "trial", "lower", "median", "upper")


# Combine the two meta popualtions 
combined.data = bind_rows(tropics.data.quantiles, north.data.quantiles, .id = "source")
combined.data$source[which(combined.data$source == 2)] = "north"
combined.data$source[which(combined.data$source == 1)] = "tropics"


combined.data %>%
  filter(metric %in% viral.fitness.metrics) %>%
  ggplot(aes(x = success, y = median, group = trial, color = source)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("purple", "orange")) +
  facet_wrap(~metric, scales = "free") +
  labs(x = "Successful Transition", y = "Number of Individuals") +
  ggtitle("Viral Fitness Metrics 10 Trials, 40 mil people") -> viral.fitness.plot

save_plot(filename = paste0(exploratory.figures, "viral.fitness.plot.pdf"), viral.fitness.plot, base_height = 8, base_aspect_ratio = 1.5)

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

