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
north.folder = "../data/north/"

exploratory.figures = "../analysis/exploratory.figs/"

## Single Geo Analysis 
#Read and combine files
success.criteria = setNames(data.frame(matrix(ncol = 2, nrow = 1)), c("freq","length.days"))
success.criteria$freq = .2
success.criteria$length.days = 90

antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")

north.data = create_meta_data_all(dir = north.folder, success.criteria)
north.data = remove_trials_data(north.data, north.correct.trial)

north.antigen.frequencies <- read_outputfiles(north.folder, "/out.antigenFrequencies.txt")
north.antigen.frequencies <- remove_trials_data(north.antigen.frequencies, north.correct.trial)

north.infected.range = calculate_max_infected(timeseries = north.timeseries)
north.data = normalize_infection(north.data, infected.range = north.infected.range)
north.data = calculate_age_of_dominant(north.data)

north.data %>%
  select(-day, -oriAntigen, -N, -R, -cases, -simDay,-min.I, -max.I) %>%
  gather(key = metric, value = value,
         -final.max, -life.length, -success, -postAntigen, -.id) -> north.data.l
north.data.l$value = as.numeric(north.data.l$value)

north.ant.density = plot_metric_density(data.l = north.data.l, metrics = antigen.specific.metrics)
save_plot(north.ant.density, 
          filename = "../analysis/exploratory.figs/N.antigenspecific.density.pdf", 
          base_height = 8, base_aspect_ratio = 2)

pop.dynamics <- c("ratio.I", "S", "I")

north.pop.dynamics.density <- plot_metric_density(north.data.l, pop.dynamics)
save_plot(north.pop.dynamics.density, 
          filename = "../analysis/exploratory.figs/N.pop.dynamics.pdf",
          base_aspect_ratio = 2)


viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma",
                           "varBeta", "varR", "varSigma", "covBetaSigma")
variance.fitness.metrics <- c("varBeta", "varR", "varSigma", "covBetaSigma")


viral.fitness.density  = plot_metric_density(north.data.l, viral.fitness.metrics)
save_plot(viral.fitness.density, 
          filename = "../analysis/exploratory.figs/N.viralfitness.density.pdf", 
          base_height = 8, base_aspect_ratio = 1.5)


######### Successful  Dynamics 
success.types = return_success_types(north.data )

### For each successtype, go into antigen frequency all, filter and full out
##### Read in the antigen frequencies and then store those in a list 
antigen.freq.success.df= filter_frequencies_success(antigen.frequencies = north.antigen.frequencies,
                                                    success.types)

#Determine the maximum number of colors going to need 
antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = n_distinct(antigentype)) -> num.transitions
max.color = sum(num.transitions$num.transitions)
myColors = set_my_colors(max.color)

## Need to first go in and fill all the missing values and then combine
ant.freq.success.l = fill_antigen_values(antigen.freq.success.df) 

## Need to first go in and fill all the missing values and then combine
north.all.lifespans = calculate_total_life_id(north.antigen.frequencies)


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

#################### Trying to calculate days between successull ---- work on this 



# Summary for Threshold Plot - Need to change this to 180 days as well

threshold.freq.levels = c(.05, .1, .15, .2, .25)
threshold.day.levels = c(90, 180, 270, 360)


tropics.antigen.summary = antigen.success.summary(antigen.frequencies = tropics.antigen.frequencies,threshold.freq.levels = threshold.freq.levels, threshold.day.levels = threshold.day.levels)
north.antigen.summary = antigen.success.summary(antigen.frequencies = north.antigen.frequencies.subset, threshold.freq.levels = threshold.freq.levels, threshold.day.levels = threshold.day.levels)

combined.antigen.summary = bind_rows(north.antigen.summary, tropics.antigen.summary, .id = "source")
combined.antigen.summary$source[which(combined.antigen.summary$source == 1)] = "north"
combined.antigen.summary$source[which(combined.antigen.summary$source == 2)] = "tropics"
colnames(combined.antigen.summary)[4:6] = c("lower", "mid", "upper")


combined.antigen.summary %>%
  gather(metric, value, -source, -freq, -day) %>%
  mutate(type = ifelse(value < 1, "percentage", "counts")) %>%
  filter(type == "counts") %>%
  spread(key = metric, value = value) %>%
  ggplot(aes(x = freq, y = median, color = source)) + geom_jitter(position=position_dodge(width=0.01)) +
  facet_wrap(~day, nrow = 1)+
  geom_errorbar(aes(ymin = min, ymax = max), position=position_dodge(width=0.01)) +
  scale_color_manual(values = c("purple", "orange")) +
  scale_x_continuous(breaks = seq(from = .05, to = .25,by = .05)) +
  labs(x = "Frequency Threshold Level", y = "Number of Success") +
  guides(color = FALSE) -> summary.plot1
  
combined.antigen.summary %>%
  gather(metric, value, -source, -freq, -day) %>%
  mutate(type = ifelse(value < 1, "percentage", "counts")) %>%
  filter(type == "percentage") %>%
  spread(key = metric, value = value) %>%
  ggplot(aes(x = freq, y = mid, color = source)) + geom_jitter(position=position_dodge(width=0.01)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.01)) +
  facet_wrap(~day, nrow = 1) +
  scale_color_manual(values = c("purple", "orange")) +
  scale_x_continuous(breaks = seq(from = .05, to = .25,by = .05)) +
  labs(x = "Frequency Threshold Level", y = "Percentage of Antigens") +
 theme(legend.position="bottom")-> summary.plot2


plot_grid(summary.plot1, summary.plot2, nrow = 2) -> summary.grid
save_plot(filename = "exploratory.figs//antigen.success.summary.pdf", 
          summary.grid, base_height = 8, base_aspect_ratio = 2)



###################################
# Testing quarter to see where it lines up - trying mid seasonality

antigen.data

tropics.timeseries = read_outputfiles(data.folder, "/out.timeseries.txt")
mid.timeseries = read_outputfiles(dir = "../data/mid/correct_mid20yr/", type = "/out.timeseries.txt")
north.timeseries = read_outputfiles(dir = "../data/north_20yrcorrect/", type = "/out.timeseries.txt")

head(tropics.timeseries)

tropics.timeseries %>%
  mutate(time.of.year = date-floor(date)) %>%
  mutate(year = floor(time.of.year)) %>%
  mutate(quarter = ifelse(time.of.year < .25, 1,
                          ifelse(time.of.year < .5, 2,
                                 ifelse(time.of.year < .75, 3,4)))) %>%
  mutate(quarter = as.factor(quarter)) %>%
  select(.id, date, year, totalI, time.of.year, quarter) -> tropics.timeseries
  
mid.timeseries %>%
  mutate(time.of.year = date-floor(date)) %>%
  mutate(year = floor(time.of.year)) %>%
  mutate(quarter = ifelse(time.of.year < .25, 1,
                          ifelse(time.of.year < .5, 2,
                                 ifelse(time.of.year < .75, 3,4)))) %>%
  mutate(quarter = as.factor(quarter)) %>%
  select(.id, date, year, totalI, time.of.year, quarter) -> mid.timeseries


north.timeseries %>%
  mutate(time.of.year = date-floor(date)) %>%
  mutate(year = floor(time.of.year)) %>%
  mutate(quarter = ifelse(time.of.year < .25, 1,
                          ifelse(time.of.year < .5, 2,
                                 ifelse(time.of.year < .75, 3,4)))) %>%
  mutate(quarter = as.factor(quarter)) %>%
  select(.id, date, year,  totalI, time.of.year, quarter) -> north.timeseries


unique(north.timeseries$.id)


quarter1.2 = data.frame(seq(.25,19.25, 1));colnames(quarter1.2) = "date"
quarter2.3 = data.frame(seq(.5, 19.5, 1)); colnames(quarter2.3) = "date"
quarter3.4 = data.frame(seq(.75,19.75,1)); colnames(quarter3.4) = "date"
quarter4.1 = data.frame(seq(0,20,1)); colnames(quarter4.1) = "date"

north.timeseries %>%
  filter(.id == "north_57") %>%
  filter(year < 10) %>%
  ggplot(aes(x=date, y=totalI*.0025)) + geom_line(size = 2) +
  coord_cartesian(xlim = c(0,10), ylim = c(0,200)) + 
  geom_vline(data = quarter1.2, aes(xintercept = date), color = "orange") +
  geom_vline(data = quarter2.3, aes(xintercept = date), color = "purple") +
  geom_vline(data = quarter3.4, aes(xintercept = date), color = "darkgrey") +
  geom_vline(data = quarter4.1, aes(xintercept = date), color = "black", lty  = 2) +
  scale_x_continuous(breaks = seq(0,10,1)) +
  labs(x = "Year", y = "Infected per 100K", title  = "Temperate Seasonal Forcing") -> temperate
  
unique(mid.timeseries$.id)
mid.timeseries %>%
  filter(.id == "mid_10") %>%
  filter(year < 10) %>%
  ggplot(aes(x=date, y=totalI*.0025)) + geom_line(size = 2) +
  coord_cartesian(xlim = c(0,10)) + 
  geom_vline(data = quarter1.2, aes(xintercept = date), color = "orange") +
  geom_vline(data = quarter2.3, aes(xintercept = date), color = "purple") +
  geom_vline(data = quarter3.4, aes(xintercept = date), color = "darkgrey") +
  geom_vline(data = quarter4.1, aes(xintercept = date), color = "black", lty  = 2) +
  labs(x = "Year", y = "Infected per 100K", title = "Mid-Level Seasonal Forcing") +
  scale_x_continuous(breaks = seq(0,10,1)) -> mid.level


tropics.timeseries %>%
  filter(.id == "tropics_20") %>%
  ggplot(aes(x=date, y=totalI*.0025)) + geom_line(size = 2) +
  coord_cartesian(xlim = c(0,10)) + 
  geom_vline(data = quarter1.2, aes(xintercept = date), color = "orange") +
  geom_vline(data = quarter2.3, aes(xintercept = date), color = "purple") +
  geom_vline(data = quarter3.4, aes(xintercept = date), color = "darkgrey") +
  geom_vline(data = quarter4.1, aes(xintercept = date), color = "black", lty  = 2) + 
  labs(x = "Year", y = "Infected per 100K", title = "No Seasonal Forcing") +
  scale_x_continuous(breaks = seq(0,10,1)) -> tropics

seasonal.forcing  = plot_grid(tropics, mid.level, temperate, ncol = 1) 
save_plot(seasonal.forcing, filename = "exploratory.figs/seasonal.forcing.pdf",
          base_height = 8,
          base_aspect_ratio = 1.5)
  

