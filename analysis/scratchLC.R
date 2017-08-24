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

left_join(north.data, north.infected.range) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> north.data
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
north.data %>%
  group_by(.id) %>%
  filter(success == "yes") %>%
  select(.id, postAntigen) -> success.types

### For each successtype, go into antigen frequency all, filter and full out
##### Read in the antigen frequencies and then store those in a list 
antigen.freq.success.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  north.antigen.frequencies %>%
    filter(.id == sim$.id[1]) %>%
    filter(antigentype %in% successful.types) -> antigen.freq.sim
  return(antigen.freq.sim)
})


#Determine the maximum number of colors going to need 
antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = n_distinct(antigentype)) -> num.transitions
max.color = sum(num.transitions$num.transitions)
myColors = set_my_colors(max.color)

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
north.all.lifespans = calculate_total_life_id(north.antigen.frequencies)

# Designate which ones are successful 
north.lifespans.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  north.all.lifespans %>%
    filter(.id == sim$.id[1]) %>%
    mutate(success = ifelse(antigentype %in% successful.types, "yes", "no")) -> all.lifespans
  return(all.lifespans)
})

north.lifespans.df %>%
  ggplot(aes(life.length/365))+
  facet_wrap(~success, scales = "free") + 
  geom_histogram(bins = 10) -> life.span.histogram

save_plot(filename = "../analysis/exploratory.figs/Nlife.span.hist.pdf",
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
  guides(col = FALSE) + guides(fill = FALSE) -> freq.plot.north

save_plot(filename = "../analysis/exploratory.figs/freq.north.plot.pdf",
          freq.plot.north,
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
  guides(col = FALSE) + guides(fill = FALSE) -> prev.north.plot

save_plot(filename = "../analysis/exploratory.figs/prev.plotN.pdf", 
          prev.north.plot,
          base_height = 8, base_aspect_ratio = 1.5)


### Combined Overview
combined.lifespans = rbind(data.frame(id = "tropics", tropics.lifespans.df),
                      data.frame(id = "north", north.lifespans.df))

combined.lifespans %>%
  filter(success == "yes") %>%
  ggplot(aes(life.length/365)) +
  geom_histogram(aes(y  = ..density.., color = id, fill = id)) +
  geom_density(aes(color = id, fill = id), alpha = .5) +
  scale_fill_manual(values = c("purple", "orange")) +
  scale_color_manual(values = c("purple", "orange"))+
  labs(x = "Years above 20% Frequency")


combined.lifespans %>%
  filter(success == "yes") %>%
  ggplot(aes(life.length/365)) +
  geom_histogram(aes(y  = ..density..), binwidth = .5) +
  geom_density(color = "blue", fill = "blue", alpha = .5) +
  facet_wrap(~id) +
  scale_fill_manual(values = c("purple", "orange")) +
  scale_color_manual(values = c("purple", "orange"))+
  labs(x = "Years above 20% Frequency") -> combined.lifespans

save_plot(filename = "exploratory.figs/combined.lifespans.pdf", combined.lifespans)


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
  
  
  
######################## Combined Overview 
combined.data = rbind(data.frame(id = "tropics", tropics.data.l),
                       data.frame(id = "north", north.data.l))
  

antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")
pop.dynamics <- c("ratio.I", "S", "I")


variance.fitness.metrics <- c("varBeta", "varR", "varSigma", "covBetaSigma")


combined.data %>%
  filter(metric == "covBetaSigma") %>%
  ggplot(aes(value, fill = success, color = success)) + 
  geom_density(alpha = .5, adjust  = 3) + 
  facet_grid(~id) +
  scale_color_manual(values = c("purple", "orange")) +
  scale_fill_manual(values = c("purple", "orange")) +
  theme(axis.text.x = element_text(size = 8)) +
#  guides(col = FALSE) + guides(fill = FALSE) +
  labs(subtitle = "Co-Variance in Beta and Sigma") -> covBetaSigma


viral.variance.plot = plot_grid(varBeta, varR, varSigma, covBetaSigma, nrow = 2)
save_plot(viral.variance.plot, filename = "exploratory.figs/viral.variance.plot.pdf",
          base_height = 8,
          base_aspect_ratio = 1.8)

pop.dynamics.plot = plot_grid(ratioI, susceptibles, infected, nrow = 2)


save_plot(pop.dynamics.plot, filename = "exploratory.figs/pop.dynamics.plot.pdf",
          base_height = 8,
          base_aspect_ratio = 1.8)

antigen.specific.plot = plot_grid(cirTypes, distance, mutLoad, nrow = 2)
save_plot(antigen.specific.plot, filename = "exploratory.figs/antigenic.specific.plot.pdf",
          base_height = 8,
          base_aspect_ratio = 1.8)


dominant.type.plot  = plot_grid(dominant.freq, age.dominant)
save_plot(dominant.type.plot, filename = "exploratory.figs/dominant.type.plot.pdf", 
          base_aspect_ratio = 2)


viral.fitness.pop.plot = plot_grid(genetic.diversity,
                                   mean.beta,
                                   mean.Load,
                                   mean.R,
                                   mean.sigma,
                                   tmrca,
                                   nrow = 2)
save_plot(viral.fitness.pop.plot, filename = "exploratory.figs/viral.fitness.pop.pdf",
          base_height = 8, base_aspect_ratio = 2.5)



frequencies <- read_delim("~/Downloads/frequencies.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


rowSums(frequencies[,2:5])


colnames(frequencies)[5] = "m171K"

frequencies %>%
  select(x, m171K) -> m171K; colnames(m171K) = c("x", "frequency")



frequencies %>%
  gather(key = clade, value = frequency, -x, -m171K) %>%
  ggplot(aes(x=x, y = frequency, fill = as.factor(clade))) +
  geom_area(color = "black", aes(color = clade, fill = clade)) +
  scale_fill_manual(values = c("orange", "purple", "grey"))+
  scale_x_continuous(breaks = seq(2014.5, 2017.5, .5)) +
  geom_line(aes(x = x, y = m171K, lty = "black"), size = 2) +
  scale_linetype(labels = c("m171K")) +
  labs(x = "Year", y = "Frequency", fill = "Clade", lty = "") +
  theme(axis.text.x = element_text(size = 10)) -> empirical.frequency

save_plot(filename = "exploratory.figs/empirical.frequency.pdf", plot = empirical.frequency,
          base_aspect_ratio = 1.5)




ant.freq.success.l %>%
  mutate(year = day/365) %>%
  filter(frequency > 0) %>%
  ggplot(aes(x = year, y = frequency, fill = as.factor(number.clusters))) +
  geom_area(color = "black", aes(color = number.clusters, fill = as.factor(number.clusters))) +
  facet_wrap(~.id) +
  scale_x_continuous(breaks = seq(0:10)) +
  scale_y_continuous(breaks = seq(.1,1,.1)) + 
  scale_color_manual(values = my.colors) + 
  scale_fill_manual(values = my.colors)+
