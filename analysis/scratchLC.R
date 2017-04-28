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
tropics.output.folder = "~/Dropbox/Projects/mutantigen/tropics/"
north.output.folder = "~/Dropbox/Projects/mutantigen/north/"

#create combined meta data for each entry of list 
timeseries.all = read.outputfiles(output.folder, "/out.timeseries.txt")
north.timeseries = read.outputfiles(north.output.folder, "/out.timeseries.txt")
north.viralfitness = read.outputfiles(north.output.folder, "/out.viralFitnessSeries.txt")

north.track.antigen = read.outputfiles(north.output.folder, "/out.trackAntigenSeries.txt")
north.track.fitness = read.outputfiles(north.output.folder, "/out.viralFitnessSeries.txt")
total.timeseries = rbind(efi.timeseries.0, efi.timeseries.100, efi.timeseries.300, efi.timeseries.400, north.timeseries, idcol = "source")

total.timeseries %>%
  gather(key = metric, value = value, -.id, -date, -source) %>%
  filter(metric == "totalS") %>%
  ggplot(aes(x = date, y = value, group = .id, color = .id)) +
  geom_line(size = 1.5)+
  facet_wrap(~source) +
  guides(color = FALSE)



correct.trials = c("north1", "north12", "north14", "north16", "north2", "north8")
incorrect.trials = c("north3", "north4", "north5", "north6", "north7", "north9",
                     "north10", "north11", "north13", "north15", "north17", 
                     "north18", "north19", "north20")



display <- function(correct.sample, incorrect.sample) { 
  north.timeseries %>%
    mutate(trial.type = ifelse(.id %in% correct.trials, "correct", "incorrect")) %>%
    gather(key = metric, value = value, -.id, -date, -trial.type) %>%
    filter((metric == "northS" | metric == "northI") & (.id == correct.sample | .id == incorrect.sample)) %>%
    ggplot(aes(x = date, y = value/100000)) +
    geom_line(size = 1.5) +
    facet_grid(metric~trial.type, scales = "free") +
    ylab(label = "Number (per 100K)") + 
    theme(strip.text.x = element_text(size = 8)) +
    theme(strip.text.y = element_text(size = 6)) -> pop.dynamics
    
  desired.antigenic.metrics = c("antigenicTypes", "meanLoad", "diversity", "antigenicDiversity")
  
  north.track.antigen %>%
    mutate(trial.type = ifelse(.id %in% correct.trials, "correct", "incorrect")) %>%
    gather(key = metric, value = value, -.id, -day, -trial.type) %>%
    filter(metric %in% desired.antigenic.metrics  & (.id == correct.sample | .id == incorrect.sample)) %>%
    ggplot(aes(x = day/365, y = value)) +
    geom_smooth(size = 1.5) +
    geom_line()+
    scale_x_continuous(limits = c(0,10), breaks  = seq(from = 0, to = 10, by = 1)) + 
    facet_grid(metric~trial.type, scales = "free") +
    theme(strip.text.x = element_text(size = 8)) +
    theme(strip.text.y = element_text(size = 6))  -> antigenic.dynamics
  
  desired.fitness.metrics = c("meanR", "meanBeta")

  north.track.fitness %>%
    mutate(trial.type = ifelse(.id %in% correct.trials, "correct", "incorrect")) %>%
    gather(key = metric, value = value, -.id, -day, -trial.type) %>%
    filter(metric %in% desired.fitness.metrics & (.id == correct.sample | .id == incorrect.sample)) %>%
    ggplot(aes(x = day/365, y = value)) +
    geom_smooth(size = 1.5)+
    geom_line(size = 1.5) +
    facet_grid(metric~trial.type, scales = "free") +
    theme(strip.text.x = element_text(size = 8)) +
    theme(strip.text.y = element_text(size = 6)) -> fitness.dynamics

  
  plot_grid(pop.dynamics, antigenic.dynamics, fitness.dynamics, ncol = 1, rel_heights = c(.9, 1.3, 1))
}


correct.sample = sample(correct.trials, 1)
incorrect.sample = sample(incorrect.trials, 1)
trial <- display(correct.sample, incorrect.sample)
save_plot(trial, filename = paste0(exploratory.figures, "trial5.pdf"), base_height = 8)


meanLoad = plot_grid(meanLoad.correct, meanLoad.incorrect)
meanR = plot_grid(meanR.correct, meanR.incorrect)
 = plot_grid(meanSigma.correct, meanSigma.incorrect)
antigenicTypes = plot_grid(antigenicTypes.correct, antigenicTypes.incorrect)
susceptibility = plot_grid(susceptibility.correct, susceptibility.incorrect)
popdynamics  = plot_grid(pop.dynamics.correct, pop.dynamics.incorrect)


north.track.fitness %>%
  mutate(trial.type = ifelse(.id %in% correct.trials, "correct", "incorrect")) %>%
  gather(key = metric, value = value, -.id, -day, -trial.type) -> north.fitness.l
north.track.antigen %>%
  mutate(trial.type = ifelse(.id %in% correct.trials, "correct", "incorrect")) %>%
  gather(key = metric, value = value, -.id, -day, -trial.type) -> north.antigen.l

north.factors = rbind(north.fitness.l, north.antigen.l)
head(north.factors)
  

desired.trials = c("north1", "north2", "north8", "north12")

north.timeseries %>%
  gather(key = metric, value = value, -.id, -date) %>%
  filter(metric == "northI" & .id %in% desired.trials) %>%
  ggplot(aes(x = date, y = value, group = .id, color = .id)) +
  geom_line()+
  scale_color_manual(values = set.my.colors(4))
  

my.colors = set.my.colors(10)

timeseries.all %>%
  gather(key = metric, value = value, -.id, -date) %>%
  filter(metric %in% desired.metrics) %>%
  ggplot(aes(x = date, y = value, group = .id, color = .id)) +
  geom_smooth() +
  scale_color_manual(values = my.colors)+
  facet_wrap(~metric, scales = "free")


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
  
###### Would also have to read in viral fitness here to see what's going on 


########## Antigen Frequencies START HERE 
source('analysis_functions.R')
output.folder = "~/Dropbox/Projects/mutantigen/north/"
output.folder = "~/Dropbox/Projects/mutantigen/tropics/"



#Read and combine files 
emergence.data = create.meta.data.all(dir = output.folder)
emergence.data.N = create.meta.data.all(dir = "~/Dropbox/Projects/mutantigen/north/")
antigen.frequencies <- read.outputfiles(output.folder, "/out.antigenFrequencies.txt")


emergence.data.N %>%
  select(-.id, -day, -oriAntigen, -N, -R, -cases, -simDay) %>%
  gather(key = metric, value = value, -final.max, -life.length, -success, -postAntigen) %>%
  filter(metric %in% antigen.specific.metrics) -> Nemergence.antigenspecific.l
Nemergence.antigenspecific.l$value = as.numeric(Nemergence.antigenspecific.l$value)

emergence.data.l = rbind(emergence.antigenspecific.l,  Nemergence.antigenspecific.l)
tail(emergence.data.l)
#######

antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq")

exploratory.figures = "../analysis/exploratory.figs/"

ggplot(aes(x = life.length, y = final.max, color = success))

emergence.data %>% 
  select(life.length, final.max, success) %>%
  ggplot(aes(x = life.length, y = final.max, color = success)) + geom_point() + 
  scale_color_manual(values = c("purple", "orange")) + geom_vline(xintercept = 180) + geom_hline(yintercept = .1)


emergence.data %>%
  select(-.id, -day, -oriAntigen, -N, -R, -cases, -simDay) %>%
  gather(key = metric, value = value, -final.max, -life.length, -success, -postAntigen) %>%
  filter(metric %in% antigen.specific.metrics) -> emergence.antigenspecific.l
emergence.antigenspecific.l$value = as.numeric(emergence.antigenspecific.l$value)

emergence.antigenspecific.l %>%
  ggplot(aes(x = value, y = life.length, color = success)) + geom_point(alpha = .5) + 
  facet_wrap(~metric, scales = "free_x") +
  scale_color_manual(values = c("purple", "orange")) -> scatter.antigenspecific.plot
save_plot(scatter.antigenspecific.plot, filename = paste0(exploratory.figures, "T.freq.scatter.antigenspecific.pdf"),base_height = 8,  base_aspect_ratio = 1.8)

emergence.antigenspecific.l %>%
  ggplot(aes(value, fill = success, color = success)) + 
  geom_density(alpha = .5) + 
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = c("purple", "orange")) +
  scale_fill_manual(values = c("purple", "orange")) +
  theme(axis.text.x = element_text(size = 8)) -> antigenspecific.density
save_plot(antigenspecific.density, filename = paste0(exploratory.figures, "Tantigenspecific.density.pdf"), 
          base_height = 8, base_aspect_ratio = 2)


pop.dynamics <- c("S", "I")

emergence.data %>%
  select(-.id, -day, -oriAntigen, -N, -R, -cases, -simDay) %>%
  gather(key = metric, value = value, -final.max, -life.length, -success, -postAntigen) %>%
  filter(metric %in% pop.dynamics) -> emergence.popdynamics.l
emergence.popdynamics.l$value = as.numeric(emergence.popdynamics.l$value)
emergence.popdynamics.l %>%
  ggplot(aes(x = value, y = life.length, color = success)) + geom_point(alpha = .5) + 
  facet_wrap(~metric, scales = "free_x") +
  scale_color_manual(values = c("purple", "orange")) +
  scale_x_continuous(labels = fancy_scientific) -> scatter.popdynamics.plot
save_plot(scatter.popdynamics.plot, filename = paste0(exploratory.figures, "Tlife.scatter.popdynamics.pdf"), base_aspect_ratio = 1.8)


emergence.popdynamics.l %>%
  ggplot(aes(value, fill = success, color = success)) + 
  geom_density(alpha = .5) + 
  facet_wrap(~metric, scales = "free_x") +
  scale_color_manual(values = c("purple", "orange")) +
  scale_fill_manual(values = c("purple", "orange")) +
  scale_x_continuous(labels = fancy_scientific) +
  theme(axis.text.x = element_text(size = 8)) -> pop.dynamics.density
save_plot(pop.dynamics.density, filename = paste0(exploratory.figures, "Tpop.dynamics.density.pdf"), base_aspect_ratio = 2)
                                                                     



viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")

emergence.data %>%
  select(-.id, -day, -oriAntigen, -N, -R, -cases, -simDay) %>%
  gather(key = metric, value = value, -final.max, -life.length, -success, -postAntigen) %>%
  filter(metric %in% viral.fitness.metrics) -> emergence.viral.fitness.l
emergence.viral.fitness.l$value = as.numeric(emergence.viral.fitness.l$value)

emergence.viral.fitness.l %>%
  ggplot(aes(x = value, y = life.length, color = success)) + geom_point(alpha = .5) + 
  facet_wrap(~metric, scales = "free_x") +
  scale_color_manual(values = c("purple", "orange")) -> scatter.viral.fitness.plot
save_plot(scatter.viral.fitness.plot, filename = paste0(exploratory.figures, "Tlife.scatter.viral.fitness.pdf"),base_height = 8,  base_aspect_ratio = 1.8)


emergence.viral.fitness.l %>%
  ggplot(aes(value, fill = success, color = success)) + geom_density(alpha = .5) + 
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = c("purple", "orange")) +
  scale_fill_manual(values = c("purple", "orange")) +
  theme(axis.text.x = element_text(size = 8)) -> density.viral.fitness
save_plot(density.viral.fitness, filename = paste0(exploratory.figures, "Tviralfitness.density.pdf"), base_height = 8, base_aspect_ratio = 1.5)


emergence.antigen.metrics.l$value = as.numeric(emergence.antigen.metrics.l$value)





######### Successful  Dynamics 
emergence.data %>%
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
desired.trials = c("north1", "north12", "north14", "north2", "north16")

ant.freq.success.l %>%
  filter(.id %in% desired.trials) -> ant.freq.success.l

# Calculate life spans
north.life.spans = calculate.life.spans(ant.freq.success.l)
tropics.life.spans = calculate.life.spans(ant.freq.success.l) 
lifespans = bind_rows(north.life.spans, tropics.life.spans, .id = "id")

lifespans %>%
  ggplot(aes(years, color = id, fill = id)) + geom_density(alpha = .5) + 
  scale_color_manual(values = c("purple", "orange")) +
  scale_fill_manual(values = c("purple", "orange"))
  

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
  guides(col = FALSE) + guides(fill = FALSE) -> prev.plot.north

save_plot(filename = paste0(exploratory.figures, "prev.plotT.pdf"), prev.plot.tropics,
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
