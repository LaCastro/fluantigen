#### What is a successful antigen
rm(list=ls())
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(broom)


##### Given a specific criteria for success, calculate the percentage of infections it accounts for the first year, second year, etc. 
tropics.folder = "../data/tropics/"

source('analysis_functions.R')
source('plotting_functions.R')


success.criteria = as.data.frame(matrix( nrow = 1, ncol = 2, data = c(180, .1)))
colnames(success.criteria) = c("length.days", "freq")

## Single Geo Analysis 
# Read and combine files 
# creates a data frame of snapshot of what the population looked like when it emerged 
tropics.data = create_meta_data_all(dir = tropics.folder, success.criteria)

days.above.thres = calculate_days_above_thres(tropics.antigen.frequencies, threshold = .4)
tropics.data %>% left_join(days.above.thres, by = c("postAntigen" = "antigentype", ".id" = ".id")) -> tropics.data

already_lost = which(is.na(tropics.data$days.above))
tropics.data$days.above[already_lost] = 0

tropics.data %>%
  mutate(success = ifelse(days.above > 365, "yes", "no")) -> tropics.data

# Get the successful types
success.types = return_success_types(tropics.data)

### For each successtype, go into antigen frequency all, filter and full out
##### Read in the antigen frequencies and then store those in a list 
tropics.antigen.frequencies <- read_outputfiles(tropics.folder, "/out.antigenFrequencies.txt")

antigen.freq.success.df = filter_frequencies_success(success.types = success.types,
                                                     antigen.frequencies = tropics.antigen.frequencies)
# Procedure for calculating the proportion of cases over a year attributable 

year.total = calculate_year_total(antigen.freq.success.df)

year.total %>%
  group_by(.id) %>%
  summarize(distinct.clusters = length(unique(number.clusters))) -> distinct.clusters 
my.colors = brewer.pal(name = "Paired", max(distinct.clusters$distinct.clusters))

year.total %>%
  group_by(.id, year) %>%
  filter(year != 0) %>%
  mutate(prop = total/sum(total)) %>%
  ggplot(aes(x = as.factor(year), y = prop, fill = as.factor(number.clusters))) + 
  geom_bar(stat = "identity") +
  facet_wrap(~.id) + 
  scale_color_manual(values = my.colors, guide = FALSE ) +
  scale_fill_manual(values = my.colors, guide = FALSE) +
  labs(x = "Year", y = "Proportion Infected With Cluster",
       title = "Success = 180 days above 10%", fill = "Antigen Cluster") -> percent.infect.yr

save_plot(filename = "exploratory.figs/percent.infect.yr.10.pdf",plot = percent.infect.yr, base_height = 8, base_aspect_ratio = 1.5)


ant.freq.success.l = left_join(x = ant.freq.success.l, y =cluster.keypair, by=c(".id" = ".id", "antigentype" = "unique.clusters"))


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
  labs(y = "Frequency", x = "Years", fill = "Unique Cluster", title = "Success Criteria: 20%")  -> antigen.dynamics

save_plot(filename = "exploratory.figs/antigen.dyanmics.20.pdf", antigen.dynamics, base_height = 8, base_aspect_ratio = 1.5)

antigen.freq.success.df %>%
  filter(.id == "tropics_6") %>%
  group_by(antigentype) %>%
  ggplot(aes(x = day, y = frequency, group = antigentype)) + geom_line() + facet_wrap(~antigentype)





############### 2 D plots #######################
tropics.data = create_meta_data_all(dir = tropics.folder, success.criteria)

# if one run doesn't look correct 
tropics.data = remove_trials_data(tropics.data, tropics.correct.trials)

# For each trial, calculate the max and min infection 
# Normalize to help identify troughs and peaks
tropics.timeseries = read_outputfiles(tropics.folder, "/out.timeseries.txt")
tropics.infected.range = calculate_max_infected(tropics.timeseries)
tropics.data = normalize_infection(tropics.data, tropics.infected.range)

# Calculate the age of the dominant cluster at the time the new mutation has emerged 
tropics.data = calculate_age_of_dominant(tropics.data)
interaction.3 = plot_scatterplot(tropics.data, variable1 = "distance", variable2 = "antigenicTypes")
save_plot(plot = interaction.3, filename = "exploratory.figs/interaction.4.pdf",
          base_aspect_ratio = 1.8, base_height = 8)


########################
# Plot General Differences between mutations that were 
# successful and those that weren't 
# Taking out first three years 

tropics.data %>%
  filter(as.numeric(day) > 3*365) %>%
  select(-day, -oriAntigen, -N, -R, -cases, -simDay, -min.I, -max.I) %>%
  gather(key = metric, value = value,
         -final.max, -life.length, -success, -postAntigen, -.id) -> tropics.data.l

tropics.data.l$value = as.numeric(tropics.data.l$value)

antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq", "age")
tropics.ant.density = plot_metric_density(data.l = tropics.data.l, metrics = antigen.specific.metrics)
save_plot(tropics.ant.density, 
          filename = "../analysis/exploratory.figs/antigenspecific.burnout360_2.pdf", 
          base_height = 8, base_aspect_ratio = 2)

pop.dynamics <- c("ratio.I", "S", "I")
tropics.pop.dynamics.density <- plot_metric_density(tropics.data.l, pop.dynamics)

save_plot(tropics.pop.dynamics.density, 
          filename = "../analysis/exploratory.figs/pop.dynamics.burnout360_2.pdf",
          base_aspect_ratio = 2)


viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma", "varSigma", "varBeta", "covBetaSigma")

viral.fitness.density  = plot_metric_density(tropics.data.l, viral.fitness.metrics)

save_plot(viral.fitness.density, 
          filename = "../analysis/exploratory.figs/viralfitness.burnout360_2.pdf", 
          base_height = 8, base_aspect_ratio = 1.5)

plot_metric_histogram(data.l = tropics.data.l, metrics = viral.fitness.metrics)


#### Plot by trial 8/23
plot_metric_id(tropics.data.l, "antigenicTypes")
samples = sample(unique(tropics.data.l$.id), size = 5)

plot1 = plot_id(tropics.data.l, trial = samples[1], metrics = viral.fitness.metrics)
plot2 = plot_id(tropics.data.l, trial = samples[2], metrics = viral.fitness.metrics)
plot3 = plot_id(tropics.data.l, trial = samples[3], metrics = viral.fitness.metrics)
plot4 = plot_id(tropics.data.l, trial = samples[4], metrics = viral.fitness.metrics)
plot5 = plot_id(tropics.data.l, trial = samples[5], metrics = viral.fitness.metrics)

trial.viral.metrics = plot_grid(plot1, plot2,  plot3, plot4, plot5, ncol = 1)
save_plot(trial.viral.metrics, filename = "exploratory.figs/trial.viral.metrics.pdf",
          base_height = 10, base_aspect_ratio = 1.6)

#################################### 8/24/17
thresholds = c(.15, .2)
days = seq(30,90,20)


tropics.timeseries = read_outputfiles(tropics.folder, "/out.timeseries.txt")
tropics.antigen.frequencies <- read_outputfiles(tropics.folder, "/out.antigenFrequencies.txt")


# changing create_meta_data so that it doesn't assign criteria for success 

for(thres in thresholds) {
  for(day.length in days) {
    
    tropics.data = create_meta_data_all(dir = tropics.folder, success.criteria)
    days.above.thres = calculate_days_above_thres(tropics.antigen.frequencies, threshold = thres)
    tropics.data %>% left_join(days.above.thres, by = c("postAntigen" = "antigentype", ".id" = ".id")) -> tropics.data
    
    already_lost = which(is.na(tropics.data$days.above))
    tropics.data$days.above[already_lost] = 0
    
    tropics.data %>%
      mutate(success = ifelse(days.above > day.length, "yes", "no")) -> tropics.data
    
    #### Plotting
    
    # Antigen Dynamics
    success.types = return_success_types(tropics.data )
    antigen.freq.success.df = filter_frequencies_success(success.types = success.types, 
                                                         antigen.frequencies = tropics.antigen.frequencies)
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
    
    ant.freq.success.l %>%
      mutate(year = day/365) %>%
     # mutate(prevalence = frequency*infected) %>%
      mutate(prevalence = frequency) %>%
      #
      filter(prevalence > 0) %>%
      ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
      geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
      #geom_line(aes(x = year, y = infected), color = "black") + 
      facet_wrap(~.id, scales = "free_y") +
      scale_color_manual(values = myColors) + 
      scale_fill_manual(values = myColors) +
      labs(y = "Frequency", x = "Years") +
      scale_x_continuous(breaks = seq(1:10)) +
      guides(col = FALSE) + guides(fill = FALSE) -> prev.plot
    
    save_plot(filename = paste0("exploratory.figs/prev.day", day.length, "thres", thres,".pdf"), prev.plot,
              base_height = 8, base_aspect_ratio = 1.6)
  }
}