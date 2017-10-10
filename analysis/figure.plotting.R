rm(list=ls())
source('analysis_functions.R')
source('plotting_functions.R')

library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(broom)

############################################# 
#### On one trial, figure out time points 
find_day_at_freq <- function(sim.dir, trial.meta.data, surveillance.freq, type) {
  # filter data to look at either successful or unsuccesful depending on desire

  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  if(nrow(selected.antigens) ==0) {
    return()
  } else { 
  # read in and extract information from antigen.frequencies 
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  
  col.name = paste0('surv_', surveillance.freq)
  # First Day that antigen is above surveillance threshold 
  day.freq.selected = adply(.data = selected.antigens, .margins = 1, function(antigen) {
    antigen.frequencies %>%
      filter(antigentype == antigen$postAntigen) %>%
      filter(frequency > surveillance.freq) %>%
      summarize(first.day = min(day)) -> df
  })
  names(day.freq.selected)[names(day.freq.selected) == "first.day"] = col.name
  return(day.freq.selected)
  }
}
find_day_at_freq_all <- function(dir, correct.trials, meta.data, surveillance.freq, type) {
  meta.data %>%
    dplyr::select(-cases, -simDay, -dominant.type) %>%
    filter(final.max > surveillance.freq) %>%  
    gather(key = metric, value = emergence.value, -.id, -postAntigen, -success) %>%
    mutate(postAntigen = as.numeric(postAntigen)) %>%
    arrange(postAntigen, .id) %>%
    filter(metric != "N") -> tidy.data 
  
  freq.data <- lapply(correct.trials, function(trial) {
    tidy.data %>%
      filter(.id == trial) -> trial.tidy.data

    print(paste0("made it here_", trial))    
    meta.data = find_day_at_freq(sim.dir = paste0(dir,trial), 
                                 surveillance.freq = surveillance.freq, 
                                 trial.meta.data = trial.tidy.data,
                                 type = type)
    meta.data$.id = trial
    return(meta.data)
  })
  
  empty.lists=which(sapply(freq.data,is.null))
  if(length(empty.lists > 0)) { 
    freq.data = freq.data[-empty.lists]
  }
  freq.data.all = rbindlist(freq.data, fill = TRUE)
  return(freq.data.all)
}
day_at_freq <- function(dir, correct.trials, surveillance.freq, meta.data) {

  freq.transient = find_day_at_freq_all(dir = dir, correct.trials = correct.trials,
                                        surveillance.freq =  surveillance.freq, meta.data = meta.data, type = "Transient")
  freq.yes = find_day_at_freq_all(dir = dir, correct.trials =correct.trials,
                                  surveillance.freq = surveillance.freq, meta.data = meta.data, type = "Est.")
  
  freq.both = rbind(data.frame(success = "Est.", freq.yes),
                    data.frame(success = "Transient", freq.transient))
  return(freq.both)
}
##############################################

# Designate Transient 
data.folder = "../data/tropics/tropics_20yr/"
antigen.data = create_meta_data_all(dir = data.folder)

timeseries = read_outputfiles(data.folder, "/out.timeseries.txt")

### For north, going to calculate the timeseries relative to what it's like below 10 years
infected.range = calculate_max_infected(timeseries)
antigen.data = normalize_infection(meta.data = antigen.data, infected.range)
thres = .2
antigen.frequencies = read_outputfiles(data.folder, "/out.antigenFrequencies.txt")
days.above.thres = calculate_days_above_thres(antigen.frequencies, threshold = thres)
antigen.data %>% left_join(days.above.thres, by = c("postAntigen" = "antigentype", ".id" = ".id")) -> antigen.data
already_lost = which(is.na(antigen.data$days.above))
antigen.data$days.above[already_lost] = 0

antigen.data %>%
  mutate(success = ifelse(days.above > 45, "Est.", ifelse(final.max > .1, "Transient", "no"))) -> antigen.data
         
correct.trials = c("tropics_20", "tropics_100", "tropics_14", "tropics_56", "tropics_53")
correct.trials.2 = c("tropics_20", "tropics_100", "tropics_14", "tropics_56", "tropics_53",
                     "tropics_trial", "tropics_93", "tropics_2", "tropics_70", "tropics_76", "tropics_80",
                     "tropics_84", "tropics_69")

correct.trials = unique(antigen.data$.id)
#antigen.data %>%
#  filter(.id %in% correct.trials.2) -> antigen.data.sub

two.emerge = day_at_freq(dir=data.folder, correct.trials, surveillance.freq = .02, meta.data = antigen.data)
five.emerge = day_at_freq(dir=data.folder, correct.trials, surveillance.freq = .05, meta.data = antigen.data)
ten.emerge = day_at_freq(dir = data.folder, correct.trials, surveillance.freq = .1, meta.data = antigen.data)

##########
# Have to combine what the emergence day was 

antigen.data %>%
  select(.id, day, postAntigen, success) %>%
  mutate(postAntigen = as.numeric(postAntigen)) %>%
  #filter(.id %in% correct.trials.2) %>%
  filter(success != "no") %>% 
  left_join(two.emerge, by = c(".id" = ".id", "postAntigen" = "postAntigen", "success" = "success")) %>%
  left_join(five.emerge, by = c(".id" = ".id", "postAntigen" = "postAntigen", "success" = "success"))-> data.l# %>%
  #left_join(ten.emerge, by = c(".id" = ".id", "postAntigen" = "postAntigen", "success" = "success")) -> data.l

########################
# Now sample a trial, and pick one of each type of mutation -- plot 
fig.path = "exploratory.figs/tropics_figures/"
plot_transient_success <- function(sample.trial, samples, timeseries, antigen.data, variable.set ) { 
 
   antigen.data.sub %>%
    filter(.id == sample.trial) %>%
    filter(success == "Est." | success == "Transient") %>%
    select(postAntigen, success) %>%
    mutate(cluster.number = 2:(n()+1)) -> success.types
  
  success.types[nrow(success.types)+1, ] = c(0, "Est.", 1)
  success.types %>%
    arrange(as.numeric(postAntigen)) %>%
    select(postAntigen) -> postAntigen.order
  
  success.types$postAntigen = as.numeric(success.types$postAntigen)
  success.types$cluster.number=as.numeric(success.types$cluster.number)
  
  # SAve only antigen frequencies that are part success 
  antigen.frequencies <- read.table(paste0(data.folder, sample.trial, "/out.antigenFrequencies.txt"), header = TRUE)
  antigen.frequencies %>%
    filter(antigentype %in% success.types$postAntigen) -> antigen.freq.sim
  antigen.freq.sim$antigentype=as.numeric(antigen.freq.sim$antigentype)
  #Determine the maximum number of colors going to need 
  ## setting cluster numbers 
  left_join(x = antigen.freq.sim, y = success.types, by = c("antigentype" = "postAntigen")) -> antigen
  
  
  max.color = length(unique(antigen$cluster.number))
  myColors = set_my_colors(max.color)
  
  # this is needed for plotting the lines when something happens 
  color.chart = data.frame(cbind(myColors, cluster.number = seq(1:length(myColors))))
  
  ## Need to first go in and fill all the missing values and then combine
  antigen %>%
    distinct(day, cluster.number, .keep_all = TRUE) -> part1 
  part1 %>%
    select(-cluster.number) %>%
    spread(key = antigentype, value = frequency, fill = 0) -> part2 
  part2 %>%
    gather(key = antigen.type, value = frequency, -day, -infected, -success) -> antigen.freq.long
  
  antigen.freq.long$antigen.type = as.numeric(antigen.freq.long$antigen.type)
  antigen.freq.long$success=as.factor(antigen.freq.long$success)
  success.types$postAntigen=as.numeric(success.types$postAntigen)
  success.types$success=as.factor(success.types$success)
  left_join(antigen.freq.long, success.types, by = c("antigen.type" = "postAntigen", "success" = "success")) -> antigen.freq.long
  
  antigen.freq.long %>%
    mutate(year = day/365) %>%
    mutate(prevalence = frequency*infected) %>%
    #mutate(prevalence = frequency) %>%
    filter(prevalence > 0) %>%
    ggplot(aes(x = year, y = prevalence, fill = as.factor(antigen.type))) + # Needs to be Factor 
    #  geom_line(aes(x= year, y = infected), color = "black") + 
    geom_area(color = "gray32", aes(color = as.factor(antigentype), fill = as.factor(antigen.type), alpha = success)) +
    scale_fill_manual(values = myColors) +
    labs(y = "Prevalence", x = "Years") +
    #  coord_cartesian(xlim = c(0,10)) +
    scale_alpha_discrete(range = c(1, .35)) +
    guides(fill = FALSE)  + guides(alpha = FALSE) +
    scale_x_continuous(breaks = seq(1,20,1)) -> base.plot 
  
  base.plot + geom_line(data = antigen.freq.long, aes(x = day/365, y = infected), color = "black") -> antigen.frequency
  
  ### Making additional plot that is just based on a close up view of the clusters in question
  sample.cluster = samples$cluster.number
  line.color = color.chart$myColors[sample.cluster]
  
  samples %>%
    gather(key = time.point, value = day, day, surv_0.02, surv_0.05) %>%
    mutate(day = as.numeric(day)) -> samples.l
  
  max.time = max(samples.l$day)/365
  min.time = min(samples.l$day)/365
  cluster.transient = samples$cluster.number[which(samples$success=="Transient")]
  cluster.success = samples$cluster.number[which(samples$success=="Est.")]     
  
  base.plot + geom_vline(data = samples.l, aes(xintercept = day/365, lty = time.point, color = success), size = 1.5) +
    scale_color_manual(values = as.character(line.color), 
                       labels = c(paste0("Est. - ", cluster.success),paste0("Transient. - ", cluster.transient))) + 
    coord_cartesian(xlim = c(min.time-1, max.time + 1)) +
    guides(lty = FALSE) +
    labs(col = "Antigen Fate")  +
    theme(legend.position = "bottom") -> antigen.frequency.close
  
  ########
  
  viralFitness <- read.table(paste0(data.folder, sample.trial, "/out.viralFitnessSeries.txt"), header = TRUE)
  trackAntigen <- read.table(paste0(data.folder, sample.trial, "/out.trackAntigenSeries.txt"), header = TRUE)
  timeSeries <- read.table(paste0(data.folder, sample.trial, "/out.timeSeries.txt"), header = TRUE)
  
  min.I = min(timeSeries$totalI); max.I = max(timeSeries$totalI)
  
  timeSeries %>%
    mutate(day = date*365,
           normalize.I =(totalI-min.I)/(max.I-min.I)) -> timeSeries.trial
  
  if(variable.set == "timeSeries") {
    variables = timeSeries.trial
    col.names = c("totalS","netau", "normalize.I")
  } else if(variable.set == "trackAntigen") {
    variables = trackAntigen
    col.names = c("antigenicDiversity", "antigenicTypes", "diversity", "meanLoad", "tmrca")
  } else if(variable.set == "viralFitness") {
    variables = viralFitness
    col.names = c("covBetaSigma", "meanBeta", "meanR", "meanSigma",
                  "varBeta", "varR", "varSigma")
  }
  
  samples.l %>%
    group_by(success) %>%
    summarize(min.day = min(day),
              max.day = max(day)) -> plot.day.break
  
  variables %>%
    gather(key = variable, value = value, -day) %>%
    filter(variable %in% col.names) %>%
    filter(day > (plot.day.break$min.day[which(plot.day.break$success =="Transient")]-100) & 
             day < (plot.day.break$max.day[which(plot.day.break$success =="Transient")]+100)) -> variables.transient
  
  
  variables %>%
    gather(key = variable, value = value, -day) %>%
    filter(variable %in% col.names) %>%
    filter(day > (plot.day.break$min.day[which(plot.day.break$success =="Est.")]-100) & 
             day < (plot.day.break$max.day[which(plot.day.break$success =="Est.")]+100)) -> variables.yes

    variables.sub = rbind(data.frame(success = "Transient", variables.transient),
                    data.frame(success = "Est.", variables.yes))
    variables.sub %>%
      ggplot(aes(x=day/365, y = value)) + 
      geom_line()+facet_grid(variable~success, scales = "free", ) +
      geom_smooth() +
      labs(x = "Years") +
      #coord_cartesian(xlim = c(min.time-1, max.time + 1)) +
      theme(text = element_text(size = 8),
          axis.text = element_text(size = 8)) -> variable.plot
  
  variable.plot + geom_vline(data = samples.l, aes(xintercept = day/365, lty = time.point, color = success), size = 1.5) +
    scale_color_manual(values = as.character(line.color)) + 
    guides(lty = FALSE) + guides(color = FALSE) -> variable.plot
  
  first.column <- plot_grid(antigen.frequency,antigen.frequency.close, nrow = 2)
  plot_grid(first.column, variable.plot, rel_widths = c(.9, 1.1))
  
}
variable.full = c("timeSeries", "trackAntigen", "viralFitness")

for (i in 1:10) { 
  sample.trial = sample(correct.trials, 1)
  data.l  %>%
    filter(.id == sample.trial) %>%
    mutate(cluster.number = 2:(n()+1)) %>%
    #filter(cluster.number == 13 | cluster.number == 15) -> samples
    group_by(success) %>%
    mutate(day = as.numeric(day)) %>%
    filter(day > 3*365) %>% # Ignoring the burn-in period
    sample_n(size = 1) -> samples 
  for(variable.set in variable.full)
    plot = plot_transient_success(sample.trial, samples, timeseries, antigen.data, variable.set = variable.set)
    save_plot(plot, filename = paste0(fig.path, variable.set, i, "_metrics_Test2.pdf"), base_height = 8, base_aspect_ratio = 2)
  }
}


#### Getting average time differences 
    data.l  %>%
      mutate(day = as.numeric(day)) %>%
      mutate(first.growth = surv_0.02-day,
             second.growth = surv_0.05-surv_0.02,
             third.growth = surv_0.1 - surv_0.05) %>%
      gather(key = growth.phase, value = days, first.growth, second.growth, third.growth) %>%
      ggplot(aes(x = success, y = days/365)) + facet_wrap(~growth.phase) + geom_boxplot() +
      labs(x = "Antigen Fate", y = "Length of Growth Phase (Years)") +
      scale_y_continuous(breaks = seq(0,3.5, .25))-> growth.plot1
    
    data.l %>%
      mutate(day = as.numeric(day)) %>%
      mutate(first.growth = surv_0.02-day,
             second.growth = surv_0.05-surv_0.02,
             third.growth = surv_0.1-surv_0.05) %>%
      gather(key = growth.phase, value = days, first.growth, second.growth, third.growth) %>%
      ggplot(aes(x = growth.phase, y = days/365)) + facet_wrap(~.id, scales = "free") + geom_boxplot(aes(color = success)) +
      scale_color_manual(values = c("purple", "orange")) +
      scale_x_discrete(labels = c("First", "Second", "Third")) + 
      labs(x = "Growth Phase", y = "Length of Growth Phase (Years)", color = "Antigen Fate") +
      theme(strip.text.x = element_text(size = 8),
            axis.text.x = element_text(size = 8)) -> growth.plot2
    
    grid.growth2 = plot_grid(growth.plot1, growth.plot2, rel_widths = c(.8, 1.2), )
    
    save_plot(filename = "exploratory.figs/tropics_figures/growth.plot2.pdf",
              plot = growth.plot2, base_height = 8, base_aspect_ratio = 1.8)
  
plot_antigen_frequencies <- function(sample.trial, antigen.data) { 
  
  antigen.data %>%
    filter(.id == sample.trial) %>%
    filter(success == "Est." | success == "Transient") %>%
    select(postAntigen, success) %>%
    mutate(cluster.number = 2:(n()+1)) -> success.types
  
  success.types[nrow(success.types)+1, ] = c(0, "Est.", 1)
  success.types %>%
    arrange(as.numeric(postAntigen)) %>%
    select(postAntigen) -> postAntigen.order
  
  success.types$postAntigen = as.numeric(success.types$postAntigen)
  success.types$cluster.number=as.numeric(success.types$cluster.number)
  # Save only antigen frequencies that are part success 
  antigen.frequencies <- read.table(paste0(data.folder, sample.trial, "/out.antigenFrequencies.txt"), header = TRUE)
  antigen.frequencies %>%
    filter(antigentype %in% success.types$postAntigen) -> antigen.freq.sim
  antigen.freq.sim$antigentype=as.numeric(antigen.freq.sim$antigentype)
  #Determine the maximum number of colors going to need 
  ## setting cluster numbers 
  left_join(x = antigen.freq.sim, y = success.types, by = c("antigentype" = "postAntigen")) -> antigen
  
  
  max.color = length(unique(antigen$cluster.number))
  myColors = set_my_colors(max.color)
  
  # this is needed for plotting the lines when something happens 
  color.chart = data.frame(cbind(myColors, cluster.number = seq(1:length(myColors))))
  
  ## Need to first go in and fill all the missing values and then combine
  antigen %>%
    distinct(day, cluster.number, .keep_all = TRUE) -> part1 
  part1 %>%
    select(-cluster.number) %>%
    spread(key = antigentype, value = frequency, fill = 0) -> part2 
  part2 %>%
    gather(key = antigen.type, value = frequency, -day, -infected, -success) -> antigen.freq.long
  
  antigen.freq.long$antigen.type = as.numeric(antigen.freq.long$antigen.type)
  antigen.freq.long$success=as.factor(antigen.freq.long$success)
  success.types$postAntigen=as.numeric(success.types$postAntigen)
  success.types$success=as.factor(success.types$success)
  left_join(antigen.freq.long, success.types, by = c("antigen.type" = "postAntigen", "success" = "success")) -> antigen.freq.long
  
  antigen.freq.long %>%
    mutate(year = day/365) %>%
    mutate(prevalence = frequency*infected) %>%
    #mutate(prevalence = frequency) %>%
    filter(prevalence > 0) %>%
    ggplot(aes(x = year, y = prevalence, fill = as.factor(antigen.type))) + # Needs to be Factor 
    #  geom_line(aes(x= year, y = infected), color = "black") + 
    geom_area(color = "gray32", aes(color = as.factor(antigentype), fill = as.factor(antigen.type), alpha = success)) +
    scale_fill_manual(values = myColors) +
    labs(y = "Prevalence", x = "Years") +
    #  coord_cartesian(xlim = c(0,10)) +
    scale_alpha_discrete(range = c(1, .35)) +
    guides(fill = FALSE)  + guides(alpha = FALSE) +
    scale_x_continuous(breaks = seq(1,20,1)) -> base.plot 
  
  base.plot + geom_line(data = antigen.freq.long, aes(x = day/365, y = infected), color = "black") -> antigen.frequency
  save_plot(antigen.frequency, filename = paste0(fig.path, sample.trial, "antigen.frequency.pdf"),
            base_height = 8, base_aspect_ratio = 1.8)
}
correct.trials.2 = c("tropics_100", "tropics_14", "tropics_20", "tropics_53", "tropics_56", "tropics_69",
                     "tropics_70", "tropics_76", "tropics_80", "tropics_84", "tropics_93")

for(trial in correct.trials.2) {
  sample.trial = trial
  plot_antigen_frequencies(sample.trial, antigen.data)
}


#######################################################
# Similar to the above plots, but looking at how different metrics fluctuate 

plot_normalized_trajectories <- function(sample.trial, samples, antigen.data, variable.set) { 
  
  antigen.data.sub %>%
    filter(.id == sample.trial) %>%
    filter(success == "Est." | success == "Transient") %>%
    select(postAntigen, success) %>%
    mutate(cluster.number = 2:(n()+1)) -> success.types
  
  success.types[nrow(success.types)+1, ] = c(0, "Est.", 1)
  success.types %>%
    arrange(as.numeric(postAntigen)) %>%
    select(postAntigen) -> postAntigen.order
  
  success.types$postAntigen = as.numeric(success.types$postAntigen)
  success.types$cluster.number=as.numeric(success.types$cluster.number)
  
  # SAve only antigen frequencies that are part success 
  antigen.frequencies <- read.table(paste0(data.folder, sample.trial, "/out.antigenFrequencies.txt"), header = TRUE)
  antigen.frequencies %>%
    filter(antigentype %in% success.types$postAntigen) -> antigen.freq.sim
  antigen.freq.sim$antigentype=as.numeric(antigen.freq.sim$antigentype)
  #Determine the maximum number of colors going to need 
  ## setting cluster numbers 
  left_join(x = antigen.freq.sim, y = success.types, by = c("antigentype" = "postAntigen")) -> antigen
  
  
  max.color = length(unique(antigen$cluster.number))
  myColors = set_my_colors(max.color)
  
  # this is needed for plotting the lines when something happens 
  color.chart = data.frame(cbind(myColors, cluster.number = seq(1:length(myColors))))
  ## Need to first go in and fill all the missing values and then combine
  antigen %>%
    distinct(day, cluster.number, .keep_all = TRUE) -> part1 
  part1 %>%
    select(-cluster.number) %>%
    spread(key = antigentype, value = frequency, fill = 0) -> part2 
  part2 %>%
    gather(key = antigen.type, value = frequency, -day, -infected, -success) -> antigen.freq.long
  
  antigen.freq.long$antigen.type = as.numeric(antigen.freq.long$antigen.type)
  antigen.freq.long$success=as.factor(antigen.freq.long$success)
  success.types$postAntigen=as.numeric(success.types$postAntigen)
  success.types$success=as.factor(success.types$success)
  left_join(antigen.freq.long, success.types, by = c("antigen.type" = "postAntigen", "success" = "success")) -> antigen.freq.long
  
  antigen.freq.long %>%
    mutate(year = day/365) %>%
    mutate(prevalence = frequency*infected) %>%
    #mutate(prevalence = frequency) %>%
    filter(prevalence > 0) %>%
    ggplot(aes(x = year, y = prevalence, fill = as.factor(antigen.type))) + # Needs to be Factor 
    #  geom_line(aes(x= year, y = infected), color = "black") + 
    geom_area(color = "gray32", aes(color = as.factor(antigentype), fill = as.factor(antigen.type), alpha = success)) +
    scale_fill_manual(values = myColors) +
    labs(y = "Prevalence", x = "Years") +
    #  coord_cartesian(xlim = c(0,10)) +
    scale_alpha_discrete(range = c(1, .35)) +
    guides(fill = FALSE)  + guides(alpha = FALSE) +
    scale_x_continuous(breaks = seq(1,20,1)) -> base.plot 
  
  base.plot + geom_line(data = antigen.freq.long, aes(x = day/365, y = infected), color = "black") -> antigen.frequency
  
  ### Making additional plot that is just based on a close up view of the clusters in question
  sample.cluster = samples$cluster.number
  line.color = color.chart$myColors[sample.cluster]
  
  samples %>%
    gather(key = time.point, value = day, day, surv_0.02, surv_0.05) %>%
    mutate(day = as.numeric(day)) -> samples.l
  
  max.time = max(samples.l$day)/365
  min.time = min(samples.l$day)/365
  cluster.transient = samples$cluster.number[which(samples$success=="Transient")]
  cluster.success = samples$cluster.number[which(samples$success=="Est.")]     
  
  base.plot + geom_vline(data = samples.l, aes(xintercept = day/365, lty = time.point, color = success), size = 1.5) +
    scale_color_manual(values = c("orange", "purple"), 
                       labels = c(paste0("Est. - ", cluster.success),paste0("Transient. - ", cluster.transient))) + 
  #  coord_cartesian(xlim = c(min.time-1, max.time + 1)) +
    guides(lty = FALSE) +
    labs(col = "Antigen Fate")  +
    theme(legend.position = "bottom") -> antigen.frequency.close
  
  ########
  
  viralFitness <- read.table(paste0(data.folder, sample.trial, "/out.viralFitnessSeries.txt"), header = TRUE)
  trackAntigen <- read.table(paste0(data.folder, sample.trial, "/out.trackAntigenSeries.txt"), header = TRUE)
  timeSeries <- read.table(paste0(data.folder, sample.trial, "/out.timeSeries.txt"), header = TRUE)
  
  min.I = min(timeSeries$totalI); max.I = max(timeSeries$totalI)
  
  timeSeries %>%
    mutate(day = date*365,
           normalize.I =(totalI-min.I)/(max.I-min.I)) -> timeSeries.trial
  
  if(variable.set == "timeSeries") {
    variables = timeSeries.trial
    col.names = c("netau", "normalize.I")
    time.variables = c("date", "day")
    variables$netau[!is.finite(variables$netau)] <- NA
  } else if(variable.set == "trackAntigen") {
    variables = trackAntigen
    col.names = c("antigenicDiversity", "antigenicTypes", "diversity", "meanLoad", "tmrca")
    time.variables = "day"
  } else if(variable.set == "viralFitness") {
    variables = viralFitness
    col.names = c("covBetaSigma", "meanBeta", "meanR", "meanSigma",
                  "varBeta", "varR", "varSigma")
    time.variables = c("day", "simDay")
  }
  
  samples.l %>%
    group_by(success) %>%
    summarize(min.day = min(day),
              max.day = max(day)) -> plot.day.break
  

  variables.Scaled = variables
  variables.Scaled %>%
    filter(day > 1095) -> variables.Scaled 
  time.variables = which(colnames(variables.Scaled) %in% time.variables)
  variables.Scaled[,-time.variables] <- lapply(variables.Scaled[,-time.variables],scale)
 
   variables.Scaled %>%
    gather(key = variable, value = value, -time.variables) %>%
    filter(variable %in% col.names) %>%
    ggplot(aes(x=day/365, y = value, group = variable)) + 
    geom_smooth(color = "black") + facet_wrap(~variable) + geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Years") +
    coord_cartesian(xlim = c(3,20)) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8)) -> variable.plot
  
  variable.plot + geom_vline(data = samples.l, aes(xintercept = day/365, lty = time.point, color = success), size = 1) +
    scale_color_manual(values = c("orange", "purple")) + 
    guides(lty = FALSE) -> variable.plot
  
  #first.column <- plot_grid(antigen.frequency,antigen.frequency.close, nrow = 2)
  plot_grid(antigen.frequency.close, variable.plot,nrow = 1, rel_widths =  c(.9, 1.1))
}

variable.full = c("timeSeries", "trackAntigen", "viralFitness")

for (sample.trial in correct.trials.2) { 
    data.l  %>%
    filter(.id == sample.trial) %>%
    mutate(cluster.number = 2:(n()+1)) %>%
    #filter(cluster.number == 13 | cluster.number == 15) -> samples
    group_by(success) %>%
    mutate(day = as.numeric(day)) %>%
    filter(day > 3*365) %>% # Ignoring the burn-in period
    sample_n(size = 1) -> samples 
  for(variable.set in variable.full) {
    plot = plot_normalized_trajectories(sample.trial, samples, antigen.data, variable.set = variable.set)
    save_plot(plot, filename = paste0(fig.path, variable.set, sample.trial, "_normalizedMetrics.pdf"), base_height = 8, base_aspect_ratio = 2)
  }
}



#### reading in data set 
quantitative.traits = read.csv("~/Dropbox/current_fluantigen/qualitative.traits.csv")

quantitative.traits %>%
  filter(Metric == "Normalize.I")

 <- within(data, {
  success <- factor(success, levels=c("yes", "no"), labels = c(1,0))
  .id <- factor(.id)

quantitative.traits %>%
  select(-Est.2, -Trans.2) %>%
  gather(key = success, value = behavior, Est.5, Trans.5) -> qdata.l 

qdata.l <- within(qdata.l, {
  success <- factor(success)
  behavior <- factor(behavior)
})

head(qdata.l)
qdata.l %>%
  group_by(Metric, success, behavior) %>%
  summarize(count.behavior = n()) %>%
  filter(count.behavior == max(count.behavior)) -> behavior.summary
  arrange(desc(Metric, count.behavior)) 

  
  
  ################################ Code to plot metrics, not 
  
  # Read in all files of that metric type
  
  # For each metric
  
  # For each trial find 
  
  data.l %>%
    filter(.id == sample.trial) %>%
    mutate(day = as.numeric(day)) %>%
    gather(key = surv, value = surv.day, surv_0.02, surv_0.05, surv_0.1) -> sample.data.l
 
  if(variable.set == "timeSeries") {
    
    variables = read.table(paste0(data.folder, sample.trial, "/out.timeSeries.txt"), header = TRUE)
    min.I = min(variables$totalI); max.I = max(variables$totalI)
    variables %>%
      mutate(day = date*365,
             normalize.I =(totalI-min.I)/(max.I-min.I)) -> timeSeries.trial
    col.names = c("totalS","netau", "normalize.I")
    time.variables = c("date", "day")
  } else if(variable.set == "trackAntigen") {
    variables = read.table(paste0(data.folder, sample.trial, "/out.trackAntigenSeries.txt"), header = TRUE)
    col.names = c("antigenicDiversity", "antigenicTypes", "diversity", "meanLoad", "tmrca")
    time.variables = "day"
  } else if(variable.set == "viralFitness") {
    variables = read.table(paste0(data.folder, sample.trial, "/out.viralFitnessSeries.txt"), header = TRUE)
    col.names = c("covBetaSigma", "meanBeta", "meanR", "meanSigma",
                  "varBeta", "varR", "varSigma")
    time.variables = c("day", "simDay")
  }
  
  variables %>%
    gather(key = metric, value = value, -day) %>%
    filter(metric %in% col.names) %>%
    ggplot(aes(x=day/365, y = value)) + facet_wrap(~metric, scale = "free") + geom_line() +
    labs(x = "Year") -> variable.plot
  

  sample.data.l %>%
    filter(success == "Est.") -> est.samples
  sample.data.l %>%
    filter(success == "Transient") -> trans.samples
    
  variable.plot + geom_vline(data = est.samples, aes(xintercept = surv.day/365, lty = surv),
                             color =  "orange", size = 1, alpha = .4) +
    guides(lty = FALSE) -> variable.plot.success 
  
  variable.plot + geom_vline(data = trans.samples, aes(xintercept = surv.day/365, lty = surv),
                             color =  "purple", size = 1, alpha = .4) +
    guides(lty = FALSE) -> variable.plot.trans
  
  variable.grid = plot_grid(variable.plot.success, variable.plot.trans, labels = c("Est.", "Transient"))
  variable.grid
  
  
  
  
  