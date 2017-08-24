rm(list=ls())
source('analysis_functions.R'); source('plotting_functions.R')

# Tropics Analysis 
exploratory.figures = "../analysis/exploratory.figs/"
tropics.folder = "../data/tropics/"

success.criteria = as.data.frame(matrix( nrow = 1, ncol = 2, data = c(180, .1)))
colnames(success.criteria) = c("length.days", "freq")

## Single Geo Analysis 

#Read and combine files 
# creates a data frame of snapshot of what the population looked like when it emergeged 
# separates mutations based on the success criteria provided 
tropics.data = create_meta_data_all(dir = tropics.folder, success.criteria)


# if one run doesn't look correct 
tropics.data = remove_trials_data(tropics.data, tropics.correct.trials)


# For each trial, calculate the max and min infection 
tropics.infected.range = calculate_max_infected(tropics.timeseries)

# Normalize the infections -- to help contextualize if it's in a peak or a trough
left_join(tropics.data, tropics.infected.range) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> tropics.data

# Calculate the age of the dominant cluster at the time the new mutation has emerged 
tropics.data = calculate_age_of_dominant(tropics.data)

########
tropics.antigen.frequencies <- read_outputfiles(tropics.folder, "/out.antigenFrequencies.txt")
tropics.antigen.frequencies <- remove_trials_data(tropics.antigen.frequencies, tropics.correct.trials)

# Plot General Differences between mutations that were successful and those that weren't 
tropics.data %>%
  select(-day, -oriAntigen, -N, -R, -cases, -simDay, -min.I, -max.I) %>%
  gather(key = metric, value = value,
         -final.max, -life.length, -success, -postAntigen, -.id) -> tropics.data.l

tropics.data.l$value = as.numeric(tropics.data.l$value)

antigen.specific.metrics <- c("distance", "mutLoad", "antigenicTypes", "dominant.freq", "age")
tropics.ant.density = plot_metric_density(data.l = tropics.data.l, metrics = antigen.specific.metrics)
save_plot(tropics.ant.density, 
          filename = "../analysis/exploratory.figs/Tantigenspecific.density.pdf", 
          base_height = 8, base_aspect_ratio = 2)

pop.dynamics <- c("ratio.I", "S", "I")
tropics.pop.dynamics.density <- plot_metric_density(tropics.data.l, pop.dynamics)

save_plot(tropics.pop.dynamics.density, 
          filename = "../analysis/exploratory.figs/T.pop.dynamics.pdf",
          base_aspect_ratio = 2)


viral.fitness.metrics <- c("diversity", "tmrca","meanR", "meanLoad", "meanBeta", "meanSigma")
viral.fitness.density  = plot_metric_density(tropics.data.l, viral.fitness.metrics)
save_plot(viral.fitness.density, 
          filename = "../analysis/exploratory.figs/Tviralfitness.density.pdf", 
          base_height = 8, base_aspect_ratio = 1.5)


#################### Just looking at  successful  Dynamics 
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
  filter(frequency > 0) %>%
  ggplot(aes(x = year, y = frequency, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id) +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  scale_x_continuous(breaks = seq(1:10))+
  guides(col = FALSE) + guides(fill = FALSE) -> freq.plot.tropics

save_plot(filename = "../analysis/exploratory.figs/freq.tropics.plot.pdf",
          freq.plot.tropics,
          base_height = 8, base_aspect_ratio = 1.5)

ant.freq.success.l %>%
  filter(.id != "tropics_11") %>%
  mutate(year = day/365) %>%
  mutate(prevalence = infected*frequency*.0025) %>%
  filter(prevalence > 0) %>%
  ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors)+
  labs(y = "Frequency", x = "Years")  +
  scale_x_continuous(breaks = seq(1:10))+
  guides(col = FALSE) + guides(fill = FALSE) -> prev.tropics.plot

save_plot(filename = "../analysis/exploratory.figs/prev.tropics.plot.pdf", 
          prev.tropics.plot,
          base_height = 8, base_aspect_ratio = 1.5)



# Calculating Life Span differences between the two 
tropics.lifespans = calculate_total_life_id(tropics.antigen.frequencies)

# Designate which ones are successful 
tropics.lifespans.df = ddply(.data = success.types, .variables = ".id", function(sim) {
  successful.types = sim$postAntigen 
  successful.types = c(0, successful.types)
  tropics.lifespans %>%
    filter(.id == sim$.id[1]) %>%
    mutate(success = ifelse(antigentype %in% successful.types, "yes", "no")) -> all.lifespans
  return(all.lifespans)
})

tropics.lifespans.df %>%
  ggplot(aes(life.length/365))+
  facet_wrap(~success, scales = "free") + 
  geom_histogram(bins = 10) +
  labs(x = "Years above threshold") -> life.span.histogram

save_plot(filename = "../analysis/exploratory.figs/T.life.span.hist.pdf",
          life.span.histogram, base_aspect_ratio = 1.5)

