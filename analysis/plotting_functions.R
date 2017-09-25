##### Plotting Functions

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

plot_timeseries_id <- function(timeseries) {
  #plot the infected timeseries, wrapping by trial
  #eventually will change this to a sample with more trials

  timeseries %>%
    dplyr::filter(date < 20) %>% # just when it's longer
    ggplot(aes(x = date, y = totalI*.0025)) + geom_line() +
    labs(x = "Year", y = "Infected per 100,000") +
    facet_wrap(~.id, scales = "free") +
    scale_x_continuous(breaks = seq(from = 1, to = 20, by = 2)) +
    theme(axis.text.x = element_text(size = 8))
}

set_my_colors <- function(length.of.success) {
  myColors = colorRampPalette(brewer.pal(8, "Accent"))(length.of.success)
  return(myColors)
}

plot_metric_density <- function(data.l, metrics) {
  # plot the differences in variables based on success
  # wraps for different metrics 
  data.l %>%
    filter(metric %in% metrics) %>%
    ggplot(aes(value, fill = success, color = success)) + 
    geom_density(alpha = .5, adjust  = 2) + 
    facet_wrap(~metric, scales = "free") +
    scale_color_manual(values = c("purple", "orange")) +
    scale_fill_manual(values = c("purple", "orange")) +
    theme(axis.text.x = element_text(size = 8)) 
}

plot_metric_histogram <- function(data.l, metrics) {
  # Plots differences in variables according to success
  # histogram
  data.l %>%
    filter(metric %in% metrics) %>%
    ggplot(aes(value)) + 
    geom_histogram(bins = 10, aes(fill = success), alpha =.8) + 
    facet_grid(success~metric, scales = "free") +
  #  scale_color_manual(values = c("purple", "orange")) +
    scale_fill_manual(values = c("purple", "orange")) +
    theme(axis.text.x = element_text(size = 8), legend.position = "none") + 
    labs(x = "Count", y  = "Value")
}

plot_full_output <- function(timeseries, trial, variable.set) {
  # Plot the timeseries of a particluar metric set for a particular run 
  
  timeseries %>%
    filter(.id == trial) %>%
    #ggplot(aes(x = day/365, y = totalI)) + geom_line() +
    ggplot(aes(x = date, y = totalI)) + geom_line() +
    labs(x = "Year", y = "Infected", title = trial) +
    scale_x_continuous(breaks = seq(1:20))-> time.series
  
  success.criteria = as.data.frame(matrix( nrow = 1, ncol = 2, data = c(180, .1)))
  colnames(success.criteria) = c("length.days", "freq")
  
  # Plot succesful dynamics
  tropics.data = create_meta_data(sim.dir = paste0(tropics.folder, trial))
  thres = .2
  tropics.antigen.frequencies = read_outputfiles(tropics.folder, "/out.antigenFrequencies.txt")
  days.above.thres = calculate_days_above_thres(tropics.antigen.frequencies, threshold = thres)
  tropics.data %>% left_join(days.above.thres, by = c("postAntigen" = "antigentype")) -> tropics.data
  
  already_lost = which(is.na(tropics.data$days.above))
  tropics.data$days.above[already_lost] = 0
  
  tropics.data %>%
    mutate(success = ifelse(days.above > 45, "yes", "no")) -> tropics.data
  
  
  tropics.data %>%
    filter(success == "yes") %>%
    select(postAntigen) -> success.types
  
  tropics.antigen.frequencies <- read.table(paste0(tropics.folder, trial, "/out.antigenFrequencies.txt"), header = TRUE)
  
  tropics.antigen.frequencies %>%
    filter(antigentype %in% success.types$postAntigen | antigentype == 0) -> antigen.freq.sim
  #Determine the maximum number of colors going to need 
  num.transitions = unique(antigen.freq.sim$antigentype)
  myColors = colorRampPalette(brewer.pal(8, "Accent"))(length(num.transitions))
  
  ## Need to first go in and fill all the missing values and then combine
  antigen.freq.sim %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2, - infected)  -> antigen.freq.long
  antigen.freq.long$antigentype = as.factor(antigen.freq.long$antigentype)
  
  antigen.freq.long %>%
    mutate(year = day/365) %>%
    mutate(prevalence = infected*frequency) %>% #*.0025
    filter(prevalence > 0) %>%
    ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
    geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
    scale_color_manual(values = myColors) + 
    scale_fill_manual(values = myColors)+
    labs(y = "Infected", x = "Year")  +
    scale_x_continuous(breaks = seq(from = 1, to = 20, by =2)) +
    guides(col = FALSE) + guides(fill = FALSE) -> prev.plot
  
  
  # Determine the columns needed depending on the variable set 
  if(variable.set == "timeSeries") {
    variables = timeseries
    col.names = c(".id", "totalS","netau")
  } else if(variable.set == "trackAntigen") {
    variables = trackAntigen
    col.names = c(".id", "antigenicDiversity", "antigenicTypes", "diversity", "meanLoad", "tmrca")
  } else if(variable.set == "viralFitness") {
    variables = viralFitness
    col.names = c("covBetaSigma", "meanBeta", "meanR", "meanSigma",
                  "varBeta", "varR", "varSigma")
  }
  
  # Plot the variable set 
  variables %>%
    gather(key = variable, value = value, -.id, -day) %>%
    filter(variable %in% col.names) %>%
    filter(.id == trial) %>%
    ggplot(aes(x=day/365, y = value)) + 
    geom_line()+facet_wrap(~variable, nrow = ifelse(variable.set == "timeseries", 2,3), scales = "free") +
    labs(x = "Year") +
    scale_x_continuous(breaks = seq(from = 1, to = 20, by =2)) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8)) -> variable.plot
  
  first.column <- plot_grid(time.series, prev.plot, nrow = 2)
  plot_grid(first.column, variable.plot, rel_widths = c(.9, 1.1))
}


plot_scatterplot <- function(meta.data, variable1, variable2, n.sample) {
  # Plot relationship between two variables 
  # wraps for difference samples 
  meta.data$variable1 = as.numeric(as.character(meta.data$variable1))
  meta.data$variable2 = as.numeric(as.character(meta.data$variable2))
  
  no.indices = sample(which(meta.data$success == "no"), size = n.sample, replace = FALSE)
  sub.indices = c(which(meta.data$success == "yes"), no.indices)
  
  sub.meta.data = meta.data[sub.indices,]
  
  sub.meta.data %>%
    select(.id, variable1,variable2, success) %>%
    ggplot(aes(x = variable1, y = variable2, col = success)) +
    geom_point(alpha = .85) + facet_wrap(~.id) +
    scale_color_manual(values = c("purple", "orange"), labels = c("Lost", "Est.")) + 
    labs(x = eval(parse(text = variable1)), 
         y = eval(parse(text = variable2))) -> interaction
  return(interaction)
}

plot_metric_id <- function(data, m) {
  # Plots a single metric but wraps for different trials 
  data.l %>%
    filter(metric == m) %>%
    ggplot(aes(value, fill = success, color = success)) + 
    geom_density(alpha = .5, adjust  = 2) + 
    facet_wrap(~.id, scales = "free") +
    scale_color_manual(values = c("purple", "orange")) +
    scale_fill_manual(values = c("purple", "orange")) +
    theme(axis.text.x = element_text(size = 8)) +
    labs(title = m) 
}

plot_id <- function(tropics.data.l, trial, metrics) {
  # Plot all the specificed metrics for a particlular trial 
  data.l %>%
    filter(.id == trial) %>%
    filter(metric %in% metrics) %>%
    ggplot(aes(value, fill = success, color = success)) + 
    geom_density(alpha = .5, adjust  = 2) + 
    facet_wrap(~metric, scales = "free",nrow = 1) +
    scale_color_manual(values = c("purple", "orange")) +
    scale_fill_manual(values = c("purple", "orange")) +
    theme(axis.text.x = element_text(size = 8), legend.position = "none") +
    labs(title = trial) 
}

fill_antigen_values <- function(meta.data.df) {
  # fill in data frame to help with plotting antigen dynamics 
  ant.freq.success.l = ddply(.data = meta.data.df, .variables = ".id", function(sim) {
    sim %>%
      distinct(day, antigentype, .keep_all = TRUE) %>%
      spread(key = antigentype, value = frequency, fill = 0) %>%
      gather(key = antigentype, value = frequency, -1, -2, - infected)  -> antigen.freq.long
    return(antigen.freq.long)
  })
  ant.freq.success.l$antigentype = as.factor(ant.freq.success.l$antigentype)
  return(ant.freq.success.l)
}

#display <- function(correct.sample, incorrect.sample) { 

#north.timeseries %>%
#  mutate(trial.type = ifelse(.id %in% north.correct.trial, "correct", "incorrect")) %>%
#  gather(key = metric, value = value, -.id, -date, -trial.type) %>%
#  filter((metric == "northI") & (.id == correct.sample | .id == incorrect.sample)) %>%
#  ggplot(aes(x = date, y = value *.0025)) +
#  geom_line(size = 1.5) +
#  facet_grid(metric~trial.type, scales = "free") +
#  ylab(label = "Number (per 100K)") + 
#  theme(strip.text.x = element_text(size = 8)) +
#  theme(strip.text.y = element_text(size = 6)) -> pop.dynamics

#desired.antigenic.metrics = c("antigenicTypes", "meanLoad", "diversity", "antigenicDiversity")

#north.track.antigen %>%
#  mutate(trial.type = ifelse(.id %in% north.correct.trial, "correct", "incorrect")) %>%
#  gather(key = metric, value = value, -.id, -day, -trial.type) %>%
#  filter(metric %in% desired.antigenic.metrics  & (.id == correct.sample | .id == incorrect.sample)) %>%
#  ggplot(aes(x = day/365, y = value)) +
#  geom_smooth(size = 1.5) +
#  geom_line()+
#  scale_x_continuous(limits = c(0,10), breaks  = seq(from = 0, to = 10, by = 1)) + 
#  facet_grid(metric~trial.type, scales = "free") +
#  theme(strip.text.x = element_text(size = 8)) +
#  theme(strip.text.y = element_text(size = 6))  -> antigenic.dynamics

#desired.fitness.metrics = c("meanR", "meanBeta")

#north.track.fitness %>%
#  mutate(trial.type = ifelse(.id %in% north.correct.trial, "correct", "incorrect")) %>%
#  gather(key = metric, value = value, -.id, -day, -trial.type) %>%
#  filter(metric %in% desired.fitness.metrics & (.id == correct.sample | .id == incorrect.sample)) %>%
#  ggplot(aes(x = day/365, y = value)) +
#  geom_smooth(size = 1.5)+
#  geom_line(size = 1.5) +
#  facet_grid(metric~trial.type, scales = "free") +
#  theme(strip.text.x = element_text(size = 8)) +
#  theme(strip.text.y = element_text(size = 6)) -> fitness.dynamics

#plot_grid(pop.dynamics, antigenic.dynamics, fitness.dynamics, ncol = 1, rel_heights = c(.9, 1.3, 1))

