##### Plotting Functions

set.my.colors <- function(length.of.success) {
  myColors = colorRampPalette(brewer.pal(8, "Accent"))(length.of.success)
  return(myColors)
}


plot.successful.frequency <- function(antigen.frequencies, successful.types) {
  # Plot the frequencies of successful antigens through time 
  myColors <- colorRampPalette(brewer.pal(8, "Accent"))(length(successful.types))
  antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)
  
  antigen.frequencies %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    filter(antigentype %in% successful.types) %>%
    select(-infected) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2) %>%
    mutate(year = day/365) %>%
    ggplot(aes(x = year, y = frequency)) + 
    geom_area(aes(color = antigentype, fill = antigentype)) + 
    guides(col = FALSE) + guides(fill = FALSE) +
    scale_color_manual(values = myColors) +
    scale_fill_manual(values = myColors) +
    labs(y = "Frequency", x = "Years") 
}


plot.successful.infections <- function(antigen.frequencies, successful.types) {
  # Plot the frequencies of successful antigens through time 
  myColors <- colorRampPalette(brewer.pal(8, "Accent"))(length(successful.types))
  antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)
  
  antigen.frequencies %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    filter(antigentype %in% successful.types) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2, -3) %>%
    mutate(year = day/365) %>%
    mutate(prevalence = infected*frequency) %>%
    ggplot(aes(x = year, y = prevalence)) + 
    geom_area(aes(color = antigentype, fill = antigentype)) + 
    #guides(col = FALSE) + guides(fill = FALSE) +
    scale_color_manual(values = myColors) +
    scale_fill_manual(values = myColors) +
    labs(y = "Infecteds", x = "Years") 
}

plot_scatterplot <- function(long.data) {
  ## THIS DOESN'T DO ANYTHING YET 
emergence.data.l %>% 
  select(life.length, final.max, success, source) %>%
  ggplot(aes(x = life.length, y = final.max, color = success)) + geom_point() + 
  facet_wrap(~source) +
  scale_color_manual(values = c("purple", "orange")) + geom_vline(xintercept = 180) + geom_hline(yintercept = .1)
}


plot_metric_density <- function(data.l, metrics) {
  data.l %>%
    filter(metric %in% metrics) %>%
    ggplot(aes(value, fill = success, color = success)) + 
    geom_density(alpha = .5, adjust  = 3) + 
    facet_wrap(~metric, scales = "free") +
    scale_color_manual(values = c("purple", "orange")) +
    scale_fill_manual(values = c("purple", "orange")) +
    theme(axis.text.x = element_text(size = 8)) 
}


display <- function(correct.sample, incorrect.sample) { 
  north.timeseries %>%
    mutate(trial.type = ifelse(.id %in% correct.trials, "correct", "incorrect")) %>%
    gather(key = metric, value = value, -.id, -date, -trial.type) %>%
    filter((metric == "northS" | metric == "northI") & (.id == correct.sample | .id == incorrect.sample)) %>%
    ggplot(aes(x = date, y = value *.0025)) +
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
