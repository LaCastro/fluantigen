rm(list=ls())
source('analysis_functions.R'); source('plotting_functions.R')

# Plotting Timeseries Analysis 
tropics.folder = "../data/tropics/"

########### Tropics #######################################
tropics.timeseries = read_outputfiles(tropics.folder, "/out.timeseries.txt")

tropics.timeseries %>%
  group_by(.id) %>%
  mutate(ratio = totalI/totalS,
         normalized.ratio = (ratio-min(ratio))/(max(ratio)-min(ratio)),
         normalized.I = (totalI-min(totalI))/(max(totalI)-min(totalI)),
         day = date*365) -> tropics.timeseries

tropics.trials = unique(tropics.timeseries$.id)
tropics.correct.trials = tropics.trials[-which(tropics.trials == "tropics_11")]


################################## Different metrics
trackAntigen <- read_outputfiles(tropics.folder, "/out.trackAntigenSeries.txt")
viralFitness <- read_outputfiles(tropics.folder, "/out.viralFitnessSeries.txt")

variable.sets = c("timeSeries", "trackAntigen", "viralFitness")

plot_full_output <- function(timeseries, trial, variable.set) {
 
   # Plot the timeseries 
  timeseries %>%
    filter(.id == trial) %>%
    ggplot(aes(x = day/365, y = totalI)) + geom_line() +
    labs(x = "Year", y = "Infected", title = trial) +
    scale_x_continuous(breaks = seq(1:10))-> time.series
  
  success.criteria = as.data.frame(matrix( nrow = 1, ncol = 2, data = c(180, .1)))
  colnames(success.criteria) = c("length.days", "freq")
  
  # Plot succesful dynamics
  tropics.data = create_meta_data(sim.dir = paste0(tropics.folder, trial), success.criteria)
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
    scale_x_continuous(breaks = seq(1:10))+
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
    scale_x_continuous(breaks = seq(1:10)) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8)) -> variable.plot
  
  first.column <- plot_grid(time.series, prev.plot, nrow = 2)
  plot_grid(first.column, variable.plot, rel_widths = c(.9, 1.1))
}

trials = unique(tropics.timeseries$.id)
trials = trials[-which(trials == "tropics_11")]

variable.sets = c("timeSeries", "trackAntigen", "viralFitness")

a_ply(.data = trials, .margins = 1, function(trial) {
  for(set in variable.sets) {
    plot.full.output = plot_full_output(timeseries = tropics.timeseries, trial = trial, variable.set = set)
    save_plot(filename = paste0("exploratory.figs/full.output.", trial, ".", set, ".pdf"),
              plot = plot.full.output,
              base_height = 8, base_aspect_ratio = 1.6)
  }
})
  










