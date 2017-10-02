trials = c("tropics_100", "tropics_20")

freq.02 %>%
  filter(.id %in% trials) -> freq.2.sub
freq.05 %>%
  filter(.id %in% trials) -> freq.5.sub
freq.10 %>%
  filter(.id %in% trials) -> freq.10.sub
freq.15 %>%
  filter(.id %in% trials ) -> freq.15.sub
freq.sub = rbind(data.frame(frequency = "2",  freq.2.sub),
                 data.frame(frequency = "5", freq.5.sub),
                 data.frame(frequency = "10", freq.10.sub), 
                 data.frame(frequency = "15", freq.15.sub))

antigen.data %>%
  filter(.id %in% trials) -> antigen.sub


###### Now get the antigen data so can plot 
success.types = return_success_types(antigen.data)
antigen.freq.success.df = filter_frequencies_success(success.types = success.types, 
                                                     antigen.frequencies = antigen.freq.sub)
success.types %>%
  filter(.id %in% trials) -> success.sub

antigen.frequencies %>%
  filter(.id %in% trials) -> antigen.freq.sub


antigen.freq.success.df %>%
  filter(day < 7300) -> antigen.freq.sucess.df

antigen.freq.success.df %>%
  group_by(.id) %>%
  summarize(num.transitions = n_distinct(antigentype)) -> num.transitions

antigen.freq.success.df = ddply(antigen.freq.success.df,.variables = ".id", function(antigen) {
  trial = antigen$.id[1]
  num.distinct.cluster = num.transitions[which(num.transitions$.id == trial), "num.transitions"]
  unique.antigens = unique(antigen$antigentype)
  cluster.pair = data.frame(cbind(unique.antigens, seq(1:num.distinct.cluster$num.transitions)))
  left_join(x = antigen, y = cluster.pair, by = c("antigentype" = "unique.antigens")) -> antigen
  colnames(antigen)[6] = "cluster.number"
  return(antigen)
})

antigen.freq.success.df %>%
  select(.id, antigentype, cluster.number) -> cluster.key

freq.sub$.id = as.factor(freq.sub$.id)
freq.sub$postAntigen = as.factor(freq.sub$postAntigen)

cluster.key 
  
indices = unique(cluster.key[,c('.id', 'antigentype')])
str(cluster.key)
cluster.key %>%
  distinct(.id, antigentype, .keep_all = TRUE) -> cluster.key

cluster.key$.id = as.factor(cluster.key$.id)
cluster.key$antigentype=as.factor(cluster.key$antigentype)

freq.sub %>%
  inner_join(cluster.key, by = c("postAntigen" = "antigentype", ".id" = ".id")) -> freq.sub



max.color = max(num.transitions$num.transitions)
myColors = set_my_colors(max.color)

## Need to first go in and fill all the missing values and then combine
ant.freq.success.l = ddply(.data = antigen.freq.success.df, .variables = ".id", function(sim) {
  sim %>%
    distinct(day, cluster.number, .keep_all = TRUE)%>%
    dplyr::select(-antigentype) -> part1
  part1 %>%
    spread(key = cluster.number, value = frequency, fill = 0) -> step2 
  step2 %>%
    gather(key = cluster.number, value = frequency, -1, -2, - infected) -> antigen.freq.long
  return(antigen.freq.long)
})

ant.freq.success.l$cluster.number = as.factor(ant.freq.success.l$cluster.number)


freq.sub %>%
  filter(cluster.number == 4) %>%
  select(frequency, postAntigen, day, .id) -> cluster.four
cluster.four$day=as.numeric(cluster.four$day)

ant.freq.success.l %>%
  mutate(year = day/365) %>%
  #mutate(prevalence = frequency*infected) %>%
  mutate(prevalence = frequency) %>%
  filter(prevalence > 0) %>%
  filter(.id %in% trials ) %>%
  ggplot(aes(x = year, y = prevalence, fill = cluster.number)) +
  geom_area(color = "black", aes(color = antigentype, fill = cluster.number)) +
  #geom_line(aes(x = year, y = infected), color = "black") + 
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") +
  coord_cartesian(xlim = c(0,10)) + 
  #scale_x_continuous(breaks = seq(1:20)) +
  guides(col = FALSE) + guides(fill = FALSE) +
  geom_vline(data = cluster.four, aes(xintercept = day/365, linetype = frequency), size = 1.5)
#-> prev.plot

ant.freq.success.l %>%
  mutate(year = day/365) %>%
  mutate(prevalence = frequency*infected) %>%
  #mutate(prevalence = frequency) %>%
  filter(prevalence > 0) %>%
  filter(.id %in% trials ) %>%
  ggplot(aes(x = year, y = prevalence, fill = cluster.number)) +
  geom_area(color = "black", aes(color = antigentype, fill = cluster.number)) +
  #geom_line(aes(x = year, y = infected), color = "black") + 
  facet_wrap(~.id, scales = "free_y") +
  scale_color_manual(values = myColors) + 
  scale_fill_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") +
  coord_cartesian(xlim = c(0,10)) + 
  #scale_x_continuous(breaks = seq(1:20)) +
  guides(col = FALSE) + guides(fill = FALSE) +
  geom_vline(data = cluster.four, aes(xintercept = day/365, linetype = frequency), size = 1.5)

prev.plot + ggplot(cluster.two, aes(day)) + geom_vline()


########

antigen.eliminated %>%
  filter(success == "transient") %>%
  rbind(antigen.success) -> antigen.relevant

correct.trials = unique(antigen.relevant$.id)

type = "transient"


find_data_at_freq <- function(sim.dir, trial.meta.data, surveillance.freq, type) {
  # filter data to look at either successful or unsuccesful depending on desire
  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  
  # read in and extract information from antigen.frequencies 
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  
  # First Day that antigen is above surveillance threshold 
  day.freq.selected = adply(.data = selected.antigens, .margins = 1, function(antigen) {
    antigen.frequencies %>%
      filter(antigentype == antigen$postAntigen) %>%
      filter(frequency > surveillance.freq) %>%
      summarize(first.day = min(day))
  }) 
  
  # If there are no antigens that reached threshold but were not successful, move on 
  if(nrow(day.freq.selected) ==0) {
    return()
  } else { 
    # read in and extract from the viral fitness
    fitness <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE)
    fit.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      fitness %>%
        filter(day == antigen$first.day) %>%
        dplyr::select(-day) -> state
      antigen.state = cbind(antigen, state)
    })
    
    # read in and extract from the timeseries 
    timeseries = read.table(paste0(sim.dir, "/out.timeseries.txt"), header = TRUE)
    
    # Gonna need the timeseries modifications in here 
    timeseries.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      timeseries %>%
        mutate(day.difference = date*365-antigen$first.day) %>% 
        filter(abs(day.difference) <= 5) %>% # helps link up the date and day
        dplyr::select(1:11)
    })
    
    left_join(
      fit.summary, timeseries.summary, by = c("postAntigen", "first.day"))  -> full.freq.summary
    full.freq.summary$first.day = as.character(full.freq.summary$first.day)
    
    ## find dominant type 
    dominant.types <- find_dominant_types_at_emerge(antigen.frequencies)
    dominant.types %>%
      filter(day %in% full.freq.summary$first.day) -> dominant.types
    colnames(dominant.types)[2] = "dominant.type"; colnames(dominant.types)[3] = "dominant.freq"
    dominant.types$day = as.character(dominant.types$day)
    
    full.freq.summary %>%
      left_join(dominant.types, by = c("first.day" = "day")) -> full.freq.summary
    
    # read in and extract track antigen closest to that day
    track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
    
    track.antigen.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      track.antigen %>%
        mutate(day.difference = abs(day-antigen$first.day)) %>%
        filter(day.difference == min(day.difference)) %>%
        dplyr::select(1, 12:14) -> track.antigen.summary
      return(track.antigen.summary)
    })
    ### find the closest to the day 
    track.antigen.summary$first.day = as.character(track.antigen.summary$first.day)
    
    full.freq.summary %>%
      left_join(track.antigen.summary, by = c("first.day", "postAntigen")) -> full.freq.summary
    
    colnames(full.freq.summary)[2] = "day"
    return(full.freq.summary)
  }
}

find_day_freq <- function(sim.dir, surveillance.freq, trial.meta.data, type) { 
  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  
  # read in and extract information from antigen.frequencies 
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  
  # First Day that antigen is above surveillance threshold 
  day.freq.selected = adply(.data = selected.antigens, .margins = 1, function(antigen) {
    antigen.frequencies %>%
      filter(antigentype == antigen$postAntigen) %>%
      filter(frequency > surveillance.freq) %>%
      summarize(first.day = min(day))
  }) 
  return(day.freq.selected)
}


#find_data_at_freq_all <- function(dir, correct.trials, meta.data, surveillance.freq, type) {
surveillance.freq = .01 

### Choosing transient and successful 
antigen.relevant %>%
    dplyr::select(-cases, -simDay, -dominant.type) %>%
    filter(final.max > surveillance.freq) %>%  
    gather(key = metric, value = emergence.value, -.id, -postAntigen, -success) %>%
    mutate(postAntigen = as.numeric(postAntigen)) %>%
    arrange(postAntigen, .id) %>%
    filter(metric == "day") -> tidy.data 


transient.freq.data <- adply(.data = correct.trials, .margins = 1,.id = NULL, function(trial) {
  tidy.data %>%
    filter(.id == trial) -> trial.tidy.data
  
    meta.data = find_day_freq(sim.dir = paste0(dir,trial), 
                            surveillance.freq = surveillance.freq, 
                            trial.meta.data = trial.tidy.data, type = "transient")
    meta.data$trial = trial
  return(meta.data)
})
successful.freq.data <- adply(.data = correct.trials, .margins = 1,.id = NULL, function(trial) {
  tidy.data %>%
    filter(.id == trial) -> trial.tidy.data
  
  meta.data = find_day_freq(sim.dir = paste0(dir,trial), 
                            surveillance.freq = surveillance.freq, 
                            trial.meta.data = trial.tidy.data, type = "yes")
  meta.data$trial = trial
  return(meta.data)
})

relevant.freq.data = rbind(data.frame(est = "success", successful.freq.data),
                           data.frame(est = "transient", transient.freq.data))

#### Will First Sample and ID , then will divide into transient and successful clusters
## Will sample one of those, and then plot 
sample.trial = sample(correct.trials, size = 1)
sample.dir = paste0(data.folder, sample.trial)
  
relevant.freq.data %>%
  filter(trial == sample.trial) %>%
  group_by(est) %>%
  sample_n(size = 1) %>%  
  left_join(tidy.data, by = c("postAntigen" = "postAntigen", "trial" = ".id")) -> sample.antigens

antigen.frequencies <- read.table(paste0(sample.dir, "/out.antigenFrequencies.txt"), header = TRUE)

antigen.frequencies %>%
  filter(antigentype %in% relevant.freq.data$postAntigen | antigentype == 0) -> antigen.freq.sim

unique(antigen.freq.sim$antigentype)
unique(relevant.freq.data$postAntigen)
relevant.freq.data %>%
  select(est, postAntigen) %>%
  right_join(antigen.freq.sim, by = c("postAntigen" = "antigentype"))



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

  