#### 
# Designate Transient 


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
>>>>>>> origin/master

antigen.data %>%
  mutate(success = ifelse(days.above > 45, "yes", ifelse(final.max > .15, "transient", "no"))) -> antigen.data
         
################################### 
#### On one trial, figure out time points 

find_day_at_freq <- function(sim.dir, trial.meta.data, surveillance.freq, type) {
  # filter data to look at either successful or unsuccesful depending on desire
  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  
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
                                  surveillance.freq =  surveillance.freq, meta.data = meta.data, type = "transient")
  freq.yes = find_day_at_freq_all(dir = dir, correct.trials =correct.trials,
                                   surveillance.freq = surveillance.freq, meta.data = meta.data, type = "yes")
  
  freq.both = rbind(data.frame(success = "yes", freq.yes),
                    data.frame(success = "transient", freq.transient))
  return(freq.both)
}

correct.trials = c("tropics_20", "tropics_100")

two.emerge = day_at_freq(dir=tropics.folder, correct.trials, surveillance.freq = .02, meta.data = antigen.data)
five.emerge = day_at_freq(dir=tropics.folder, correct.trials, surveillance.freq = .05, meta.data = antigen.data)

##########
# Have to combine what the emergence day was 
antigen.data %>%
  select(.id, day, postAntigen, success) %>%
  mutate(postAntigen = as.numeric(postAntigen)) %>%
  filter(.id %in% correct.trials) %>%
  filter(success != "no") %>% 
  left_join(two.emerge, by = c(".id" = ".id", "postAntigen" = "postAntigen", "success" = "success")) %>%
  left_join(five.emerge, by = c(".id" = ".id", "postAntigen" = "postAntigen", "success" = "success")) -> data.l


########################
# Now sample a trial, and pick one of each type of mutation -- plot 

sample.trial = sample(correct.trials, 1)

data.l %>%
  filter(.id == sample.trial) %>%
  group_by(success)%>%
  sample_n(size = 1) -> samples 

data.l %>%
  filter(.id == sample.trial) %>%
  left_join(cluster.pair, by = c("postAntigen" = "unique.antigens")) -> cluster.l
colnames(cluster.l)[colnames(cluster.l)=="V2"] = "cluster.number"

#########################
# Plot Timeseries 

timeseries %>%
  filter(.id == sample.trial) %>%
  #ggplot(aes(x = day/365, y = totalI)) + geom_line() +
  ggplot(aes(x = date, y = totalI)) + geom_line() +
  labs(x = "Year", y = "Infected", title = sample.trial) +
  scale_x_continuous(breaks = seq(1:20)) -> time.series

antigen.data %>%
  filter(.id == sample.trial) %>%
  filter(success == "yes" | success == "transient") %>%
  select(postAntigen, success) -> success.types

antigen.frequencies <- read.table(paste0(tropics.folder, sample.trial, "/out.antigenFrequencies.txt"), header = TRUE)
antigen.frequencies %>%
  filter(antigentype %in% success.types$postAntigen | antigentype == 0) -> antigen.freq.sim

#Determine the maximum number of colors going to need 
success.types %>%
  summarize(num.transitions = n_distinct(postAntigen)) -> num.transitions
### setting cluster numbers 
num.distinct.cluster = num.transitions$num.transitions + 1
unique.antigens = unique(antigen.freq.sim$antigentype)
cluster.pair = data.frame(cbind(unique.antigens, seq(1:num.distinct.cluster)))

left_join(x = antigen.freq.sim, y = cluster.pair, by = c("antigentype" = "unique.antigens")) -> antigen
antigen$.id = sample.trial
colnames(antigen)[5] = "cluster.number"

max.color = max(antigen$cluster.number)
myColors = set_my_colors(max.color)

## Need to first go in and fill all the missing values and then combine
  antigen %>%
    distinct(day, cluster.number, .keep_all = TRUE)%>%
    dplyr::select(-antigentype, -.id) -> part1
  part1 %>%
    spread(key = cluster.number, value = frequency, fill = 0) -> step2 
  step2 %>%
    gather(key = cluster.number, value = frequency, -1, -2, - infected) -> antigen.freq.long

antigen.freq.long$cluster.number = as.factor(antigen.freq.long$cluster.number)
antigen.freq.long$cluster.number=as.numeric(antigen.freq.long$cluster.number)

data.l %>%
  filter(.id == sample.trial) %>%
  left_join(cluster.pair, by = c("postAntigen" = "unique.antigens")) -> cluster.l
colnames(cluster.l)[colnames(cluster.l)=="V2"] = "cluster.number"

antigen.freq.long %>%
  left_join(cluster.l, by = c("cluster.number" = "cluster.number")) -> antigen.freq.long
antigen.freq.long$success[is.na(antigen.freq.long$success)] = "yes"

antigen.freq.long %>%
  mutate(year = day.x/365) %>%
  mutate(prevalence = infected*frequency) %>% #*.0025
  filter(prevalence > 0) %>%
  ggplot(aes(x = year, y = prevalence, fill = cluster.number)) +
  geom_area(color = "black", aes(color = cluster.number, fill = success)) +
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



  
  
  
  
  ###### Ploting Other Timeseries
  head(time.series)
  
  
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
