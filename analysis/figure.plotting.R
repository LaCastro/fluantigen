#### 
# Designate Transient 


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
  scale_fill_manual(values = c("orange", "purple")) +
  labs(y = "Infected", x = "Year")  +
  scale_x_continuous(breaks = seq(from = 1, to = 20, by =2)) +
  guides(col = FALSE) + guides(fill = FALSE) + 
  
  
  
  
  
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

