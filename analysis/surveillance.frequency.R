north.data %>%
  select(-min.I, -cases, -simDay, -max.I, -dominant.type) %>%
  gather(key = metric, value = emergence.value, -.id, -postAntigen, -success) %>%
  mutate(postAntigen = as.numeric(postAntigen)) %>%
  arrange(postAntigen, .id) %>%
  filter(metric != "N") -> north.data.tidy
unique(north.data.tidy$metric)




north.data %>%
  group_by(.id) %>%
  summarize(max.I = max(max.I),
            min.I = min(min.I)) -> summary.I


tropics.data %>%
  group_by(.id) %>%
  summarize(max.I = max(max.I),
            min.I = min(min.I)) -> summary.tropics

tropics.data %>%
  select(-min.I, -cases, -simDay, -max.I, -dominant.type) %>%
  gather(key = metric, value = emergence.value, -.id, -postAntigen, -success) %>%
  mutate(postAntigen = as.numeric(postAntigen)) %>%
  arrange(postAntigen, .id) %>%
  filter(metric != "N") -> tropics.data.tidy
unique(north.data.tidy$metric)


## have to make data tidy before this 
find_data_at_freq <- function(sim.dir, trial.meta.data, surveillance.freq, type) {
  
  # filter data to look at either successful or unsuccesful depending on desire
  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  
   # read in and extract information from antigen.frequencies 
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  day.freq.selected = adply(.data = selected.antigens, .margins = 1, function(antigen) {
    antigen.frequencies %>%
      filter(antigentype == antigen$postAntigen) %>%
      filter(frequency > surveillance.freq) %>%
      summarize(first.day = min(day))
        }) 
 
  # select just those post antigens that ever reach the surveillance frequency 
  if(type == "no") {
    day.freq.selected %>%
      filter(!is.na(first.day)) -> day.freq.selected
  }
  
  # read in and extract from the viral fitness
  fitness <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE)
  fit.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
    fitness %>%
      filter(day == antigen$first.day) %>%
      select(-day) -> state
    antigen.state = cbind(antigen, state)
  })
  
  # read in and extract from the timeseries 
  timeseries = read.table(paste0(sim.dir, "/out.timeseries.txt"), header = TRUE)
  #day = seq(from = 10, to = 10*365, by = 10) 
  #timeseries$day= day
  timeseries.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
    timeseries %>%
      mutate(day.difference = date*365-antigen$first.day) %>%
      filter(abs(day.difference) <= 5) %>%
      select(1:11)
  })
  
  left_join(
    fit.summary, timeseries.summary, by = "postAntigen")  -> full.freq.summary
  
  ## find dominant type 
  dominant.types <- find_dominant_types_at_emerge(antigen.frequencies)
  days.at.freq.thres <- full.freq.summary$first.day.x
  dominant.types %>%
    filter(day %in% days.at.freq.thres) -> dominant.types
  colnames(dominant.types)[2] = "dominant.type"; colnames(dominant.types)[3] = "dominant.freq"
  dominant.types$day = as.character(dominant.types$day)
  full.freq.summary$first.day.x = as.character(full.freq.summary$first.day.x)
  
  full.freq.summary %>%
    left_join(dominant.types, by = c("first.day.x" = "day")) -> full.freq.summary
  
  # read in and extract track antigen 
  track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
  track.antigen.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
    track.antigen %>%
      mutate(day.difference = abs(day-antigen$first.day)) %>%
      filter(day.difference == min(day.difference)) %>%
      select(12:14)
  })
   ### find the closest to the day 
  track.antigen.summary$first.day = as.character(track.antigen.summary$first.day)
  full.freq.summary %>%
    left_join(track.antigen.summary, by = c("first.day.x" = "first.day", "postAntigen")) -> full.freq.summary
  
  colnames(full.freq.summary)[2] = "day"
  return(full.freq.summary)
}


find_data_at_freq_all <- function(dir, correct.trials, tidy.data, surveillance.freq, type) {
  freq.data <- lapply(correct.trials, function(trial) {
    tidy.data %>%
      filter(.id == trial) -> trial.tidy.data
    meta.data = find_data_at_freq(sim.dir = paste0(dir,trial), 
                                  surveillance.freq = surveillance.freq, 
                                  trial.meta.data = trial.tidy.data,
                                  type = type)
    return(meta.data)
  })
  names(freq.data) = correct.trials
  freq.data.all = rbindlist(freq.data, idcol = TRUE)
  return(freq.data.all)
}


data.freq.ten.no = find_data_at_freq_all(dir = tropics.folder, correct.trials = tropics.correct.trials,
                                       surveillance.freq = .1, tidy.data = tropics.data.tidy, type = "no")
data.freq.ten.yes = find_data_at_freq_all(dir = tropics.folder, correct.trials = tropics.correct.trials,
                                          surveillance.freq = .1, tidy.data = tropics.data.tidy, type = "yes")


# This is hard coded, would need to change
colnames(data.freq.five.no)[19:23] = c("N", "S", "I", "R", "cases")
combined.freq.five = rbind(data.frame(success = "no", data.freq.five.no),
                           data.frame(success = "yes", data.freq.five.yes))

left_join(combined.freq.five, summary.tropics) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> combined.freq.five

combined.freq.five %>%
  select(-first.day.y) %>%
  gather(key = metric, value = surveillance.05, -.id, - postAntigen, -success) -> combined.freq.five.tidy

metrics = unique(combined.freq.five.tidy$metric)
metrics = metrics[-which(metrics == "day")]
metrics = metrics[-which(metrics == "simDay")]
metrics = metrics[-which(metrics == "R")]
metrics = metrics[-which(metrics == "dominant.type")]
metrics = metrics[-which(metrics == "date")]


combined.freq.five.tidy$surveillance.05 = as.numeric(combined.freq.five.tidy$surveillance.05)
combined.freq.five.tidy %>%
  #filter(.id == north.correct.trial[6]) %>%
  filter(metric %in% metrics[21:24]) %>%
  #gather(key = time.point, value = value, -.id, -postAntigen, -success,-metric) %>%
  # group_by(.id) %>%
  ggplot(aes(surveillance.05, color = success, fill = success)) + geom_density(alpha = .5) + 
  facet_wrap(~metric,scales = "free") -> fifth.combined.at.five

save_plot(fifth.combined.at.five, 
          filename = "../analysis/exploratory.figs/fifth.combined.at.five.tropics.pdf",
          base_height = 8,
          base_aspect_ratio = 1.8)



## 10% frequency
data.freq.five.tropics= find_data_at_freq_all(dir = tropics.folder, correct.trials = tropics.correct.trials,
                                      surveillance.freq = .05, tidy.data =tropics.data.tidy, type = "yes")


left_join(data.freq.five.tropics, summary.tropics) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> data.freq.five.tropics

data.freq.five.tropics %>%
  select(-first.day.y) %>%
  gather(key = metric, value = surveillance.05, -.id, - postAntigen) ->  data.freq.five.tropics.tidy
left_join(tropics.data.tidy, data.freq.five.tropics.tidy) -> sample.merge



## 10% frequency
data.freq.ten.tropics = find_data_at_freq_all(dir = tropics.folder, correct.trials = tropics.correct.trials,
                                              surveillance.freq = .1, tidy.data =tropics.data.tidy, type = "yes")
colnames(data.freq.ten.tropics)[19:23] = c("N", "S", "I", "R", "cases")
left_join(data.freq.ten.tropics,summary.tropics) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> data.freq.ten.tropics


data.freq.ten.tropics %>%
  select(-first.day.y) %>%
  gather(key = metric, value = surveillance.10, -.id, - postAntigen) -> data.freq.ten.tropics.tidy
left_join(sample.merge, data.freq.ten.tropics.tidy) -> sample.merge


data.freq.fifteen.tropics = find_data_at_freq_all(dir = tropics.folder, correct.trials = tropics.correct.trials,
                                                  surveillance.freq = .15, tidy.data =tropics.data.tidy, type = "yes")

colnames(data.freq.fifteen.tropics)[19:23] = c("N", "S", "I", "R", "cases")

left_join(data.freq.fifteen.tropics, summary.tropics) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> data.freq.fifteen.tropics

data.freq.fifteen.tropics %>%
  select(-first.day.y) %>%
  gather(key = metric, value = surveillance.15, -.id, - postAntigen) -> data.freq.fifteen.tropics.tidy
left_join(sample.merge, data.freq.fifteen.tropics.tidy) -> sample.merge

data.freq.twenty.tropics = find_data_at_freq_all(dir = tropics.folder, correct.trials = tropics.correct.trials,
                                                 surveillance.freq = .2, tidy.data =tropics.data.tidy, type = "yes")
colnames(data.freq.twenty.tropics)[19:23] = c("N", "S", "I", "R", "cases")

left_join(data.freq.twenty.tropics, summary.tropics) %>%
  mutate(ratio.I = (infected-min.I)/(max.I-min.I)) -> data.freq.twenty.tropics

data.freq.twenty.tropics %>%
  select(-first.day.y) %>%
  gather(key = metric, value = surveillance.20, -.id, - postAntigen) -> data.freq.twenty.tropics.tidy
left_join(sample.merge,  data.freq.twenty.tropics.tidy) -> sample.merge


excluded.metrics = c("oriAntigen", "distance", "mutLoad", "final.max", "life.length")

sample.merge %>%
  filter(!(metric %in% excluded.metrics)) %>%
  filter(success == "yes") %>%
  distinct(.id, postAntigen, success, metric, .keep_all = TRUE) -> surveillance.subset


surveillance.subset$surveillance.05=as.numeric(surveillance.subset$surveillance.05); surveillance.subset$surveillance.10=as.numeric(surveillance.subset$surveillance.10)
surveillance.subset$surveillance.15=as.numeric(surveillance.subset$surveillance.15); surveillance.subset$surveillance.20=as.numeric(surveillance.subset$surveillance.20)
surveillance.subset$emergence.value = as.numeric(surveillance.subset$emergence.value)


surveillance.subset %>%
  filter(metric == "day") %>%
  mutate(time.to.five = surveillance.05 - emergence.value,
         time.to.ten = surveillance.10 - surveillance.05,
         time.to.fifteen = surveillance.15-surveillance.10,
         time.to.twenty = surveillance.20-surveillance.15) %>%
  select(.id, postAntigen, success, time.to.five, time.to.ten, time.to.fifteen, time.to.twenty) -> surveillance.time.differences
  
colnames(surveillance.time.differences)[4:7] = c("surveillance.05", "surveillance.10", "surveillance.15", "surveillance.20")

surveillance.time.differences$metric = "time.difference"
surveillance.time.differences$emergence.value = NA
rbind(surveillance.subset, surveillance.time.differences) -> surveillance.subset


metrics = unique(surveillance.subset$metric)
metrics = metrics[-which(metrics == "day")]
metrics = metrics[-which(metrics == "R")]


surveillance.subset %>%
  #filter(.id == north.correct.trial[6]) %>%
  filter(metric %in% metrics[21:22]) %>%
  gather(key = time.point, value = value, -.id, -postAntigen, -success,-metric) %>%
 # group_by(.id) %>%
  ggplot(aes(value, color = time.point, fill = time.point)) + geom_density(alpha = .5) + 
  facet_wrap(~metric,scales = "free") -> full.surveillance.five.tropics

save_plot(full.surveillance.five.tropics,
          filename = "../analysis/exploratory.figs/full.surveillance.five.tropics.pdf",
          base_height = 8,
          base_aspect_ratio = 1.8)


surveillance.subset %>%
  filter(metric == "time.difference") %>%
  gather(key = time.point, value = value, -.id, -postAntigen, -success,-metric) %>%
  ggplot(aes(value/365)) + geom_histogram(binwidth = 45/365) + 
  facet_grid(time.point~metric) 


surveillance.subset %>%
  filter(metric == "time.difference") %>%
  select(-emergence.value) %>%
  gather(key = time.point, value = value, -.id, -postAntigen, -success,-metric) %>%
  ggplot(aes(value)) + geom_histogram(binwidth = 45) + 
  facet_grid(time.point~metric) +
  scale_x_continuous(breaks = seq(0,600,45)) +
  theme(strip.text.y = element_text(size = 8)) -> time.difference.plot


save_plot(time.difference.plot, filename = "exploratory.figs/time.difference.plot.pdf",
          base_height = 8,
          base_aspect_ratio = 1.8)
