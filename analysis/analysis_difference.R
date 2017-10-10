#### Function to calculate the difference between two time points 

## Two Emerge First Trial -- go by sample and ifirst postAntigen
get_surveillance_difference <- function(trial.id, selected.antigen, time.point1, time.point2, data.l, variable.set) {
 # browser()
   data.l %>%
    filter(.id == trial.id & postAntigen == selected.antigen) -> sample.1
# read in and extract from the data set 
  viralFitness <- read.table(paste0(data.folder, trial.id, "/out.viralFitnessSeries.txt"), header = TRUE)
  trackAntigen <- read.table(paste0(data.folder, trial.id, "/out.trackAntigenSeries.txt"), header = TRUE)
  timeSeries <- read.table(paste0(data.folder, trial.id, "/out.timeSeries.txt"), header = TRUE)
  
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
  
  #Cutting the timeseries to what I wanted 
  variables %>%
    mutate(day = as.numeric(day)) %>%
    filter(day > as.numeric(sample.1[, time.point1]) & day < as.numeric(sample.1[,time.point2])) %>%
    select(eval(col.names)) %>%
    slice(c(1,n())) -> variables.sub
  if(nrow(variables.sub) == 0) {
    return()
  } else {

  difference = as.data.frame(t(apply(variables.sub, 2, diff)))
  difference = cbind(selected.antigen, trial.id, time.point2, difference)
  
  return(difference)
  }
}
get_surveillance_sd <- function(trial.id, selected.antigen, time.point1, time.point2, data.l, variable.set) {
  # browser()
  data.l %>%
    filter(.id == trial.id & postAntigen == selected.antigen) -> sample.1
  # read in and extract from the timeseries 
  viralFitness <- read.table(paste0(data.folder, trial.id, "/out.viralFitnessSeries.txt"), header = TRUE)
  trackAntigen <- read.table(paste0(data.folder, trial.id, "/out.trackAntigenSeries.txt"), header = TRUE)
  timeSeries <- read.table(paste0(data.folder, trial.id, "/out.timeSeries.txt"), header = TRUE)
  
  min.I = min(timeSeries$totalI); max.I = max(timeSeries$totalI)
  
  timeSeries %>%
    mutate(day = date*365,
           normalize.I =(totalI-min.I)/(max.I-min.I)) -> timeSeries.trial
  
  if(variable.set == "timeSeries") {
    variables = timeSeries.trial
    col.names = c("totalS","netau", "normalize.I")
   variables$netau[!is.finite(variables$netau)] <- NA
  } else if(variable.set == "trackAntigen") {
    variables = trackAntigen
    col.names = c( "antigenicDiversity", "antigenicTypes", "diversity", "meanLoad", "tmrca")
  } else if(variable.set == "viralFitness") {
    variables = viralFitness
    col.names = c("covBetaSigma", "meanBeta", "meanR", "meanSigma",
                  "varBeta", "varR", "varSigma")
  }
  
  #Cutting the timeseries to what I wanted 
  variables %>%
    mutate(day = as.numeric(day)) %>%
    filter(day > as.numeric(sample.1[, time.point1]) & day < as.numeric(sample.1[,time.point2])) %>%
    select(eval(col.names)) -> variables.sub
    
  sd = as.data.frame(t(apply(variables.sub, 2, function(x) sd(x, na.rm=TRUE))))
  sd = cbind(selected.antigen, trial.id, time.point2,sd)
  return(sd)
}

### need to do this for each antige for each trial

antigen.data %>%
  filter(success == "Est." | success == "Transient") -> antigens.analyze

data.l.trials = unique(data.l$.id)
antigens.analyze %>%
  filter(.id %in% data.l.trials) -> antigens.analyze
time.point1 = "surv_0.02"
time.point2 = "surv_0.05"

time.Series.change.dataset2 = ddply(.data = antigens.analyze, .variables = ".id", function(trial) {
  trial.id = unique(trial$.id)
  postAntigen.list = unique(trial$postAntigen)
  postAntigen.list= as.numeric(as.character(postAntigen.list))
  
  diff.trial = adply(.data = postAntigen.list,.margins = 1, .id = NULL, function(selected.antigen) {
    print(paste0(trial.id, selected.antigen))
    diff.data = get_surveillance_difference(trial.id, selected.antigen, time.point1, time.point2, data.l, variable.set = "timeSeries")
    if(is.null(diff.data)) {
      return()
    } else {
      sd.data = get_surveillance_sd(trial.id, selected.antigen, time.point1, time.point2, data.l, "timeSeries")
      meta.data = rbind(data.frame(metric = "diff", diff.data),
                        data.frame(metric = "sd", sd.data))
      return(meta.data)
    }})
  print(trial.id)
  return(diff.trial)
  })
trackAntigen.change.dataset2 = ddply(.data = antigens.analyze, .variables = ".id", function(trial) {
  trial.id = unique(trial$.id)
  postAntigen.list = unique(trial$postAntigen)
  postAntigen.list= as.numeric(as.character(postAntigen.list))
  
  diff.trial = adply(.data = postAntigen.list,.margins = 1, .id = NULL, function(selected.antigen) {
    print(paste0(trial.id, selected.antigen))
    diff.data = get_surveillance_difference(trial.id, selected.antigen, time.point1, time.point2, data.l, variable.set = "trackAntigen")
    if(is.null(diff.data)) {
      return()
    } else {
      sd.data = get_surveillance_sd(trial.id, selected.antigen, time.point1, time.point2, data.l, "trackAntigen")
      meta.data = rbind(data.frame(metric = "diff", diff.data),
                        data.frame(metric = "sd", sd.data))
      return(meta.data)
    }})
  print(trial.id)
  return(diff.trial)
})
viralFitness.change.dataset2 = ddply(.data = antigens.analyze, .variables = ".id", function(trial) {
  trial.id = unique(trial$.id)
  postAntigen.list = unique(trial$postAntigen)
  postAntigen.list= as.numeric(as.character(postAntigen.list))
  
  diff.trial = adply(.data = postAntigen.list,.margins = 1, .id = NULL, function(selected.antigen) {
    print(paste0(trial.id, selected.antigen))
    diff.data = get_surveillance_difference(trial.id, selected.antigen, time.point1, time.point2, data.l, variable.set = "viralFitness")
    if(is.null(diff.data)) {
      return()
    } else {
      sd.data = get_surveillance_sd(trial.id, selected.antigen, time.point1, time.point2, data.l, "viralFitness")
      meta.data = rbind(data.frame(metric = "diff", diff.data),
                        data.frame(metric = "sd", sd.data))
      return(meta.data)
    }})
  print(trial.id)
  return(diff.trial)
})


time.series.columns = c(".id", "metric", "selected.antigen", "time.point2", "date", "diversity", "tmrca", 
                        "netau", "serialInterval", "antigenicDiversity", "totalS", "totalI", "totalCases",
                        "normalize.I", "day.y", "surv_0.02", "surv_0.05", "success")

mydivide <- function(x){x[1,]/x[2,]}


time.Series.change.dataset %>%
  left_join(data.l, by = c(".id" = ".id", "selected.antigen" = "postAntigen")) %>%
  left_join(trackAntigen.change.dataset) %>%
  left_join(viralFitness.change.dataset)  %>%
  select(-trial.id) %>%
  group_by(.id, selected.antigen) %>%
  gather(key = variable, value = value, -.id, -metric, -selected.antigen, -surv_0.02, -surv_0.05,-time.point2,-success,-day) %>%
  spread(key = metric, value = value) %>%
  mutate(diff = as.numeric(diff),
         sd = as.numeric(sd),
         ratio = diff/sd) -> first.growth


growth.data = rbind(data.frame(phase = "first", first.growth),
                    data.frame(phase = "second", second.growth))

growth.data %>%
  arrange(selected.antigen) %>%
  mutate(var.id = paste0(phase, "_", variable)) %>%
  gather(key = summary, value = value, diff, sd, ratio) %>%
  mutate(var.id = paste0(var.id, "_", summary)) %>%
  select(.id, selected.antigen, success, var.id, value) -> growth.data.l
