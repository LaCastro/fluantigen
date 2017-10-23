#### If only looking at transient -- should identify which ones those are right away before 
#### collecting the rest of the data on them


# Step 1. Put together list trials of transient and successful antigens ; this will be based on max.freq, and days above
trial.dirs = dir(tropics.folder)
antigen.frequencies = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.antigenFrequencies.txt", header = TRUE)) 
names(antigen.frequencies) = trial.dirs
name(antigen.frequencies[1])

######## Step 2: For each list entry need to calculate max frequency and days above
thres = .2

(antigen.frequencies[[1]])

# Calculating the two criteria for determining success 
antigen.frequencies %>%
  map(find_max_frequency) %>%
  map(function(x) mutate_at(x, "antigentype", as.character)) %>%
  map(function(x) mutate(trial.dir = name(x)-> max.frequencies

max.frequencies[[1]]

antigen.frequencies %>%
  map(function(x) calculate_days_above(x,threshold = .2)) -> days.above

# Creating the combined data frame and using the criteria to assign a label 
full.data = map2(max.frequencies, days.above, left_join) %>%
  map(replace_na_zeros) %>% 
  map(determine_success_labels)
names(full.data) = trial.dirs

####### Step 3: Subsetting for just transient and yes  
subset.data = map(full.data, filter_out_loss) 
subset.df = do.call("rbind", subset.data)
subset.df$name = rep(trial.dirs, sapply(subset.data, nrow))
head(subset.df)
####### Step 4: Create Meta Data for these 

create_meta_data <- function(sim.dir) {
  browser()
  # Combines all the output files for novel antigens
  # reading in and cleaning up the console file 
  console.file <- read.table(paste0(sim.dir, "/out.console.txt"), header = TRUE, fill = TRUE) %>% clean_console()
  
  # find the time when novel antigens emerge
  novel.types = find_antigen_emergence(console.file) %>% mutate_at("day", as.character)
  
  # read in output file that tracks state of simulation when antigens emerge
  track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE) %>%
    mutate_at(.vars = "day", as.character)
  meta.data  = left_join(x = novel.types, y = track.antigen, "day") 
  
  # record the exact days when new novel types are generated 
  days.of.emergence = meta.data$day
  
  ## Read in viral fitness, dplyr::select on the days that correspond to emergence 
  viral.fitness.emergence <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE) %>%
    filter(day %in% days.of.emergence) %>%
    mutate(day = as.character(day)) %>%
    dplyr::select(-simDay) %>%
    filter(!duplicated(day)) 
  
  meta.data %>%
    left_join(viral.fitness.emergence, by = "day") -> meta.data
  
  ## combine dominant frequency circulating at time of emergence
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE) 
  dominant.types = find_dominant_types_at_emerge(antigen.frequencies) %>%
    filter(day %in% days.of.emergence) %>%
    dplyr::rename(dominant.type = antigentype, dominant.freq = frequency) %>%
    ungroup() %>%
    mutate_at(.vars = "day", as.character)
  
  meta.data %>%
    left_join(dominant.types, by = "day") -> meta.data
  
  # combine maximum frequency the strain itself ever achieved
  maximum.freq.type <- find_max_frequency(antigen.frequencies) %>% mutate_at("antigentype", as.character)
  meta.data %>%
    mutate_at("postAntigen", as.character) %>%
    left_join(maximum.freq.type, by = c("postAntigen" = "antigentype")) -> meta.data
  
  life.span <- calculate_total_life(antigen.frequencies) %>% mutate_at("antigentype", as.character)
  meta.data %>% left_join(life.span, by = c("postAntigen" = "antigentype")) -> meta.data
  
  # Differentiate whether it was sucessful or not
  #days.above.thres = _above_thres(antigen.frequencies, success.criteria$threshold)
  #days.above.thres$antigentype = as.character(days.above.thres$antigentype)
  #meta.data %>% left_join(days.above.threshold, by = c("postAntigen" = "antigentype")) -> meta.data
  #meta.data %>%
  #  mutate(success = ifelse((final.max > success.criteria$freq & life.length > success.criteria$length.days), "yes", "no")) -> meta.data
  
  return(meta.data)
}




## Calculates the number of days an antigen type is present above a certain frequency
calculate_days_above <- function(antigen.frequencies, threshold) { 
  antigen.frequencies %>%
    group_by(antigentype) %>%
    dplyr::filter(frequency > threshold) %>%
    summarize(day.emerge = day[1],
            last.day = tail(day)[1],
            frequency.emerge = frequency[1],
            frequency.down = tail(frequency)[1]) %>%
    mutate(days = last.day-day.emerge) %>%
    dplyr::select(antigentype, days) %>%
    mutate_at("antigentype", as.character) 
}
replace_na_zeros = function(x) {
  #browser()
  already_lost = which(is.na(x$days))
  x$days[already_lost] = 0
  return(x)
}
determine_success_labels = function(x) {
  x %>%
    mutate(success = ifelse(days > 45, "Est.", 
                          ifelse(final.max > .1, "Transient", "no")))
}
filter_out_loss = function(x) {
  x %>%
    filter(success != "no")
}