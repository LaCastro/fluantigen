##### Script for combining the different outputs to get one database 

output.folder = "../04-04-2017_03-52-57/"

antigenic.mutations = read.table(paste0(output.folder, "out.console.txt"), header = TRUE, fill = TRUE)
antigenic.mutations = antigenic.mutations[1:(nrow(antigenic.mutations)-2),]

# Function for getting the day when the mutation happened

find.antigen.emergence <- function(antigen.data){
  antigen.data %>%
    filter(!duplicated(postAntigen)) %>%
    filter(postAntigen != 0) -> novel.types 
  return(novel.types)
}

novel.types = find.antigen.emergence(antigenic.mutations)


track.antigen <- read.table(paste0(output.folder, "out.trackAntigenSeries.txt"), header = TRUE)
track.antigen$day = as.character(track.antigen$day)

days.of.emergence = meta.data$day

left_join(x = novel.types, y = track.antigen, "day") -> meta.data


## Read in viral fitness-get days of emergence
viral.fitness <- read.table(paste0(output.folder, "out.viralFitnessSeries.txt"), header = TRUE)
head(viral.fitness)

viral.fitness %>%
  filter(day %in% days.of.emergence) %>%
  mutate(day = as.character(day)) %>%
  select(-simDay) %>%
  filter(!duplicated(day)) -> viral.fitness.emergence

meta.data %>%
  left_join(viral.fitness.emergence, by = "day") -> meta.data


## combine dominant and frequency
antigen.frequencies <- read.table(paste0(output.folder, "out.antigenFrequencies.txt"), header = TRUE)

find.dominant.types.at.emerge <- function(antigen.frequencies) {
  antigen.frequencies %>%
    group_by(day) %>%
    slice(which.max(frequency)) 
}

dominant.types <- find.dominant.types.at.emerge(antigen.frequencies)
dominant.types %>%
  filter(day %in% days.of.emergence) -> dominant.types

colnames(dominant.types)[2] = "dominant.type"
colnames(dominant.types)[3] = "dominant.freq"
dominant.types$day = as.character(dominant.types$day)

meta.data %>%
  left_join(dominant.types, by = "day") -> meta.data

successful.types = find.successful.types(antigen.frequencies, threshold = .1)

successful.types