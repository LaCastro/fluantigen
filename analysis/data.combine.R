##### Script for combining the different outputs to get one database 
source('analysis_functions.R')

library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)


output.folder = "../04-05-2017_09-25-21/"
output.folder = "~/Dropbox/Projects/mutantigen/tropics/"

antigenic.mutations = read.table(paste0(output.folder, "out.console.txt"), header = TRUE, fill = TRUE)
antigenic.mutations = antigenic.mutations[1:(nrow(antigenic.mutations)-2),]

### now need to find NA

#### Turn this into a big function 
## could have 

create.meta.data <- function(sim.dir) {
  
  # reading in the console file 
  console.file <- read.table(paste0(sim.dir, "/out.console.txt"), header = TRUE, fill = TRUE)
  console.file = console.file[1:(nrow(console.file)-2),]
  
  # finding the time when novel antigens emerge
  novel.types = find.antigen.emergence(console.file)
  
  # reading in output file for when antigens emerge
  track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
  track.antigen$day = as.character(track.antigen$day)
  
  meta.data  = left_join(x = novel.types, y = track.antigen, "day") 
  days.of.emergence = meta.data$day
  
  ## Read in viral fitness-get days of emergence
  viral.fitness <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE)
  
  viral.fitness %>%
    filter(day %in% days.of.emergence) %>%
    mutate(day = as.character(day)) %>%
    select(-simDay) %>%
    filter(!duplicated(day)) -> viral.fitness.emergence
  
  meta.data %>%
    left_join(viral.fitness.emergence, by = "day") -> meta.data
  
  ## combine dominant and frequency
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  
  dominant.types <- find.dominant.types.at.emerge(antigen.frequencies)
  dominant.types %>%
    filter(day %in% days.of.emergence) -> dominant.types
  
  colnames(dominant.types)[2] = "dominant.type"
  colnames(dominant.types)[3] = "dominant.freq"
  dominant.types$day = as.character(dominant.types$day)
  
  meta.data %>%
    left_join(dominant.types, by = "day") -> meta.data
  
  # Differentiate whether it was sucessful or not
  meta.data$success = NA
  successful.types = find.successful.types(antigen.frequencies, threshold = .1)
  
  meta.data %>%
    mutate(success = ifelse(postAntigen %in% successful.types, "yes", "no")) -> meta.data
  
  return(meta.data)
}

trial.1 = create.meta.data(sim.dir = paste0(output.folder, "north3"))


novel.types = find.antigen.emergence(antigenic.mutations)
track.antigen <- read.table(paste0(output.folder, "out.trackAntigenSeries.txt"), header = TRUE)
track.antigen$day = as.character(track.antigen$day)


left_join(x = novel.types, y = track.antigen, "day") -> meta.data
days.of.emergence = meta.data$day

## Read in viral fitness-get days of emergence
viral.fitness <- read.table(paste0(output.folder, "out.viralFitnessSeries.txt"), header = TRUE)

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

##### meta.data, differentiate whether it was sucessful or not
meta.data$success = NA

meta.data %>%
  mutate(success = ifelse(postAntigen %in% successful.types, "yes", "no")) -> meta.data
  


