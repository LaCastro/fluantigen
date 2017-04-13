###### Functions

read.outputfiles <- function(dir, type ) {
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    output.file <- read.table(paste0(dir,.file, type), header = TRUE)
    output.file
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}



read.console.files <- function(dir) {
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    output.file <- read.table(paste0(dir,.file, "/out.console.txt"), header = TRUE, fill = TRUE)
    output.file = output.file[1:(nrow(output.file)-2),]
    output.file
  })
  names(population.data) = file.list
  #population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}


find.antigen.emergence <- function(antigen.data){
  antigen.data %>%
    filter(!duplicated(postAntigen)) %>%
    filter(postAntigen != 0) -> novel.types 
  return(novel.types)
}



find.successful.types <- function(frequencies, threshold) {
  frequencies %>%
    group_by(antigentype) %>%
    summarize(max.freq = max(frequency)) %>%
    filter(max.freq > threshold) -> successful.types
    
    successful.types.id = successful.types$antigentype
  
    return(successful.types.id)
}


find.dominant.types.at.emerge <- function(antigen.frequencies) {
  antigen.frequencies %>%
    group_by(day) %>%
    slice(which.max(frequency)) 
}


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}



create.meta.data <- function(sim.dir) {
  ### Combines all the output files for novel antigens
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
  
  colnames(dominant.types)[2] = "dominant.type"; colnames(dominant.types)[3] = "dominant.freq"
  dominant.types$day = as.character(dominant.types$day)
  
  meta.data %>%
    left_join(dominant.types, by = "day") -> meta.data
  
  # Differentiate whether it was sucessful or not
  meta.data$success = NA
  successful.types = find.successful.types(antigen.frequencies, threshold = .15)
  meta.data %>%
    mutate(success = ifelse(postAntigen %in% successful.types, "yes", "no")) -> meta.data
  
  return(meta.data)
}


create.meta.data.all <- function(dir) {
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    meta.data = create.meta.data(paste0(output.folder, sim.dir = .file))
    return(meta.data)
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}

