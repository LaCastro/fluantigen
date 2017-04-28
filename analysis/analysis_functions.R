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
    # Clean console file here 
    output.file = clean.console(output.file)
    output.file
  })
  names(population.data) = file.list
  #population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}


clean.console <- function(console.file) {
  na.rows = output.file[is.na(output.file$postAntigen),]
  na.indicies = which(is.na(output.file$postAntigen))
  modified.rows = output.file[na.indices-1,]
  
  for(i in 1:length(na.indices)) {
    index = na.indices[i]
    output.file[index, 2:6] = unname(na.rows[i,1:5])
    output.file[index, 1] = modified.rows[i,6]
    distance =  na.rows[i, 3]
    output.file[index, 4] = as.character(distance)
    output.file[index, 1] = modified.rows[i,6]
  }
  output.file = output.file[-(na.indicies-1), ]
  return(output.file)
}


find.antigen.emergence <- function(antigen.data){
  antigen.data %>%
    filter(!duplicated(postAntigen)) %>%
    filter(postAntigen != 0) -> novel.types 
  return(novel.types)
}



find.successful.types <- function(frequencies, threshold, length) {
  frequencies %>%
    filter(frequency > threshold) %>%
    group_by(antigentype) %>%
    summarize(n = n(), 
              max.freq  = max(frequency)) %>%
    filter(n > length) -> successful.types
    successful.types.id = successful.types$antigentype
  
    return(successful.types.id)
}



find.dominant.types.at.emerge <- function(antigen.frequencies) {
  antigen.frequencies %>%
    group_by(day) %>%
    slice(which.max(frequency)) 
}


find.max.frequency <- function(antigen.frequencies) {
  antigen.frequencies %>%
    group_by(antigentype) %>%
    summarize(final.max = max(frequency))
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

  # reading in and cleaning up the console file 
  console.file <- read.table(paste0(sim.dir, "/out.console.txt"), header = TRUE, fill = TRUE)
  console.file = console.file[1:(nrow(console.file)-2),]
  
  
  # find the time when novel antigens emerge
  novel.types = find.antigen.emergence(console.file)
  
  # read in output file that tracks state of simulation when antigens emerge
  track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
  track.antigen$day = as.character(track.antigen$day)
  meta.data  = left_join(x = novel.types, y = track.antigen, "day") 
  
  # record the exact days when new novel types are generated 
  days.of.emergence = meta.data$day
  
  ## Read in viral fitness, select on the days that correspond to emergence 
  viral.fitness <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE)
  viral.fitness %>%
    filter(day %in% days.of.emergence) %>%
    mutate(day = as.character(day)) %>%
    select(-simDay) %>%
    filter(!duplicated(day)) -> viral.fitness.emergence
  meta.data %>%
    left_join(viral.fitness.emergence, by = "day") -> meta.data
  
  ## combine dominant frequency circulating at time of emergence
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  dominant.types <- find.dominant.types.at.emerge(antigen.frequencies)
  
  dominant.types %>%
    filter(day %in% days.of.emergence) -> dominant.types
  colnames(dominant.types)[2] = "dominant.type"; colnames(dominant.types)[3] = "dominant.freq"
  dominant.types$day = as.character(dominant.types$day)
  meta.data %>%
    left_join(dominant.types, by = "day") -> meta.data
  
  # combine maximum frequency the strain itself ever achieved
  maximum.freq.type <- find.max.frequency(antigen.frequencies)
  meta.data %>% left_join(maximum.freq.type, by = c("postAntigen" = "antigentype")) -> meta.data
  
  life.span <- calculate.total.life(antigen.frequencies)
  meta.data %>% left_join(life.span, by = c("postAntigen" = "antigentype")) -> meta.data
  # Differentiate whether it was sucessful or not
  #meta.data$success = NA
  #successful.types = find.successful.types(antigen.frequencies, threshold = .1, length = 180)
  
  
  meta.data %>%
    mutate(success = ifelse((final.max > .1 & life.length > 180), "yes", "no")) -> meta.data
  
  return(meta.data)
}


create.meta.data.all <- function(dir) {
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    meta.data = create.meta.data(paste0(dir, sim.dir = .file))
    return(meta.data)
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}



count.n.antigens <- function(antigen.record) {
  antigen.record %>%
    group_by(.id) %>%
    summarize(unique.antigens = n_distinct(antigentype)) -> n.antigens
}

count.success.antigens <- function(antigen.record, threshold) {
  antigen.record %>%
    group_by(.id, antigentype) %>%
    summarize(max.freq = max(frequency)) %>%
    ungroup() %>%
    group_by(.id) %>%
    summarize(above.threshold = sum(max.freq > threshold)) 
}

calculate.success.rate <- function(n.unique.antigens, n.success.antigens) {
  success.rate = n.success.antigens$above.threshold/n.unique.antigens$unique.antigens
  quantile(success.rate, probs = c(0.025, .5, 0.975))*100
}

antigen.success.summary <- function(threshold.levels, antigen.frequencies) {
  
  n.unique.antigens = count.n.antigens(antigen.frequencies)
  
  antigen.success.quantile = adply(.data = threshold.levels, .margins = 1, function(thres) {
    success.antigens = count.success.antigens(antigen.frequencies, threshold = thres)
    quantiles = calculate.success.rate(n.unique.antigens, n.success.antigens = success.antigens)
  })
  num.antigen.strains = adply(.data = threshold.levels, .margins = 1, function(thres) {
    success.antigens = count.success.antigens(antigen.frequencies, threshold = thres)
    success.antigens %>%
      summarize(min = min(above.threshold),
                mean = mean(above.threshold),
                max = max(above.threshold))
  })
  antigen.success.summary.df = left_join(antigen.success.quantile, num.antigen.strains)
  antigen.success.summary.df$X1 = threshold.levels 
  colnames(antigen.success.summary.df)[1] = "threshold"
  return(antigen.success.summary.df)
}



calculate.life.spans.success <- function(ant.freq.success.l) {
  ant.freq.success.l %>%
    group_by(.id, antigentype) %>%
    filter(frequency > .1) %>%
    summarize(n =n())%>%
    mutate(years = n/365)
}

calculate.total.life <- function(antigen.frequencies) {
  antigen.frequencies %>%
    group_by(antigentype) %>%
    summarize(day.emerge = day[1],
              last.day = tail(day)[6]) -> birth.death.days 
   
  birth.death.days %>% 
    mutate(life.length = ifelse(is.na(last.day), 0, last.day-day.emerge)) %>%
    select(antigentype, life.length)
}