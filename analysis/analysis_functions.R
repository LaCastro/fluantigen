###### Functions



read.outputfiles <- function(dir, type ) {
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    output.file <- read.table(paste0(analysis.dir,.file, type), header = TRUE)
    output.file
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}



find.successful.types <- function(frequencies, threshold) {
  frequencies %>%
    group_by(antigentype) %>%
    summarize(max.freq = max(frequency)) %>%
    filter(max.freq > threshold) -> successful.types
    
    successful.types.id = successful.types$antigentype
  
    return(successful.types.id)
}