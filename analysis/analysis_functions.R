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

