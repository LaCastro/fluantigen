##### Plotting Functions

set.my.colors <- function(length.of.success) {
  myColors = colorRampPalette(brewer.pal(8, "Accent"))(length.of.success)
  return(myColors)
}


plot.successful.frequency <- function(antigen.frequencies, successful.types) {
  # Plot the frequencies of successful antigens through time 
  myColors <- colorRampPalette(brewer.pal(8, "Accent"))(length(successful.types))
  antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)
  
  antigen.frequencies %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    filter(antigentype %in% successful.types) %>%
    select(-infected) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2) %>%
    mutate(year = day/365) %>%
    ggplot(aes(x = year, y = frequency)) + 
    geom_area(aes(color = antigentype, fill = antigentype)) + 
    guides(col = FALSE) + guides(fill = FALSE) +
    scale_color_manual(values = myColors) +
    scale_fill_manual(values = myColors) +
    labs(y = "Frequency", x = "Years") 
}


plot.successful.infections <- function(antigen.frequencies, successful.types) {
  # Plot the frequencies of successful antigens through time 
  myColors <- colorRampPalette(brewer.pal(8, "Accent"))(length(successful.types))
  antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)
  
  antigen.frequencies %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    filter(antigentype %in% successful.types) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2, -3) %>%
    mutate(year = day/365) %>%
    mutate(prevalence = infected*frequency) %>%
    ggplot(aes(x = year, y = prevalence)) + 
    geom_area(aes(color = antigentype, fill = antigentype)) + 
    #guides(col = FALSE) + guides(fill = FALSE) +
    scale_color_manual(values = myColors) +
    scale_fill_manual(values = myColors) +
    labs(y = "Infecteds", x = "Years") 
}
