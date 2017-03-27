#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)


### Timeseries
analysis.dir <- "../03-27-2017_09-26/"

timeseries <- read.table("../03-27-2017_09-26/out.timeseries.txt", header = TRUE)
hundred.thousand <- read.table("../03-27-2017_01-55/out.timeseries.txt", header = TRUE)

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


timeseries %>% select(date, totalS, totalI) %>%
  gather(key = state, value = value, -date) %>%
  ggplot(aes(x = date, y = value)) + geom_line()+
  facet_wrap(~state, scales = "free_y") +
  scale_y_continuous(labels = fancy_scientific) +
  labs(x = "Date (years)", y = "Individuals")

hundred.thousand %>% select(date, totalS, totalI, diversity) %>%
  gather(key = variable, value = value, -date) %>%
  ggplot(aes(x = date, y = value)) + geom_line()+
  facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(labels = fancy_scientific) +
  labs(x = "Date (years)", y = "Value")


####### Antigenic Change
# one line graph, number of times they occur
# histogram, the size
library(readr)
antigenic.mutations <- read_delim("~/Documents/projects/fluantigen/03-27-2017_09-26/out.console.txt", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
antigenic.mutations$distance = as.numeric(antigenic.mutations$distance)

# Plotting timeseries of antigenic mutations 
antigenic.mutations %>% ggplot(aes(x = day, y = distance)) + 
  geom_point() + facet_wrap(~oriAntigenType)

antigenic.mutations %>% ggplot(aes(x = distance)) + geom_histogram() + geom_vline(xintercept = .012)

antigenic.mutations %>% filter(oriAntigenType == "0") %>%
  ggplot(aes(x = distance)) + geom_histogram() + geom_vline(xintercept = .012)

antigenic.mutations %>% group_by(day) %>%
  summarise(num.mutations = length(day)) %>%
  ggplot(aes(x = day, y = num.mutations)) + geom_line()


##### Plotting frequencies of different antigenic 

out_antigenFrequencies <- read_delim("~/Documents/projects/fluantigen/03-27-2017_09-26/out.antigenFrequencies.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

out_antigenFrequencies$antigentype = as.factor(out_antigenFrequencies$antigentype)
out_antigenFrequencies$frequency = as.numeric(out_antigenFrequencies$frequency)

## Filter out Antigen Frequencies that are greater than 5%
out_antigenFrequencies %>%
  group_by(antigentype) %>%
  summarize(max.freq = max(frequency)) -> max.frequency

quantile(x = max.frequency$max.freq, probs = )
length(which(max.frequency$max.freq < .01))/3321 # Over 99% of novel antigenic mutations don't reach over 10%
%>%
  ggplot(aes(max.freq)) + geom_histogram(bins = 100) + ylim(0, 500) + geom_vline(xintercept = .05)
  
max.frequency %>% 
  filter(max.freq > .05) -> threshold.frequencies

selected.antigens = threshold.frequencies$antigentype

out_antigenFrequencies %>%
  filter(antigentype %in% selected.antigens) %>%
  mutate(year = day/365) %>%
  ggplot(aes(x = year, y = frequency, fill = antigentype)) +
  geom_area() + guides(fill = FALSE) -> ci.frequency

save_plot(filename = "ci.frequency.pdf", ci.frequency)




antigen_series %>% 
  gather(key = variable, value = value, -day) %>%
  ggplot(aes(x = day, y = value)) + geom_line() + 
  facet_wrap(~variable, scales = "free_y") -> summary.plot
save_plot(filename = "summary.plot.pdf", summary.plot, base_height = 8, base_aspect_ratio = 1.5)



#### Variance
out_viralFitnessSeries <- read_delim("~/Documents/projects/fluantigen/03-27-2017_09-26/out.viralFitnessSeries.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

out_viralFitnessSeries %>% 
  gather(key = variable, value = value, -date) %>%
  ggplot(aes(x = date, y = value)) + geom_line() + 
  facet_wrap(~variable, scales = "free_y") -> viralFitnessSeries
save_plot(filename = "summary.plot.pdf", summary.plot, base_height = 8, base_aspect_ratio = 1.5)
