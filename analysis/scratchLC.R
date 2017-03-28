#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)

### Timeseries
analysis.dir <- "~/Dropbox/Projects/mutantigen/"
variation.list = list.dirs(analysis.dir, full.names = TRUE, recursive = FALSE)

all.antigen <- lapply(variation.list, function(.file) {
  file.list = list.files(.file)
  out.timeseries <- read.table(paste0(.file,"/out.trackAntigenSeries.txt"), header = TRUE)
  out.timeseries
})


variation.timeseries = rbindlist(all.timeseries, idcol = TRUE)
variation.timeseries$.id = as.factor(variation.timeseries$.id)

myColors <- colorRampPalette(brewer.pal(8, "Dark2"))(10)

variation.timeseries %>% 
  ggplot(aes(x = date, y = diversity, group = .id, col = .id)) + 
  geom_line(alpha = .5) + guides(col = FALSE) +
  scale_color_manual(values = myColors)


variation.timeseries %>% 
  ggplot(aes(x = date, y = diversity)) +
  facet_wrap(~.id) + geom_line(alpha = .5) 



variation.frequencies = rbindlist(all.frequencies, idcol = TRUE)
variation.frequencies$.id = as.factor(variation.frequencies$.id)

### Visualzing the frequency of the wildtype

variation.frequencies %>%
  filter(antigentype == 0) %>%
  ggplot(aes(x = simTime, y = frequency)) + geom_area() +
  facet_wrap(~.id) + 
  geom_area(color  = "black", fill = "purple", alpha = .5) +
  guides(fill = FALSE) -> wild.type.frequency


variation.fitness  = rbindlist(all.fitness, idcol = TRUE)
variation.fitness$.id = as.factor(variation.fitness$.id)


variation.antigen = rbindlist(all.antigen, idcol = TRUE)
variation.antigen$.id = as.factor(variation.antigen$.id)

variation.antigen %>%
  gather(key = metric, value = value, -day, -.id) %>%
  filter((metric != "R") & (metric != "N")) %>%
  ggplot(aes(x = day, y = value, group = .id, col = .id)) + geom_line() +
  facet_wrap(~metric, scales = "free_y") + guides(color = FALSE)


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


##### Population
analysis.dir <- "~/Dropbox/Projects/mutantigen/population_test/"
population.list = list.dirs(analysis.dir, full.names = FALSE, recursive = FALSE)

population.timeseries <- lapply(population.list, function(.file) {
  out.timeseries <- read.table(paste0(analysis.dir,.file, "/out.timeseries.txt"), header = TRUE)
  out.timeseries
})

names(population.timeseries) = population.list
population.timeseries = rbindlist(population.timeseries, idcol = TRUE)

population.timeseries %>% 
  ggplot(aes(x = date, y = totalI)) + geom_line()+
  facet_wrap(~.id, scales = "free_y")


################# Antigen Series
population.antigenFrequencies <- lapply(population.list, function(.file) {
  out.timeseries <- read.table(paste0(analysis.dir,.file, "/out.antigenFrequencies.txt"), header = TRUE)
  out.timeseries
})

names(population.antigenFrequencies) = population.list
population.antigenFrequencies = rbindlist(population.antigenFrequencies, idcol = TRUE)

population.antigenseries %>% 
  ggplot(aes(x = day, y = cumulativeTypes)) + geom_line()+
  facet_wrap(~.id, scales = "free_y")

population.antigenFrequencies$antigentype = as.factor(population.antigenFrequencies$antigentype)

population.antigenFrequencies %>% 
  filter(.id == "ten_sixth") %>%
  ggplot(aes(x = simTime, y = frequency, group = antigentype, col = antigentype)) + 
  geom_line() + guides(col = FALSE)
  


head(population.antigenseries)


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


############ Library phylotate
tree_sixth <- read_annotated(filename = paste0(analysis.dir,"/ten_sixth/out.trees.txt"), format = "newick")
tree_third <- read_annotated(filename = paste0(analysis.dir,"/ten_third/out.trees.txt"), format = "newick")
tree_fourth <- read_annotated(filename = paste0(analysis.dir,"/ten_fourth/out.trees.txt"), format = "newick")
tree_fifth <- read_annotated(filename = paste0(analysis.dir,"/ten_fifth/out.trees.txt"), format = "newick")

plot(tree_sixth, show.tip.label = FALSE)
