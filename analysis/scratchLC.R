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
analysis.dir <- "~/Dropbox/Projects/mutantigen/04-04-2017_09-20-07/"
analysis.dir <- "~/Documents/projects/fluantigen/04-03-2017_12-32-14/"
variation.list = list.dirs(analysis.dir, full.names = TRUE, recursive = FALSE)


variation.list <- lapply(variation.list, function(.file) {
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


variation.antigen = rbindlist(
  , idcol = TRUE)
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

population.antigenFrequencies <- lapply(population.list, function(.file) {
  out.timeseries <- read.table(paste0(analysis.dir,.file, "/out.antigenFrequencies.txt"), header = TRUE)
  out.timeseries
})

names(population.antigenFrequencies) = population.list
population.antigenFrequencies = rbindlist(population.antigenFrequencies, idcol = TRUE)

population.antigenFrequencies %>% 
  ggplot(aes(x = date, y = totalI)) + geom_line()+
  facet_wrap(~.id, scales = "free_y")


################# Antigen Series


population.timeseries <- lapply(population.list, function(.file) {
  out.timeseries <- read.table(paste0(analysis.dir,.file, "/out.timeseries.txt"), header = TRUE)
  out.timeseries
})

names(population.timeseries) = population.list
population.timeseries = rbindlist(population.timeseries, idcol = TRUE)

#population.timeseries %>% 
#  ggplot(aes(x = day, y = cumulativeTypes)) + geom_line()+
#  facet_wrap(~.id, scales = "free_y")

timeseries  = read.table(paste0(analysis.dir, "/out.timeseries.txt"), header = TRUE)
desired.metrics = c("northI", "antigenicDiversity", "northTmrca", "diversity")

timeseries %>%
  gather(key = metric, value = value, - date) %>%
  filter(metric %in% desired.metrics ) %>%
  ggplot(aes(x = date, y = value)) + geom_line() +
  facet_wrap(~metric, scales = "free")

timeseries %>%
  gather(key = metric, value = value, -date) %>%
  filter(metric == "northTmrca") -> tmrca

northTmrca
antigenFrequencies = read.table(paste0(analysis.dir, "out.antigenFrequencies.txt"), header = TRUE)
antigenFrequencies$antigentype = as.factor(antigenFrequencies$antigentype)




antigenFrequencies %>%
  group_by(antigentype) %>%
  summarize(max.freq = max(frequency)) %>%
  filter(max.freq > .1) -> successful.mutants
  
successful.mutants

successful.mutants.id = successful.mutants$antigentype

myColors <- colorRampPalette(brewer.pal(8, "Accent"))(length(successful.mutants.id))

antigenFrequencies %>%
  filter(antigentype %in% successful.mutants.id) %>%
  ggplot(aes(x = simTime, y = frequency, group = antigentype, col = antigentype)) + 
  geom_line(size = 1.2) + guides(col = FALSE) + 
  scale_color_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") -> frequency.short





population.antigenFrequencies %>% 
  filter(.id == "ten_sixth") %>%
  ggplot(aes(x = simTime, y = frequency, group = antigentype, col = antigentype)) + 
  geom_line() + guides(col = FALSE) + 
  scale_color_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") -> antigenFrequencies.plot
save_plot(antigenFrequencies.plot, filename = "antigenFrequencies.pdf", base_height = 4, base_aspect_ratio = 1.5)


population.antigenFrequencies %>%
  filter(.id == "ten_sixth") %>%
  filter(antigentype %in% successful.mutants.id) %>%
  mutate(antigen.prevalence = round(frequency * infected)) %>%
  ggplot(aes(x = simTime, y = antigen.prevalence)) +
  geom_area(aes(fill = antigentype), stat = "identity", position = "stack", na.rm = TRUE) +
  #geom_line(aes(color = antigentype, fill = antigentype), size = 1.2) + guides(col = FALSE) + 
  scale_color_manual(values = myColors.short) +
  labs(y = "Individuals Infected", x = "Years") +
  guides(fill = FALSE) -> individuals.infected.plot.successful



plot_grid(individuals.infected.plot.successful, frequency.short, ncol = 1)


population.antigenFrequencies %>% 
  filter(.id == "ten_sixth") %>%
  ggplot(aes(x = simTime, y = frequency, group = antigentype, col = antigentype)) + 
  geom_line() + guides(col = FALSE) + 
  scale_color_manual(values = myColors) +
  labs(y = "Frequency", x = "Years") -> antigenFrequencies.plot
save_plot(antigenFrequencies.plot, filename = "antigenFrequencies.pdf", base_height = 4, base_aspect_ratio = 1.5)


#extract days that are 10 and then 
population.antigenFrequencies %>%
  filter(.id == "ten_sixth") %>%
  filter((day %% 10 == 0)) %>%
  mutate(antigen.prevalence = round(frequency * infected)) %>%
  ggplot(aes(x = simTime, y = antigen.prevalence, group = antigentype, col = antigentype)) +
  geom_line() + guides(col = FALSE) + 
  scale_color_manual(values = myColors) +
  labs(y = "Individuals Infected", x = "Years") -> individuals.infected.plot

 
plot_grid(antigenFrequencies.plot, individuals.infected.plot, ncol = 1) 



population.timeseries %>%
  filter(.id == "ten_sixth") %>%
  ggplot(aes(x = date, y = totalI)) + geom_line()+
  scale_y_continuous(labels = fancy_scientific) +
  labs(x = "Years", y = "Infected") -> infected.dynamics


plot_grid(antigenFrequencies.plot, infected.dynamics, ncol = 1) -> ten.sixth.dynamics
save_plot(ten.sixth.dynamics, filename = "ten.sixth.dynamics.pdf", base_height = 8, base_aspect_ratio = 1.2)




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
population.antigenFrequencies$.id <- factor(population.antigenFrequencies$.id, levels = c("ten_third", "ten_fourth", "ten_fifth", "ten_sixth"))

population.antigenFrequencies %>%
  group_by(.id, antigentype) %>%
  summarize(max.freq = max(frequency)) %>%
  ggplot(aes(max.freq)) + geom_histogram(bins = 50)+
  facet_wrap(~.id, scales = c("free")) +
  labs(x = "Maximum Frequency of Antigen") -> antigen.histogram
save_plot(filename = "antigen.histogram.pdf", plot = antigen.histogram, base_height = 8)


population.antigenFrequencies %>%
  filter(.id == "ten_sixth") %>%
  group_by(antigentype) %>%
  summarize(max.freq = max(frequency)) %>%
  ggplot(aes(max.freq)) + geom_histogram(bins = 50) +
  coord_cartesian(ylim = c(0,100)) +
  scale_x_continuous(breaks = seq(0, 1, .1)) -> zoomed.10_6

save_plot(filename = "antigen.histogram106.pdf", plot = zoomed.10_6)


population.timeseries %>%
  filter(.id == "ten_sixth") %>%
  mutate(prevalence = totalI/100,000) %>%
  ggplot(aes(x = date, y = prevalence)) + geom_line()
  


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



### Mean Viral Fitness Series 

viral.fitness = read.table(paste0(analysis.dir, "out.viralFitnessSeries.txt"), header = TRUE)
viral.fitness = viral.fitness[!duplicated(viral.fitness),]


metric.mean = c("meanR", "varR")
metric.beta = c("meanBeta", "varBeta")
metric.sigma = c("meanSigma", "varSigma")

viral.fitness %>%
  gather(key = metric, value = value, -date) %>%
  mutate(value.exp = exp(value)) %>%
  filter(metric %in% metric.mean) %>%
  select(date, metric, value.exp) %>%
  spread(key = metric, value = value.exp) %>%
  ggplot(aes(x = date, y = meanR)) + geom_line() 

viral.fitness %>%
  gather(key = metric, value = value, -date) %>%
  mutate(value.exp = exp(value)) %>%
  filter(metric %in% metric.beta) %>%
  select(date, metric, value.exp) %>%
  spread(key = metric, value = value.exp) -> viral.fitness.beta

viral.fitness %>%
  gather(key = metric, value = value, -date) %>%
  mutate(value.exp = exp(value)) %>%
  filter(metric %in% metric.sigma) %>%
  select(date, metric, value.exp) %>%
  spread(key = metric, value = value.exp) -> viral.fitness.sigma


viral.fitness.beta %>%
  ggplot(aes(x = date, y = meanBeta)) + geom_ribbon(aes(ymin = meanBeta-varBeta, ymax = meanBeta + varBeta), fill = "lightblue") + 
  geom_line(aes(y = meanBeta))

