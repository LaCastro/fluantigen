#### Successful Mutations

library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(data.table)
library(RColorBrewer)


analysis.dir <- "~/Dropbox/Projects/mutantigen/population_test/"
figure.dir <- "exploratory.figs/"


##### Population
## Step 1 Figure out the antigen types that were successful - this should turn into a function to read in different files

antigenFrequencies = read.outputfiles(analysis.dir, type = "/out.antigenFrequencies.txt")
antigenFrequencies$antigentype = as.factor(antigenFrequencies$antigentype)

antigenFrequencies %>%
  filter(.id == "ten_sixth2") -> antigenFrequencies

successful.types = find.successful.types(antigenFrequencies, threshold = .1)

##### Step 2 Collect Data for time of emergence by individual file
antigen.emergence <- read.outputfiles(analysis.dir, type = "/out.console.txt")

antigen.emergence %>%
  filter(postAntigen %in% successful.types) %>%
  mutate(occurrence = !duplicated(postAntigen)) %>%
  filter(occurrence == "TRUE") -> parent.of.successful




##### For each parent-find out how many antigenic mutations had occurrenced on it without changing antigen type

parents = parent.of.successful$oriAntigen

### This doesn't work--doesn't show how many antigenic mutations have occurred on that specific mutatio

antigen.emergence %>% 
  filter(oriAntigen %in% parents) %>%
  filter(oriAntigen == postAntigen) %>%
  group_by(oriAntigen) %>%
  summarize(previous.mutations = length(distance))


myColors <- colorRampPalette(brewer.pal(8, "Accent"))(3321)
myColors.short <- brewer.pal(n = 7, "Accent")

population.antigenFrequencies %>%
  filter(.id == "ten_sixth") %>%
  filter(antigentype %in% successful.mutants.id) %>%
  ggplot(aes(x = simTime, y = frequency, group = antigentype, col = antigentype)) + 
  geom_line(size = 1.2) + guides(col = FALSE) + 
  scale_color_manual(values = myColors.short) +
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

