#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)


### Beford Timeseries
bedford.timeseries = read.table("~/Documents/projects/bedford_koelle/antigen/out.timeseries", header = TRUE)
head(bedford.timeseries)
bedford.timeseries %>% ggplot(aes(x = date, y = totalI)) + geom_line()



analysis.dir <- "../build/"

timeseries <- read.table("../out.timeseries.txt", header = TRUE)
timeseries %>% ggplot(aes(x = date, y = totalI)) + geom_line()





desired.variables <- c("diversity", "netau", "antigenicDiversity", "totalI", "totalS")

timeseries %>% gather(key = variable, value = value, -date) %>%
  filter(variable %in% desired.variables) -> timeseries.long

timeseries.long %>% ggplot(aes(x = date, y = value, group = variable)) + 
  geom_line(size = 1) + 
  facet_wrap(~variable, scales = "free")



example.timeseries <- read.table("../../antigen/example/out.timeseries", header = TRUE)
example.timeseries %>% gather(key = variable, value = value, -date) %>%
  filter(variable %in%  desired.variables) %>%
  ggplot(aes(x=date, y = value, group = variable))+
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free")


console %>% gather(key = variable, value = value, -day) %>%
  filter(variable == "diversity") %>%
  ggplot(aes(x = day, y = value, group = value))+ geom_line()



antigen.frequencies <- read.table("../out.antigenFrequencies.txt", header=TRUE)
antigen.frequencies$antigentype = as.factor(antigen.frequencies$antigentype)
antigen.frequencies %>% 
  mutate(individuals.I = round(frequency*infected)) -> antigen.frequencies

antigen.frequencies %>%
  ggplot(aes(x = day, y = individuals.I, group = antigentype, col = antigentype)) +
  geom_line() + guides(color = FALSE)


antigen.population <- read.table("../out.trackAntigenSeries.txt", header = TRUE)




### 
ci.pop.1000 <- ggplot2::ggplot(data = long.strains.absolute.sample,
                               aes(x = vtime, y = num.infected, fill = strain.name)) +
  facet_wrap(~iter, nrow = 3) +
  geom_area(color = "black", size = .2, alpha = .3) + guides(fill = FALSE) +
  #geom_vline(data = time.max.sample, aes(xintercept = jitter(value, .25), color = metric), size = 1.25, alpha = .85) +
  #scale_color_manual(values = c("black", "mediumblue", "purple"), guide = FALSE) + 
  #theme(strip.background = element_blank(),  strip.text.x = element_blank()) +
  labs(x = "Time", y = "Num.Infected")
