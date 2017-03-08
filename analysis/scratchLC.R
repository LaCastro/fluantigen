#####
library(plyr)
library(tidyverse)
library(ggplot2)
library(cowplot)

analysis.dir <- "../build/"

timeseries <- read.table(paste0(analysis.dir, "out.timeseries"), header = TRUE)

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
