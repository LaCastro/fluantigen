rm(list=ls())
source('loadLibraries.R')
source('analysis_functions.R')
source('plotting_functions.R')



############### Precision 
model.performance = read.csv("~/Dropbox/model.performance.csv")
model.performance$model <- factor(model.performance$model, levels = c("freq.1", "freq.3", "freq.5", 
                                                                      "gp.1","gp.2",
                                                                      "full.data", "antigen.5", "antigen.full",
                                                                      "population.5", "population.full", "ratio.5",
                                                                      "ratio.full", "viral.fitness.5", "viralFitness.full"))

levels(model.performance$model.type) = c("All Variables", "Variable Set")
myColors = colorRampPalette(brewer.pal(8, "Accent"))(nrow(model.performance))
model.performance %>% ggplot(aes(sen, spec)) + geom_point(aes(fill = model), size = 3, colour = "black", pch=21) +
  #geom_point(,size = 3,) + 
  labs(x = "True Positve Rate", y = "False Positive Rate", fill = "Model") + 
  scale_y_continuous(breaks = seq(.55, .95, .05)) + scale_x_continuous(breaks = seq(.45,.8,.05))+
  geom_abline(intercept = 0, slope = 1, color = "grey", alpha = .5) + facet_wrap(~model.type) + 
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(values = myColors, labels = c("All Variables: 1%", "All Variables: 3%", "All Variables: 5%",
                                                   "All Variables: Growth Phase 1", "All Variables: Growth Phase 2", 
                                                   "All Variables: all time", "Antigen variables: 5%",
                                                   "Antigen variables: all time", "Population variables: 5%",
                                                   "Population variables: all time", "Relative fitness variables: 5%",
                                                   "Relative fitness variables: all time", "Viral pop. variables: 5%", 
                                                   "Viral pop.variables: all time")) -> model.precision
model.precision
save_plot("exploratory.figs/logistic_regressions/model.precision.pdf", model.precision, base_height = 6, base_aspect_ratio = 1.8)


###################################################
meta.data.variables = read.csv("~/Dropbox/meta.data.variables.csv", header = TRUE)

# Lookinng at terms that were important
##### Plot 1 
# Creates data frame 
meta.data.variables %>%
  filter(independent == "yes") %>%
  group_by(variable, model.inclusion) %>%
  summarize(n = n(), variable.type = variable.type[1]) %>%
  mutate(freq = n/sum(n),
         total.opp = sum(n)) %>%
  filter(model.inclusion == "yes") %>%
  arrange(desc(freq)) -> prop.included

# Creates plot 
variable.order = prop.included$variable[rev(order(prop.included$freq))]
variable.type.order = c("relative.fitness", "antigen.level", "population.level", "viral.pop", "Other")
prop.included$variable <- factor(prop.included$variable, levels = variable.order)  
variable.colors = brewer.pal(n = length(unique(prop.included$variable.type)), name = "Paired")
prop.included$variable.type = factor(prop.included$variable.type, levels = variable.type.order)  


ggplot(prop.included, aes(x = variable, y = freq)) + geom_col(aes(fill = variable.type)) + 
  geom_text(aes(label = total.opp)) + 
  labs(x = "Variable", y = "Proportion Included in Model", fill = "Variable Type") + 
  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1)) +
  scale_fill_manual(values = variable.colors) +
  theme(legend.text = element_text(size = 10)) -> prop.included.plot

# Plot 2 Graph
meta.data.variables$variable.type = factor(meta.data.variables$variable.type, levels = variable.type.order) 
description.colors = brewer.pal(n = length(unique(meta.data.variables$description)), name = "Dark2")

meta.data.variables$description = factor(meta.data.variables$description, levels = c("freq.1", "freq.3", "freq.5",
                                                                                     "gp.1", "gp.2", "accel"))
meta.data.variables %>%
  filter(independent == "yes" & model.inclusion == "yes") %>%
  group_by(variable, description) %>%
  summarize(n = n(), variable.type = variable.type[1]) %>%
  mutate(freq = n/sum(n)) %>%
  ggplot(aes(x = variable, y = freq, fill = description)) +
  facet_wrap(~variable.type, scales = "free") + 
  geom_col() + scale_fill_manual(values = description.colors, 
                                 labels = c("Freq 1%", "Freq 3%", "Freq 5%", "EGP 1", "EGP 2", "Diff EGP")) +
  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  labs( x = "Variable", y = "Count", fill = "Data Time Point") -> variable.time.point.plot


# Combine Plots 
important.variables.plot = plot_grid(prop.included.plot, variable.time.point.plot, rel_widths = c(1, 1)) 
save_plot(plot = important.variables.plot,
          filename = "exploratory.figs/logistic_regressions/included.variables.pdf", base_height = 8,
          base_aspect_ratio = 2)

####### Relative Plot I'm not Using 
#meta.data.variables %>%
#  filter(independent == "yes") %>%
#  group_by(variable, model.inclusion) %>%
#  summarize(n = n(),
#            variable.type = variable.type[1]) %>%
#  mutate(freq = n/sum(n)) -> prop.included.relative
#prop.included.relative$variable <- factor(prop.included.relative$variable, 
#                                          levels = variable.order)  
#prop.included.relative$variable.type = factor(prop.included.relative$variable.type, 
#                                              levels = c("relative.fitness", "antigen.level", "population.level", 
#                                                         "viral.pop", "Other"))  prop.included.relative %>%
#  ggplot(aes(x = variable, y = n,  fill = model.inclusion)) + geom_col() +
#  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1)) +
#  labs( x = "Variable", y = "Count") +
#  scale_fill_manual(values = c("#330066", "#009933")) + 
#  facet_wrap(~variable.type, scales = "free") -> variable.included.count.plot


####################################### Model Rank
combined.models = read.csv("~/Dropbox/combined.models.csv")

unique.variables = combined.models %>% filter(term.type == "single")  %>% 
  filter(model.type == "full.set") %>%
  summarize(n = length(unique(variable)))

combined.models %>% 
  filter(term.type == "single") %>%
  filter(model.type == "full.set" & model != "full.model") -> freq.models


my.colors = colorRampPalette(brewer.pal(8, "YlOrRd"))(17)

freq.models$variable = factor(freq.models$variable, 
                              levels = rev(c("ratio.meanR", "ratio.mutation", "ratio.varR", "ratio.varSigma",
                                             "individual.varSigma", "individual.meanSigma",
                                             "varR", "meanBeta", "covBetaSigma", 
                                              "tmrca", "infected", "entropy",  "day", "dominant.freq", "diversity", "antigenicDiversity")))
freq.models %>% 
  ggplot(aes(x = model, y = variable, fill = as.factor(term))) + geom_tile() +
  scale_fill_manual(values = rev(my.colors)) +
  labs(x = "Model", y = "Variable", fill = "Term", title = "Snapshot Models") +
  scale_x_discrete(labels = c("Freq 1%", "Freq 3%", "Freq 5%", "EGP 1", "EGP 2")) +
  geom_text(aes(label = term)) +
  guides(fill = FALSE) -> freq.models.plot

combined.models %>% 
  filter(term.type == "single") %>%
  filter(model.type == "full.set" & model == "full.model") -> full.model 
full.model$variable = factor(full.model$variable, levels = rev(c("ratio.meanR", "ratio.mutation", "ratio.varR", "ratio.varSigma",
                                                             "individual.varSigma", "varR", "varBeta", "tmrca", "infected",
                                                             "entropy", "day")))
full.model$descriptor = factor(full.model$descriptor, levels = c("freq.1", "freq.3", "freq.5", "gp.2", "accel"))
full.model %>%
  ggplot(aes(x = descriptor, y = variable, fill = as.factor(term))) + 
  geom_tile() +
  scale_fill_manual(values = rev(my.colors)) + 
  labs(x = "Time", y = "", fill = "Term", title = "Full Data Model") +
  geom_text(aes(label = term)) +
  scale_x_discrete(labels = c("Freq 1%", "Freq 3%", "Freq 5%",
                              "EGP 1", "EGP 2", "Diff EGP")) -> full.model.plot


model.rank.plot = plot_grid(freq.models.plot, full.model.plot, rel_widths = c(1,1.1))
save_plot(filename = "exploratory.figs/logistic_regressions/model.rank.plot.pdf", model.rank.plot,
          base_height = 8, base_aspect_ratio = 1.8)



######################################## 
# Variable Zoom

full.model %>%
  arrange(term) 

freq.1.subset %>%
  select(success, name,  antigentype, ratio.mutation, infected,ratio.varSigma) %>%
  mutate_if(is.numeric, scale) %>%
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq1") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.1.variables
freq.3.subset %>%
  select(success, name, antigentype, individual.varSigma) %>%
  mutate_if(is.numeric, scale) %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq3") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.3.variables
freq.5.subset %>%
  select(success, name, antigentype, ratio.meanR, varR, ratio.mutation, ratio.varR, tmrca) %>%
  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq5") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.5.variables
growth.2.subset %>%
  select(success, name, antigentype, entropy, day, varR, varBeta, ratio.meanR, day) %>%
  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "gp.2") %>%
  mutate(id = paste0(variable, "_", time.point)) -> gp.2variables
accelerate %>%
  select(success,name, antigentype, ratio.meanR, day) %>%
  mutate_if(is.numeric, scale)   %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "accel") %>%
  mutate(id = paste0(variable, "_", time.point)) -> accel.variables


all.variables = rbind(freq.1.variables, freq.3.variables, freq.5.variables,
                      gp.2variables, accel.variables)
all.variables$id = factor(all.variables$id, levels = c("ratio.meanR_freq5", "varR_freq5", "ratio.meanR_gp.2", "ratio.mutation_freq5",
                                                       "entropy_gp.2", "day_gp.2", "ratio.varR_freq5", "individual.varSigma_freq3",
                                                       "ratio.meanR_accel", "ratio.mutation_freq1", "infected_freq1", "ratio.varSigma_freq1",
                                                       "varR_gp.2", "varBeta_gp.2", "tmrca_freq5", "day_accel"),
                          labels = c("ratio.meanR: 5%", "varR: 5%", "ratio.meanR: EGP 2", "ratio.mutation: 5%",
                                     "entropy: EGP 2", "day EGP 2", "ratio.varR: 5%", "individual.varSigma: 3%", "ratio.meanR: Diff EGP",
                                     "ratio.mutation: 1%", "infected: 1%", "ratio.varSigma: 1%", "varR: EGP 2", "varBeta: EGP 2", "tmrca: 5%",
                                     "day: Diff EGP"))
all.variables %>% 
  ggplot(aes(id, value, color = success)) + geom_boxplot() +
  scale_y_log10() + 
  scale_color_manual(values = c("orange", "purple")) + 
  labs(x = "", y = "Scaled/Centralized Value") + 
  theme(axis.text.x=element_blank()) +
  theme(legend.position = "top") -> boxplot.variables


full.model %>%
  mutate(id = paste0(variable, "_", descriptor)) -> full.model
unique(full.model$id)
full.model$id = factor(full.model$id, levels = c("ratio.meanR_freq.5", "varR_freq.5", "ratio.meanR_gp.2", "ratio.mutation_freq.5",
                                                       "entropy_gp.2", "day_gp.2", "ratio.varR_freq.5", "individual.varSigma_freq.3",
                                                       "ratio.meanR_accel", "ratio.mutation_freq.1", "infected_freq.1", "ratio.varSigma_freq.1",
                                                       "varR_gp.2", "varBeta_gp.2", "tmrca_freq.5", "day_accel"),
                          labels = c("ratio.meanR: 5%", "varR: 5%", "ratio.meanR: EGP 2", "ratio.mutation: 5%",
                                     "entropy: EGP 2", "day EGP 2", "ratio.varR: 5%", "individual.varSigma: 3%", "ratio.meanR: Diff EGP",
                                     "ratio.mutation: 1%", "infected: 1%", "ratio.varSigma: 1%", "varR: EGP 2", "varBeta: EGP 2", "tmrca: 5%",
                                     "day: Diff EGP"))

full.model.summary.plot = plot_grid(boxplot.variables, log.odds.plot, ncol = 1,
                                    align = 'v', axis = "r")
save_plot(full.model.summary.plot, filename = "exploratory.figs/logistic_regressions/full.model.summary.pdf",
          base_height = 8, base_aspect_ratio = 1.5)


full.model %>%
  ggplot(aes(x = id, y = estimate)) + geom_point() +
  geom_errorbar(aes(ymin = estimate -std.error, ymax = estimate+std.error)) +
  geom_hline(yintercept = 0, color = "grey") +
  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1))  + 
  labs(x = "Variable", y = "Log Odds ") -> log.odds.plot
  
  
#full.model %>%
#  ggplot(aes(x = id, y =exp(estimate))) + geom_point() +
#  geom_errorbar(aes(ymin = exp(estimate -std.error), ymax = exp(estimate+std.error))) +
#  geom_hline(yintercept = 1, color = "grey") +
#  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1))  + 
#  labs(x = "Variable", y = "Odds Ratio") -> odds.ratio.plot


######################## Interactions
combined.models %>% 
  filter(term.type == "int") %>%
  filter(model.type == "full.set" & model == "full.model")

int.one = c("ratio.varR: 5%", "ratio.mutation: 1%")
int.two = c("day EGP 2", "varR: EGP 2")
int.three = c("ratio.varR: 5%", "tmrca: 5%")
int.four = c("ratio.meanR: EGP 2", "ratio.varR: 5%")
int.five = c("ratio.meanR: Diff EGP", "day: Diff EGP")

all.variables %>%
  filter(id %in% int.one) %>%
  dplyr::select(success, name, antigentype, value, id) %>%
  #group_by(success, antigentype) %>% 
  spread(key = id, value = value) %>% 
  ggplot(aes(`ratio.varR: 5%`, `ratio.mutation: 1%`, color = success)) + 
  geom_smooth(method = "lm") + 
  geom_point(alpha = .8) + scale_color_manual(values = c("orange", "purple")) +
  guides(color = FALSE) + 
  annotate("text", label = "L.Os = 0.1884450", x = 4, y = 2, color = "black" )-> int.1.plot

all.variables %>%
  filter(id %in% int.two) %>%
  filter(!is.na(value)) %>% 
  dplyr::select(success, name, antigentype, value, id) %>%
  #group_by(success, antigentype) %>% 
  spread(key = id, value = value) %>% 
  ggplot(aes(`day EGP 2`, `varR: EGP 2`, color = success)) + 
  geom_point(alpha = .8) + geom_jitter() +  scale_color_manual(values = c("orange", "purple")) +
  geom_smooth(method = "lm") + 
  guides(color = FALSE) + 
  annotate("text", label = "L.Os = 0.3685605", x = 6, y = 3, color = "black" )-> int.2.plot

int.three = c("ratio.varR: 5%", "tmrca: 5%")
all.variables %>%
  filter(id %in% int.three) %>%
  filter(!is.na(value)) %>% 
  dplyr::select(success, name, antigentype, value, id) %>%
  #group_by(success, antigentype) %>% 
  spread(key = id, value = value) %>% 
  ggplot(aes(`ratio.varR: 5%`, `tmrca: 5%`, color = success)) + 
  geom_smooth(method = "lm") + 
  geom_point(alpha = .8) + geom_jitter() +  scale_color_manual(values = c("orange", "purple")) +
  guides(color = FALSE) + 
  annotate("text", label = "L.Os = 0.2018307", x = 3, y = 3, color = "black" )-> int.3.plot

int.four = c("ratio.meanR: EGP 2", "ratio.varR: 5%")
all.variables %>%
  filter(id %in% int.four) %>%
  filter(!is.na(value)) %>% 
  dplyr::select(success, name, antigentype, value, id) %>%
  #group_by(success, antigentype) %>% 
  spread(key = id, value = value) %>% 
  ggplot(aes(`ratio.meanR: EGP 2`, `ratio.varR: 5%`, color = success)) + 
  geom_smooth(method = "lm") + 
  geom_point(alpha = .8) + geom_jitter() +  scale_color_manual(values = c("orange", "purple")) +
  guides(color = FALSE) + 
  annotate("text", label = "L.Os = 0.1282845", x = 0, y = 5.5, color = "black" )-> int.4.plot


int.five = c("ratio.meanR: Diff EGP", "day: Diff EGP")
all.variables %>%
  filter(id %in% int.five) %>%
  filter(!is.na(value)) %>% 
  dplyr::select(success, name, antigentype, value, id) %>%
  #group_by(success, antigentype) %>% 
  spread(key = id, value = value) %>% 
  ggplot(aes(`ratio.meanR: Diff EGP`, `day: Diff EGP`, color = success)) + 
  geom_smooth(method = "lm") + 
  geom_point(alpha = .8) + geom_jitter() +  scale_color_manual(values = c("orange", "purple")) +
  annotate("text", label = "L.Os = 0.1485334", x = 3, y = 4, color = "black" )-> int.5.plot

interaction.plot = plot_grid(int.1.plot, int.2.plot, int.3.plot, int.4.plot, int.5.plot,
                             rel_widths = c(1, 1.3, 1))
save_plot(filename = "exploratory.figs/logistic_regressions/interaction.plot.pdf", plot = interaction.plot,
          base_height = 8, base_aspect_ratio = 1.5)




########## Need some sort of plotting function so that it takes in two variables and plots the probability 
as.factor( as.numeric(cut(dataPurged$ratio.meanR,20)))
variable.1 = "ratio.meanR"
variable.2 = "varR"

dataPurged$ratio.meanR = cut(dataPurged$ratio.meanR, breaks = seq(from = min(dataPurged$ratio.meanR),
                                                                  to = max(dataPurged$ratio.meanR), length.out = 15))
dataPurged$varR = cut(dataPurged$varR, breaks = seq(from = min(dataPurged$varR),
                                                           to = max(dataPurged$varR), length.out = 15))
dataPurged %>%
  ggplot(aes(ratio.meanR, varR, color = success)) + geom_point() -> scatter.plot



dataPurged %>%
  select(ratio.meanR, varR, success) %>%
 # filter(ratio.meanR == "(-0.496,-0.0653]") %>%
  group_by(ratio.meanR, varR, success) %>%
  summarize(num = n()) -> trial.bin

trial.bin %>%
  group_by(ratio.meanR, varR) %>%
  mutate(rel.freq = num/sum(num)) -> trial.bin
x.breaks = round(seq(-1.36, 4.24, length.out = 15), digits = 1)
y.breaks = round(seq(-1.9, 8.99, length.out = 15), digits = 1)

trial.bin %>%
  spread(key = success, value = rel.freq, fill = 0) %>%
  gather(key = success, value = rel.freq, Transient, Est.) %>%
  filter(success=="Est.") %>%
  ggplot(aes(ratio.meanR, varR, fill = rel.freq)) + geom_tile() +
  scale_x_discrete(labels = my.breaks) + scale_y_discrete(labels = y.breaks) +
  scale_fill_viridis(option="viridis") -> heat.map

plot_grid(scatter.plot, heat.map)


