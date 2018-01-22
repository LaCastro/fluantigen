rm(list=ls())
library(viridis)
source('loadLibraries.R')
source('analysis_functions.R')
source('plotting_functions.R')

##################################
two.d.plot <- function(data.set, variable.1, variable.2, num.bins, digits) {
  
  variables = c(variable.1, variable.2)
  
  all.variables %>%
    filter(id %in%variables) %>%
    select(success,name,antigentype,id,value) %>%
    spread(key = id, value = value) -> data.set
  
  
  var.1.breaks = round(seq(from = min(data.set[,`variable.1`]),to = max(data.set[,`variable.1`]), length.out = num.bins), digits = digits)
  var.2.breaks = round(seq(from = min(data.set[,`variable.2`]),to = max(data.set[,`variable.2`]), length.out = num.bins), digits = digits)
  
  
  data.set %>%
    select(success, one_of(variables)) -> data.sub
  
  data.sub %>%
    mutate(variable.1.group = cut(data.set[,`variable.1`], breaks = var.1.breaks),
           variable.2.group = cut(data.set[, `variable.2`], breaks = var.2.breaks)) %>%
    na.omit() %>%
    group_by(variable.1.group, variable.2.group, success) %>%
    summarize(num = n()) %>%
    mutate(rel.freq = num/sum(num)) %>%
    ungroup() %>%
    select(-num) %>%
    group_by(variable.1.group, variable.2.group) %>%
    spread(key = success, value = rel.freq, fill = 0) %>%
    gather(key = success, value = rel.freq, Transient, Est.) %>%
    filter(success == "Est.") %>%
    ggplot(aes(variable.1.group, variable.2.group, fill = rel.freq)) + geom_tile() +
    scale_x_discrete(labels = var.1.breaks) + scale_y_discrete(labels = var.2.breaks) +
    scale_fill_viridis(option="plasma")
}
clean.2D.plot <- function(plot, variable.name.1, variable.name.2, guide = TRUE) { 
  step1 = plot + labs(x = variable.name.1, y  = variable.name.2, fill = "Persistent \n Antigens \n") +
    theme(axis.text = element_text(size =14), axis.text.x = element_text(angle = 60,  hjust=1),
          axis.title = element_text(size = 16)) + 
    theme(legend.key.size = unit(0.7, "cm")) + theme(legend.text = element_text(size = 14), 
                                                     legend.title = element_text(size=14))
  if(guide == FALSE) {
    step1 + guides(fill = FALSE)
  }
  step1
}


################################## Precision 
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

transient.03 = read.csv("../results/transient.03.csv")

unique.variables = combined.models %>% filter(term.type == "single")  %>% 
  summarize(n = length(unique(variable)))


combined.models %>% 
  filter(term.type == "single") %>%
  filter(model.type == "full.set" & model != "full.model") -> freq.models


my.colors = colorRampPalette(brewer.pal(8, "YlOrRd"))(17)
my.colors = colorRampPalette(brewer.pal(8, "YlOrRd"))(max(transient.03$term, na.rm = TRUE))

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


transient.03 = read.csv("../results/transient.03.csv")
transient.03$time.1 = factor(transient.03$time.1,
                                 levels = c("freq.1", "freq.2", "freq.3", "gp.1", "gp.2", "accel"))
transient.03$variable1 = factor(transient.03$variable1, 
                               levels = rev(c("ratio.meanR", "ratio.mutation", "ratio.varR", "ratio.varSigma",
                                              "varR", "entropy", "infected", "day",  "individual.varSigma", 
                                              "antigenicDiversity",  "diversity" )))

my.colors = colorRampPalette(brewer.pal(8, "YlOrRd"))(18)
transient.03 %>%
  mutate(term  = row_number()) %>%
  filter(term.type == "single") %>%
  ggplot(aes(x = time.1, y = variable1, fill = as.factor(term))) + 
  geom_tile() + scale_fill_manual(values = rev(my.colors)) + 
  labs(x = "Time", y = "", fill = "Term", title = "Single Terms") +
  geom_text(aes(label = term))  + 
  scale_x_discrete(labels = c("Freq 1%", "Freq 2%", "Freq 3%", "EGP 1", "EGP 2", "Diff EGP")) -> transient.03.single.plot
transient.03.single.plot  
transient.03 %>%
  mutate(term  = row_number()) %>%
  filter(term.type == "int") %>%
  gather(key = variable.term, value = metric, variable1, variable2) %>%
  gather(key = time.description, value = descriptor, time.1, time.2) %>%
  arrange(term) 
  ggplot(aes(x = descriptor, y = metric)) + facet_wrap(~as.factor(term)) + geom_tile()




transient.01 = read.csv("~/Dropbox/current_fluantigen/model_csvs/full.model.trans01.csv")
transient.01 =transient.01[1:11,]
transient.01 %>%
  filter(term.type)


###############################################################
roc.03 = trans.03
roc.01 = read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.trans.01.csv")
roc.01 = roc.01[,-1]

roc.compare = rbind(data.frame(trans = .01, roc.01),
                    data.frame(trans = .03, roc.03))
roc.compare %>% 
  filter(thres > .49 & thres < .51) %>%
  group_by(trans) %>%
  mutate(fifty_thres = abs(.5-thres)) %>%
  slice(which.min(fifty_thres)) -> half.thres

roc.compare %>% 
  ggplot(aes(x = spec, y = sen, group = trans, color = as.factor(trans))) + geom_line(size = 1.5) +
  scale_x_reverse(breaks = seq(0,1,.1)) + labs(x = "Specificity", y  = "Sensitvity", color = "Transient Definition") +
  scale_color_manual(values = c("black", "darkgrey"), labels = c("1%", "3%"))  +
  scale_y_continuous(breaks = seq(0,1,.1)) -> roc.lines

roc.lines + geom_hline(data = half.thres, aes(yintercept = sen, color = as.factor(trans))) +
  geom_vline(data = half.thres, aes(xintercept = spec, color = as.factor(trans))) -> roc.comparison.plot



######################################## 
# Variable Zoom

full.model %>%
  arrange(term) 

freq.df.subset %>%
  filter(freq == "freq.02") -> freq.2.subset


freq.1.subset %>%
  select(success, name,  antigentype, ratio.varSigma, infected, diversity, antigenicDiversity) %>%
#  mutate_if(is.numeric, scale) %>%
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq1") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.1.variables
freq.2.subset %>%
  select(success, name, antigentype, ratio.varR, individual.varSigma) %>%
#  mutate_if(is.numeric, scale) %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq2") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.2.variables
freq.3.subset %>%
  select(success, name, antigentype, ratio.meanR, varR, ratio.mutation) %>%
#  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq3") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.3.variables
growth.1.subset %>%
  select(success, name, antigentype, ratio.meanR, ratio.mutation) %>%
#  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "gp.1") %>%
  mutate(id = paste0(variable, "_", time.point)) -> gp.1.variables
growth.2.subset %>%
  select(success, name, antigentype,ratio.meanR, entropy, day, ratio.varR, varR) %>%
#  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "gp.2") %>%
  mutate(id = paste0(variable, "_", time.point)) -> gp.2.variables
accelerate %>%
  select(success,name, antigentype, ratio.mutation, day) %>%
#  mutate_if(is.numeric, scale)   %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "accel") %>%
  mutate(id = paste0(variable, "_", time.point)) -> accel.variables


all.variables = rbind(freq.1.variables, freq.2.variables, freq.3.variables,
                      gp.1.variables, gp.2.variables, accel.variables)
unique(all.variables$id)
all.variables$id = factor(all.variables$id, levels = c("ratio.meanR_freq3", "ratio.meanR_gp.2", "varR_freq3", "ratio.mutation_freq3",
                                                       "entropy_gp.2", "ratio.varR_freq2", "ratio.meanR_gp.1", "ratio.varSigma_freq1",
                                                       "infected_freq1", "day_gp.2", "ratio.mutation_gp.1", "ratio.mutation_accel", 
                                                       "ratio.varR_gp.2", "individual.varSigma_freq2", "day_accel", "diversity_freq1", 
                                                       "antigenicDiversity_freq1", "varR_gp.2"))
                          labels = c("ratio.meanR: 3%", "ratio.meanR: EGP 2", "varR: 3%", "ratio.mutation: 3%",
                                     "entropy: EGP 2", "ratio.varR: 2%", "ratio.meanR: EGP 1", "ratio.varSigma: 1%", "infected: 1%", 
                                     "day: EGP 2", "ratio.mutation: EGP 1", "ratio.mutation: Diff EGP", "ratio.varR: EGP 2", "individual.varSigma: 2%",
                                     "day: Diff EGP", "diversity: 1%", "antigenicDiversity: 1%", "varR: EGP 2")
all.variables %>% 
  ggplot(aes(id, value, color = success)) + geom_boxplot() +
  scale_y_log10() + 
  scale_color_manual(values = c("black", "gray"), labels = c("Persistence", "Transient")) + 
  labs(x = "", y = "Scaled/Centralized Value", color = "Antigen Cluster Fate") + 
  theme(axis.text.x=element_blank()) +
  theme(legend.position = "top") -> boxplot.variables


transient.03 %>%
  filter(term.type == "single") -> single.03

single.03 %>%
  mutate(id = paste0(variable1, "_", time.1)) -> single.03
single.03$id = factor(single.03$id, levels = c("ratio.meanR_freq.3", "ratio.meanR_gp.2",  "varR_freq.3", "ratio.mutation_freq.3",  "entropy_gp.2",
                                                 "ratio.varR_freq.2",  "ratio.meanR_gp.1" , "ratio.varSigma_freq.1", "infected_freq.1", "day_gp.2" ,
                                                 "ratio.mutation_gp.1",  "ratio.mutation_accel", "ratio.varR_gp.2",  "individual.varSigma_freq.2", 
                                                 "day_accel" , "diversity_freq.1" ,  "antigenicDiversity_freq.1",  "varR_gp.2"),
                          labels = c("ratio.meanR: 5%", "ratio.meanR: EGP 2", "varR: 3%", "ratio.mutation: 3%",
                                     "entropy: EGP 2", "ratio.varR: 2%", "ratio.meanR: EGP 1", "ratio.varSigma: 1%", "infected: 1%",
                                     "day: EGP 2", "ratio.mutation: EGP 1", "ratio.mutation: Diff EGP", "ratio.varR: EGP 2", "individual.varSigma: 2%",
                                     "day: Diff EGP", "diversity: 1%", "antigenicDiversity: 1%", "varR: EGP 2"))

single.03 %>%
  ggplot(aes(x = id, y = estimate)) + geom_point() +
  geom_errorbar(aes(ymin = estimate -std.error, ymax = estimate+std.error)) +
  geom_hline(yintercept = 0, color = "grey54") +
  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1))  + 
  labs(x = "Variable", y = "Log Odds ") -> log.odds.plot

full.model.summary.plot = plot_grid(boxplot.variables, log.odds.plot, ncol = 1)
save_plot(full.model.summary.plot, filename = "exploratory.figs/full.model.summary.pdf",
          base_height = 8, base_aspect_ratio = 1.5)


  
#full.model %>%
#  ggplot(aes(x = id, y =exp(estimate))) + geom_point() +
#  geom_errorbar(aes(ymin = exp(estimate -std.error), ymax = exp(estimate+std.error))) +
#  geom_hline(yintercept = 1, color = "grey") +
#  theme(axis.text.x = element_text(size = 10, angle = 60,  hjust=1))  + 
#  labs(x = "Variable", y = "Odds Ratio") -> odds.ratio.plot


######################## Interactions

summary(model.null)

sjp.int(model.null, type = "eff")


############################# 2D Plots 
variable.1 = "freq.3.ratio.mutation"
variable.2 = "gp.2.ratio.varR"
all.variables

plot.1 = two.d.plot(data.set = all.variables, variable.1 = "diversity_freq1", variable.2 = "varR_gp.2", num.bins = 15, digits = 4)
int.1 = clean.2D.plot(plot.1, variable.name.1="Diversity: 1%", variable.name.2 = "Ratio.varR: EGP 2")

plot.2 = two.d.plot(data.set = all.variables, variable.1 = "ratio.meanR_gp.1", variable.2 = "day_gp.2", num.bins = 15, digits = 4)
int.2 = clean.2D.plot(plot.2, variable.name.1="Ratio Mutation: EGP 1", variable.name.2 = "Day: EGP 2")

plot.3 = two.d.plot(data.set = all.variables, variable.1 = "day_gp.2", variable.2 = "ratio.meanR_freq3", num.bins = 20, digits = 2)
int.3 = clean.2D.plot(plot.3, variable.name.1="Day: EGP 2", variable.name.2 = "Ratio.meanR: 3%")

plot.4 = two.d.plot(data.set = all.variables, variable.1 = "ratio.mutation_freq3", variable.2 = "ratio.varR_freq2", num.bins = 15, digits = 3)
int.4 = clean.2D.plot(plot.4, variable.name.1 = "Ratio.mutation: 3%", variable.name.2 = "ratio.varR: 2%")

plot.5 = two.d.plot(data.set = all.variables, variable.1 = "ratio.mutation_freq3", variable.2 = "ratio.varR_gp.2", num.bins = 15, digits = 4)
int.5 = clean.2D.plot(plot.5, variable.name.1 = "Ratio.mutation: 3%", variable.name.2 = "ratio.varR: EGP 2")

int.plots = cowplot::plot_grid(int.1, int.2, int.3, int.4, int.5, labels = "AUTO")
cowplot::save_plot(int.plots, filename = "exploratory.figs/int.plots.pdf",
                   base_height = 8, base_aspect_ratio = 1.9)


variable.1 = 
variable.2 = 


top.variables = two.d.plot(data.set = all.variables, variable.1 = "ratio.meanR_freq3", variable.2 = "ratio.meanR_gp.2", num.bins = 15, digits = 3)
top.5.plot.1 = clean.2D.plot(top.variables, variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")), 
                             variable.name.2 = expression(paste("Chage in ", italic("R"), " Advantage from 2-3%")))


top.variables = two.d.plot(data.set = dataPurged, variable.1 = "ratio.meanR_freq3", variable.2 = "varR_freq3", num.bins = 15, digits = 4)
top.5.plot.2 = clean.2D.plot(top.variables, variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")),
                             variable.name.2 = expression(paste("Variance in ", italic("R")," : 3%")))


top.variables = two.d.plot(data.set = dataPurged, variable.1 =  "ratio.meanR_freq3", variable.2 = "ratio.mutation_freq3", num.bins = 15, digits = 3)
top.5.plot.3 = clean.2D.plot(top.variables, variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")),
                             variable.name.2 = "Relative Deleterious Mutation Load: 3%")

top.variables = two.d.plot(data.set = dataPurged, variable.1 = "ratio.meanR_freq3", variable.2 = "entropy_gp.2", num.bins = 15, digits = 3)
top.5.plot.4 = clean.2D.plot(top.variables, variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")),
                             variable.name.2 = "Change in Shannon's Diversity Index from 2-3%")

top.5.variables = cowplot::plot_grid(top.5.plot.1, top.5.plot.2, top.5.plot.3, top.5.plot.4)

cowplot::save_plot(top.5.variables, filename = "exploratory.figs/top.5.variables.pdf", base_height = 8,
                   base_aspect_ratio = 1.3)



sample.plot = two.d.plot(data.set = dataPurged, variable.1 = "varR_freq3", variable.2 = "infected_freq1", num.bins = 15, digits = 5)
prob.plot = clean.2D.plot(sample.plot, variable.name.1 = "varR: 3%", variable.name.2 = "Infected: 3%")

variable.compare = c("varR_freq3","infected_freq1")

all.variables %>%
  filter(id %in% variable.compare) %>%
  select(-time.point,-variable) %>%
  spread(key = id, value = value) %>%
  ggplot(aes(x = varR_freq3, y = infected_freq1, color = success)) + geom_point() +
  scale_color_manual(values = c("black", "grey"))  +
  labs(x = "VarR: 3%", y  = "Infected: 1%") -> boxplot.compare
 
plot_grid(prob.plot, boxplot.compare)


############### Difference between ratio.mutation 
growth.2.subset %>% 
  select(success, ratio.varR) %>% 
  #mutate(direction =  ifelse(ratio.varR > 0, "positive", "negative")) %>%
  ggplot(aes(success, y = ratio.varR, color = success)) + geom_boxplot() + 
  scale_color_manual(values = c("black", "grey54")) + 
  labs(x = "Direction", y = "Ratio VarR: EGP 2", color = "Antigen Fate") +
  scale_x_discrete(labels = c("Decrease", "Increase")) -> ratio.varR.combined

save_plot(filename = "exploratory.figs/entropy.box.plot.pdf", plot = entropy.box.plot)


###################  viral Fitness 
viralType.frequencies = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.typeViralFitness.txt"), header = TRUE)) 
names(antigen.frequencies) = trial.dirs

meanViralFitness = map(trial.dirs, function(x) read.table(paste0(tropics.folder, x, "/out.viralFitnessSeries.txt"), header = TRUE))
names(meanViralFitness) = trial.dirs

antigen.frequencies[[3]] %>%
  distinct(day, antigentype, .keep_all = TRUE) %>%
  spread(key = antigentype, value = frequency, fill = 0) %>%
  gather(key = antigentype, value = frequency, -day, - infected) %>%
  mutate_at("antigentype", as.factor) %>%
  mutate(year = day/365 -5) %>%
  mutate(prevalence = infected*frequency *.0025) %>%
  filter(prevalence > 0) %>%
  ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
  geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
  labs(y = "Infected per 100K", x = "Year")  +
  scale_x_continuous(breaks = seq(from = 1, to = 25, by =2)) +
  guides(col = FALSE) + guides(fill = FALSE) -> prev.plot


meanViralFitness[[3]] %>%
  ggplot(aes(x = (day/365)-5, y = meanR)) + geom_line() + 
  labs(x = "Year") + 
  scale_x_continuous(breaks = seq(1,25, 2)) -> meanRplot


meanRsummary.plot = plot_grid(prev.plot, meanRplot, ncol = 1)


freq.df.subset %>%
  filter(name == "tropics_13") %>%
  group_by(antigentype) %>%
  summarize(success = success[1]) -> success.summary

meanViral.pop = meanViralFitness[[3]] %>% select(day, meanR)
colnames(meanViral.pop)[2] = "pop.meanR"


viralType.frequencies[[3]] %>%
  left_join(meanViral.pop, by = "day") %>%
  mutate_at("antigenType", as.character) %>%
  left_join(success.summary, by = c("antigenType" = "antigentype")) %>%
  filter(!is.na(success)) %>%
  ggplot(aes(x = (day/365)-5, y = exp(meanR), group = antigenType, color = antigenType)) +  geom_line(size = 1, alpha = .8) +
  #scale_color_manual(values = c("black", "grey")) +
  labs(y = "meanR", x = "Year") +
  scale_x_continuous(breaks = seq(1,25,2)) + facet_wrap(~success) + guides(color = FALSE) + 
  geom_line(aes(x=(day/365)-5, y = pop.meanR), color = "black") -> individual.meanR.plot

individual.meanR.plot

viralType.frequencies[[3]] %>%
  mutate_at("antigenType", as.character) %>%
  left_join(success.summary, by = c("antigenType" = "antigentype")) %>%
  filter(!is.na(success)) %>%
  group_by(antigenType) %>%
  mutate(days.alive = day-day[1]) %>%
  ggplot(aes(x = days.alive, y = exp(meanR), group = antigenType, color = success)) +  geom_line(size = 1) +
  scale_color_manual(values = c("black", "grey")) + labs(y = "meanR", x = "Days") +
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,50)) +
  scale_y_continuous(limits = c(0.9, 1.2)) -> meanR.days.alive.plot


viralType.frequencies[[3]] %>%
  mutate_at("antigenType", as.character) %>%
  left_join(success.summary, by = c("antigenType" = "antigentype")) %>%
  left_join(meanViral.pop, by = "day") %>%
  filter(!is.na(success)) %>%
  group_by(antigenType) %>%
  mutate(days.alive = day-day[1]) %>%  
  ggplot(aes(x = days.alive, y = exp(meanR)/pop.meanR, group = antigenType, color = success)) +  geom_line(size = 1) +
  scale_color_manual(values = c("black", "grey")) + labs(y = "ratio.R", x = "Days") +
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,50)) +
  scale_y_continuous(limits = c(0.9, 1.2)) -> ratio.meanR.days.alive.plot





viralType.frequencies[[3]] %>%
  left_join(meanViral.pop, by = "day") %>%
  mutate_at("antigenType", as.character) %>%
  left_join(success.summary, by = c("antigenType" = "antigentype")) %>%
  filter(!is.na(success)) %>%
  mutate(ratio.meanR = exp(meanR)/pop.meanR) %>%
  ggplot(aes(x = (day/365)-5, y = ratio.meanR, group = antigenType, color = antigenType)) +  geom_line(size = 1, alpha = .8) +
  #scale_color_manual(values = c("black", "grey")) +
  labs(y = "ratio.meanR", x = "Year") + facet_wrap(~success) + 
  scale_x_continuous(breaks = seq(1,25,2)) + guides(color = FALSE) + 
  geom_hline(yintercept = 1, color = "black")-> ratio.mean.plot


plot_grid(meanRsummary.plot, individual.meanR.plot, ratio.meanR.days.alive.plot,  ratio.mean.plot)


growth.1.subset %>%
  select(success, ratio.meanR) %>% 
  mutate(direction = ifelse(ratio.meanR > 0 , "positive", "negative")) %>%
  ggplot(aes(x = direction, y = ratio.meanR, color = success)) + geom_boxplot() + 
  scale_color_manual(values = c("black", "grey54"), labels = c("Decrease", "Increase")) +
  labs(x = "Direction", y = "Change in Ratio.meanR") -> boxplot.meanR



freq.df.subset %>% filter(freq == "freq.03") -> freq.3.subset

focus.data = as.data.frame(cbind(success = freq.3.subset$success, ratio.meanR = freq.3.subset$ratio.meanR, egp2 = growth.2.subset$ratio.meanR))

freq.3.subset %>%
  select(success, antigentype, name, ratio.meanR, ratio.mutation, varR) -> freq.3.variables

growth.2.subset %>%
  select(antigentype, name, ratio.meanR, entropy) -> gp2.variables

colnames(gp2.variables) = c("antigentype", "name", "gp2.ratio.meanR", "gp2.entropy")
focus.data = left_join(freq.3.variables, gp2.variables) %>% mutate_at("antigentype", as.factor)

focus.data %>%
#  filter(ratio.meanR) %>% 
  mutate(quadrant  = ifelse(gp2.entropy < 0 & ratio.meanR < 1, "lowerleft", 
                            ifelse(gp2.entropy > 0 & ratio.meanR < 1, "upperleft",
                                   ifelse(gp2.entropy> 0 & ratio.meanR > 1, "upperright", "lowerright")))) %>%
  filter(ratio.meanR > 1) -> focus.sub

focus.sub %>%
  group_by(quadrant) %>%
  summarize(num.rows = n()) %>%
  mutate(rel.freq = num.rows/sum(num.rows)) 


table.sucess = table(focus.sub$success, focus.sub$quadrant)

focus.sub %>% 
#ggplot(aes(x = success, y = gp2.ratio.meanR)) + geom_boxplot()  
  ggplot(aes(x = ratio.meanR, y = gp2.entropy, color = success)) + geom_point(alpha = .5) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 1)


variable.1 = "ratio.meanR"
variable.2 = "varR"

#variables = c(variable.1, variable.2)

var.1.breaks = seq(from = .95, to = 1.16, .02)
var.2.breaks = seq(from = 0, to =  .004, .0005)

  

focus.data %>%
    mutate(variable.1.group = cut(focus.sub[,`variable.1`], breaks = var.1.breaks),
           variable.2.group = cut(focus.sub[, `variable.2`], breaks = var.2.breaks)) %>%
    na.omit() %>%
    group_by(variable.1.group, variable.2.group, success) %>%
    summarize(num = n()) %>%
    mutate(rel.freq = num/sum(num)) %>%
    ungroup() %>%
    select(-num) %>%
    group_by(variable.1.group, variable.2.group) %>%
    spread(key = success, value = rel.freq, fill = 0) %>%
    gather(key = success, value = rel.freq, Transient, Est.) %>%
    filter(success == "Est.") %>%
    ggplot(aes(variable.1.group, variable.2.group, fill = rel.freq)) + geom_tile() +
    scale_x_discrete(labels = var.1.breaks) + scale_y_discrete(labels = var.2.breaks) +
    scale_fill_viridis(option="plasma")
}
clean.2D.plot <- function(plot, variable.name.1, variable.name.2, guide = TRUE) { 
  step1 = plot + labs(x = variable.name.1, y  = variable.name.2, fill = "Proportion of \nPersistent Antigens") +
    theme(axis.text = element_text(size =8), axis.text.x = element_text(angle = 60,  hjust=1),
          axis.title = element_text(size = 10)) + 
    theme(legend.key.size = unit(0.5, "cm")) +
    theme(legend.text = element_text(size = 8), legend.title = element_text(size=8))  
  if(guide == FALSE) {
    step1 + guides(fill = FALSE)
  }
  step1
}

