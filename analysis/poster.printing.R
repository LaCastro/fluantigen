library(hexbin)
library(viridis)

#scaled data 
input.data = read.csv("~/Dropbox/current_fluantigen/model_csvs/lower.vif.prop.csv")

############# Poster Scripts 
freq.1.subset %>%
  select(success, name,  antigentype, varR, ratio.varSigma, meanR, individual.varSigma) %>%
  #  mutate_if(is.numeric, scale) %>%
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq1") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.1.variables
freq.2.subset %>%
  select(success, name, antigentype, tmrca, antigenicTypes) %>%
  #  mutate_if(is.numeric, scale) %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq2") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.2.variables
freq.3.subset %>%
  select(success, name, antigentype, ratio.meanR, varR, varBeta, ratio.mutation, meanR) %>%
  #  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "freq3") %>%
  mutate(id = paste0(variable, "_", time.point)) -> freq.3.variables
growth.1.subset %>%
  select(success, name, antigentype, meanSigma, individual.meanSigma) %>%
  #  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "gp.1") %>%
  mutate(id = paste0(variable, "_", time.point)) -> gp.1.variables
growth.2.subset %>%
  select(success, name, antigentype, entropy, day, individual.meanSigma, meanSigma) %>%
  #  mutate_if(is.numeric, scale)  %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "gp.2") %>%
  mutate(id = paste0(variable, "_", time.point)) -> gp.2.variables
accelerate %>%
  select(success,name, antigentype, day) %>%
  #  mutate_if(is.numeric, scale)   %>% 
  gather(key = variable, value = value, -success, -antigentype, -name) %>%
  mutate(time.point = "accel") %>%
  mutate(id = paste0(variable, "_", time.point)) -> accel.variables


all.variables = rbind(freq.1.variables, freq.2.variables, freq.3.variables,
                      gp.1.variables, gp.2.variables, accel.variables)


#all.variables$id = factor(all.variables$id, levels = c("ratio.meanR_freq3", "ratio.meanR_gp.2", "varR_freq3", "ratio.mutation_freq3",
#                                                       "entropy_gp.2", "ratio.varR_freq2", "ratio.meanR_gp.1", "ratio.varSigma_freq1",
#                                                       "infected_freq1", "day_gp.2", "ratio.mutation_gp.1", "ratio.mutation_accel", 
#                                                       "ratio.varR_gp.2", "individual.varSigma_freq2", "day_accel", "diversity_freq1", 
#                                                       "antigenicDiversity_freq1", "varR_gp.2"))

all.variables$id = factor(all.variables$id, levels = c("ratio.meanR_freq3", "varR_freq3", "individual.meanSigma_gp.2",
                                                       "meanSigma_gp.2", "individual.varSigma_freq1", "entropy_gp.2", "individual.meanSigma_gp.1",
                                                       "antigenicTypes_freq2", "ratio.varBeta_freq3", "ratio.mutation_freq3", "day_gp.2",
                                                       "meanSigma_gp.1", "ratio.varSigma_freq1", "meanR_freq1", "tmrca_freq1", "varR_freq1", "day_accel"))


########## Figure 4 A 
fill = "rel.freq"

cut_data <- function(data, variables, type, bins) { 
  all.variables %>%
    filter(id %in%variables) %>%
    select(success,name,antigentype,id,value) %>%
    spread(key = id, value = value) -> data.set
  
  data.set %>%
    select(success, one_of(variables)) %>%
    mutate(variable.1.group = cut(data.set[,`variable.1`], breaks = bins),
           variable.2.group = cut(data.set[, `variable.2`], breaks = bins)) %>%
    na.omit() -> cut.data
  
  cut.data %>%   
    group_by(variable.1.group, variable.2.group, success) %>%
    summarize(num = n()) %>%
    mutate(rel.freq = num/sum(num)) -> rel.freq.data
  if (type == "Est.") {
    cut.data %>% 
      left_join(rel.freq.data) %>% 
      ungroup() %>%
      select(-num) %>%
      filter(success == "Est.") -> return.data 
  } else { 
    cut.data %>% 
      left_join(rel.freq.data) %>%
      ungroup() %>%
      select(-num) %>%
      filter(success == "Transient" & rel.freq == "1") -> return.data
  } 
  return(return.data)
}

variable.1 = "ratio.meanR_freq3"; variable.2 = "varR_freq3"
variables = c(variable.1, variable.2)

variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")) 
variable.name.2 = expression(paste("Variance in Population ", italic("R")))

est.data = cut_data(all.variables, variables, type = "Est.", bins = 15)
transient.data = cut_data(all.variables, variables, type = "Transient", bins = 15)
                                   
est.data %>%
  ggplot(aes(x = ratio.meanR_freq3, y = varR_freq3, z = rel.freq)) + stat_summary_hex(fun = mean, bins = 40)  + 
  scale_fill_viridis(option="plasma") + labs(x = variable.name.1, y  = variable.name.2, fill = "Prop. \n Successful") +
  geom_point(data = transient.data, aes(x = ratio.meanR_freq3, y = varR_freq3), color = "grey", alpha = .5) + 
  theme(legend.position = "none") -> one.vs.two


################### Figure 4 B #################################
variable.1 = "ratio.meanR_freq3"; variable.2 = "individual.meanSigma_gp.2" 
variables = c(variable.1, variable.2)

variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")) 
variable.name.2 = expression(paste("Susceptibility to Antigen: 2-3%"))

est.data = cut_data(all.variables, variables, type = "Est.", bins = 15)
transient.data = cut_data(all.variables, variables, type = "Transient", bins = 15)

est.data %>%
  ggplot(aes_string(x = variable.1, y = variable.2, z = fill)) + stat_summary_hex(fun = mean, bins = 40)  + 
  scale_fill_viridis(option="plasma") + labs(x = variable.name.1, y  = variable.name.2, fill = "Prop. \n Successful") +
  geom_point(data = transient.data, aes_string(x = variable.1, y = variable.2), color = "grey", alpha = .5) + 
  theme(legend.position = "none") -> one.vs.three

##################### Figure 4 C ################################
variable.1 = "ratio.meanR_freq3"; variable.2 = "meanSigma_gp.2"
variables = c(variable.1, variable.2)

variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")) 
variable.name.2 = expression(paste("Population Susceptibility : 2%"))

est.data = cut_data(all.variables, variables, type = "Est.", bins = 20)
transient.data = cut_data(all.variables, variables, type = "Transient", bins = 20)

est.data %>%
  ggplot(aes_string(x = variable.1, y = variable.2, z = fill)) + stat_summary_hex(fun = mean, bins = 50)  + 
  scale_fill_viridis(option="plasma") + labs(x = variable.name.1, y  = variable.name.2, fill = "Prop. \n Successful") +
  geom_point(data = transient.data, aes_string(x = variable.1, y = variable.2), color = "grey", alpha = .5) + 
  theme(legend.position = "none") -> one.vs.four

##################### Figure 4 D ################################

variable.1 = "ratio.meanR_freq3"; variable.2 = "individual.varSigma_freq1"
variables = c(variable.1, variable.2)

variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")) 
variable.name.2 = expression(paste("Susceptibility to Angtigen (Var) : 1%"))

est.data = cut_data(all.variables, variables, type = "Est.", bins = 20)
transient.data = cut_data(all.variables, variables, type = "Transient", bins = 20)

est.data %>%
  ggplot(aes_string(x = variable.1, y = variable.2, z = fill)) + stat_summary_hex(fun = mean, bins = 40)  + 
  scale_fill_viridis(option="plasma") + labs(x = variable.name.1, y  = variable.name.2, fill = "Prop. \n Successful") +
  geom_point(data = transient.data, aes_string(x = variable.1, y = variable.2), color = "grey", alpha = .5) + 
  theme(legend.position = "none") -> one.vs.five


##################### Figure 4 E  ################################
variable.1 = "ratio.meanR_freq3"; variable.2 = "entropy_gp.2"
variables = c(variable.1, variable.2)

variable.name.1 = expression(paste("Relative ", italic("R"), " Advantage: 3%")) 
variable.name.2 = expression(paste("Change in ", italic("H"),  ":2-3%"))

est.data = cut_data(all.variables, variables, type = "Est.", bins = 20)
transient.data = cut_data(all.variables, variables, type = "Transient", bins = 20)

est.data %>%
  ggplot(aes_string(x = variable.1, y = variable.2, z = fill)) + stat_summary_hex(fun = mean, bins = 40)  + 
  scale_fill_viridis(option="plasma") + labs(x = variable.name.1, y  = variable.name.2, fill = "Prop. \n Successful") +
  geom_point(data = transient.data, aes_string(x = variable.1, y = variable.2), color = "grey", alpha = .5) -> one.vs.six
one.vs.six.b = one.vs.six + theme(legend.position = "none")

legend <- get_legend(one.vs.six)
plots <- align_plots(one.vs.two, one.vs.three, one.vs.four, one.vs.five, one.vs.six.b, align = "v", axis = "l")

bottom.row = plot_grid(plots[[4]], plots[[5]], legend, nrow = 1)
first.row = plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow = 1)

final.plot = plot_grid(first.row, bottom.row, ncol = 1)
save_plot(final.plot, filename = "exploratory.figs/modulation.variables.pdf", base_height = 8, base_aspect_ratio = 2)

################################ Figure 1 -- Full Overview 
timeseries = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.timeSeries.txt"), header = TRUE))
timeseries[[5]] %>%
    ggplot(aes(x = date, y = totalI*.0025)) + geom_line() +
    labs(x = "Year", y = "Infected") +
    scale_x_continuous(breaks = seq(1, 25,2)) -> timeseries.plot

freq.1.subset %>% filter(name == "tropics_18") %>% select(antigentype) -> successful.antigens
num.transitions = nrow(successful.antigens)
myColors = colorRampPalette(brewer.pal(8, "Accent"))(num.transitions)

## Antigen Dynamics   
antigen.frequencies[[5]] %>%
    filter(antigentype %in% successful.antigens$antigentype) %>%
    distinct(day, antigentype, .keep_all = TRUE) %>%
    mutate_at("antigentype", as.factor) %>%
    spread(key = antigentype, value = frequency, fill = 0) %>%
    gather(key = antigentype, value = frequency, -1, -2, - infected)  -> antigen.freq.long
  
antigen.freq.long %>%
    mutate(year = day/365 - 5) %>%
    mutate(prevalence = infected*frequency *.0025) %>% #*.0025
    filter(prevalence > 0) %>%
    ggplot(aes(x = year, y = prevalence, fill = antigentype)) +
    geom_area(color = "black", aes(color = antigentype, fill = antigentype)) +
    scale_color_manual(values = myColors) + scale_fill_manual(values = myColors)+
    labs(y = "Infected per 100K ", x = "Year") + scale_x_continuous(breaks = seq(from = 1, to = 25, by =2)) +
    guides(col = FALSE) + guides(fill = FALSE) + scale_y_continuous(breaks = seq(from = 0, 150, by = 25)) +
    theme(axis.text = element_text(size = 18), axis.title = element_text(size  = 18)) -> antigen.dynamics 

### Figuring out the time slot
antigen.time = bind_rows(freq.1.subset, freq.2.subset,freq.3.subset)  %>%  
  filter(antigentype == "5316") %>%
  select(simDay, freq) %>%
  mutate(day = simDay*365)

antigen.freq.long%>%
  mutate(status = ifelse(antigentype == "5316", "focal", "ignore")) %>%
  ggplot(aes(x = (day/365)-5, y = frequency, group = antigentype, fill = status)) + geom_area(color = "grey54") +
  scale_fill_manual(values = c("black", "grey")) + labs(x = "Year", y = "Relative Frequency") +
  guides(fill = FALSE) +  scale_x_continuous(breaks = seq(from= 7.5, to =11.5, by=.5), limits = c(9, 11.5)) -> relative.frequency.plot

relative.frequency.plot + geom_vline(data = antigen.time, aes(xintercept = simDay, color = freq), size = 2) +
  scale_color_brewer(palette = "Paired") +
  guides(color = FALSE) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18)) +
  labs(color = "Threshold") -> plot.RF

plot.RF + annotate("text", x = 9.57, y = .50, label = "1%", size = 8) +
  annotate("text", x = 9.67, y = .6, label = "2%", size = 8) +
  annotate("text", x = 9.9, y = .7, label = "3%", size = 8) -> relative.frequency.plot

# Plot the variable set 
viral.fitness.metrics = read.table("../data/tropics/eligible/tropics_18/out.viralFitnessSeries.txt", header = TRUE)

viral.fitness.metrics %>%
    select(day, meanR) %>% 
    gather(key = variable, value = value, -day) -> viral.subset
viral.subset$variable = factor(viral.subset$variable, labels = c("Mean Reproductive Number"))

viral.subset %>% 
  ggplot(aes(x=day/365-5, y = value)) + 
    geom_line(size  =2 ) +facet_wrap(~variable, scales = "free") +
    labs(x = "Year", y = "Population Mean R") +
    scale_x_continuous(breaks = seq(from= 7.5, to =12.5, by=.5), limits = c(9 ,11.5)) +
    theme(text = element_text(size = 18)) +   
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18)) +
    theme(strip.text = element_blank())-> viral.metrics.plot 

viral.metrics.plot + geom_vline(data = antigen.time, aes(xintercept = simDay, color = freq), size = 2) + 
  scale_color_brewer(palette = "Paired", labels = c("1%", "2%", "3%")) + labs(color = "Frequency Time Point") +
  guides(color = FALSE) -> viral.metrics.plot


row1 = plot_grid(antigen.dynamics, labels = "AUTO")
row2 = plot_grid(relative.frequency.plot, viral.metrics.plot, labels = c("B", "C"))
poster.summary1 = plot_grid(row1, row2,  nrow = 2) 
save_plot(filename = "exploratory.figs/poster.summary1.pdf", plot = poster.summary1, base_height = 9, base_aspect_ratio = 1.8)
  

############################ Rank Figure 
########## Model Peformance
################################## Comparing Model Performances
full.model = read.csv("~/Dropbox/current_fluantigen/model_csvs/transient.03.csv")
pop.model = read.csv("~/Dropbox/current_fluantigen/model_csvs/trans03.poplevel.csv")
pop.3.model = read.csv("~/Dropbox/current_fluantigen/model_csvs/trans03.popleve3only.csv")

full.model %>%
  select(term.type, variable, time.1, per) %>%
  filter(term.type == "single") %>%
  mutate(model = "full") -> full.model.sub
pop.model %>%
  select(term.type, variable, time.1, per) %>%
  filter(term.type == "single") %>%
  mutate(model = "pop") -> pop.model.sub
pop.3.model %>%
  select(term.type, variable, time.1, per) %>%
  mutate(model = "pop.3") -> pop.3.model

combined.models = rbind(full.model.sub, pop.model.sub, pop.3.model)
combined.models %>%
  group_by(model) %>%
  mutate(rank = row_number()) -> combined.models

combined.models %>%
  ggplot(aes(rank, per, color = model, group = model)) + geom_line(size = 2) +
  scale_color_manual(values = c("purple", "orange", "grey"),
                     labels = c("Complete", "Real World", "Real World: 3%")) +
  scale_x_continuous(breaks = seq(1,18,2)) + labs(x = "Model Term", y = "Average AUC", color = "") +
  scale_y_continuous(breaks = seq(0.55, .95,.05)) +
  theme(legend.position = "left", legend.title = element_text(size = 14)) + 
  theme(axis.text = element_text(size = 14), text = element_text(size  = 14)) -> model.performance


roc.full = read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.trans.03.csv")
roc.03.pop3 = read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.t03.poplevel.csv") 
roc.pop = read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.transpop.03.csv") 

combined.roc = rbind(data.frame(model = "full", roc.full),
                     data.frame(model = "pop.3", roc.03.pop3),
                     data.frame(model = "pop", roc.pop))

combined.roc %>%
  ggplot(aes(x = 1-spec, y = sen, color = model, group = model)) + geom_line(size = 2) +
  labs(x = "False Positive rate \n (1-Specificity)", y = "True Positive rate \n (Sensitivity)", color = "Variable Set") +
  scale_color_manual(values = c("purple", "orange", "grey"),
                     labels = c("Complete", "Real World", "Real World: 3%")) +
  guides(color = FALSE) + 
  theme(axis.text = element_text(size = 14), text = element_text(size  = 14)) +
  scale_y_continuous(breaks = seq(0,1,.2)) + scale_x_continuous(breaks = seq(0,1,.2)) -> model.roc

#########
results.03 = read.csv("~/Dropbox/current_fluantigen/model_csvs/transient.03.csv")
results.03.pop = read.csv("~/Dropbox/current_fluantigen/model_csvs/trans03.poplevel.csv")
results.03.freq.3 = read.csv("~/Dropbox/current_fluantigen/model_csvs/trans03.popleve3only.csv")


results.03 %>%
  select(-variable2, -time.2) %>%
  filter(term.type == "single") -> results.03
results.03.pop %>%
  filter(term.type == "single") -> results.03.pop
results.03.freq.3 %>%
  filter(term.type == "single") -> results.03.freq.3

compare.results = rbind(data.frame(model = "pop", results.03.pop),
                         data.frame(model = "full", results.03),
                         data.frame(model = "pop.3", results.03.freq.3))
compare.results %>%
  group_by(model) %>%
  mutate(rank = row_number()) -> compare.results.rank
my.colors = colorRampPalette(brewer.pal(8, "YlOrRd"))(max(compare.results.rank$rank))
compare.results.rank$variable = factor(compare.results.rank$variable, 
                                       levels = rev(c("ratio.meanR", "varR", "ratio.mutation", "entropy", "ratio.varR", 
                                                      "ratio.varSigma", "infected", "day", "meanR", "individual.varSigma", "diversity", 
                                                      "antigenicDiversity")),
                                       labels = rev(c("Rel. R Advantage", "Pop. R (Var.)", "Rel.Deleterious Mutation Load", "Shannon's Diversity", 
                                                  "Rel. R Advantage (Var.)", "Rel. Susceptibility (Var.)", "Prevalence", "Days", "Pop. R",
                                                  "Focal Antigen Susceptibility (Var.)", "Diversity", "Antigenic Diversity")))
compare.results.rank$time.1 = factor(compare.results.rank$time.1, levels = c("freq.1", "freq.2", "freq.3", "gp.1","gp.2", "accel"))
compare.results.rank$model = factor(compare.results.rank$model, levels = c("full", "pop", "pop.3" ),labels = c("Complete", "Real World", "Real World: 3%"))
compare.results.rank %>% 
  #filter(model == "full") %>%
  ggplot(aes(x = time.1, y = variable, fill = as.factor(rank))) + geom_tile() +
  scale_fill_manual(values = rev(my.colors)) + 
  facet_wrap(~model) + 
  labs(x = "Data Time Point", y = "", fill = "Term Order") +
  scale_x_discrete(labels = c("RF \n 1%", "RF \n 2%", "RF \n3%", "EGP \n1", "EGP \n2", "Diff \nEGP")) +
  geom_text(aes(label = rank)) + guides(fill = FALSE) + 
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14)) -> rank.compare.plot
#save_plot(rank.compare.plot, filename = "exploratory.figs/rank.compare.plot.pdf", base_height = 8)

second.row = plot_grid(model.performance, model.roc, rel_widths = c(1.2, 1), labels = "AUTO")
combined.rank.plot = plot_grid(second.row, rank.compare.plot, ncol = 1 , rel_heights = c(.7, 1), labels = c("", "C"))
save_plot(combined.rank.plot, filename = "exploratory.figs/combined.rank.plot.pdf", base_height = 9, base_aspect_ratio = 1.8)



################### Calculating Total Antigen life
trial = trial.dirs[[3]]
## Create a function that takes in a name, reads the antigen frequencies, gets success list, and calculates average life 
test.ant = antigen.frequencies[[3]]

freq.df.subset %>%
 filter()
calculate_days_alive = function(trial) {
  success.types = freq.1.subset %>%
    filter(name == trial) %>%
    filter(success == "Est.") %>%
    select(antigentype)
  
  ant.freq = antigen.frequencies[[eval(trial)]]
  
  ant.freq %>%
    group_by(antigentype) %>%
    filter(antigentype %in% success.types$antigentype) %>%
    summarize(days.alive = day[n()]-day[1])
}

days.alive = map(trial.dirs, calculate_days_alive)
names(days.alive) = trial.dirs
days.alive.df = do.call("rbind", days.alive)
days.alive.df %>%
  ggplot(aes(days.alive/365)) + geom_histogram()

days.alive.df %>%
  summarize(median.days = median(days.alive/365),
            lower.quantile = quantile(days.alive/365, probs = .25),
            upper.quantile = quantile(days.alive/365, probs = .75))


#########################################
data.set = read.csv("../results/")



variable.1 = "ratio.meanR_freq3"; variable.2 = "varR_freq3"; variable.3 = "entropy_gp.2"

variables = c(variable.1, variable.2, variable.3)

all.variables %>%
  filter(id %in%variables) %>%
  select(success,name,antigentype,id,value) %>%
  spread(key = id, value = value) -> data.set

head(data.set)

data.set %>% 
  ggplot(aes(ratio.meanR_freq3, varR_freq3,color = entropy_gp.2, size = entropy_gp.2)) + 
  geom_point(alpha = .8) + scale_color_viridis(option="plasma") + 
  facet_wrap(~success)
