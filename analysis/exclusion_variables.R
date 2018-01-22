######## Var R 
viral.fitness[[3]] %>%
  select(varR) -> varR.trial

track.antigen[[3]] %>%
  select(day, antigenicTypes, I) -> antigenicTypes.trial
trial.data = cbind(varR.trial, antigenicTypes.trial)
trial.data %>%
  gather(key = variable, value = value, -day) -> trial.data.l

trial.data %>%
  gather(key = variable, value = value, -day) %>%
  ggplot(aes(x = day, y = value)) + geom_line() + facet_wrap(~variable, ncol = 1, scales = "free") -> varR.by.infected


######## Does variance in frequency correspond to variance in frequency or hthe dominant type

# Caclulating dominant. freq for each time step 
antigen.frequencies[[3]] %>%
  group_by(day) %>%
  summarize(dom.freq = max(frequency),
            antigentype = antigentype[which.max(frequency)]) -> dominant.freq
dominant.freq$variable =  "dom.freq"
colnames(dominant.freq) = c("day", "value", "antigentype", "variable")


dominant.freq %>%
  select(-antigentype) %>%
  bind_rows(trial.data.l) 
  spread(key = variable, value = value) %>%
  ggplot(aes(x=dom.freq, y =antigenicTypes)) + geom_point() + geom_smooth(method = "lm")

cor.test(varR.dataset$dom.freq, varR.dataset$antigenicTypes)

#  ggplot(aes(x = day/365-5, y = value)) + geom_line() +facet_wrap(~variable, scales = "free", ncol = 1)  +
#  scale_x_continuous(breaks = seq(0,25, 1))




varR.population.plot = plot_grid(antigen.dynamics, varR.plot,  ncol = 1)


# Caclulating dominant. freq for each time step 
antigen.frequencies[[3]] %>%
  group_by(day) %>%
  summarize(mean.frequency = mean(frequency)) -> mean.frequency

antigen.frequencies[[3]] %>%
  left_join(mean.frequency) %>%
  mutate(diff = (frequency-mean.frequency)^2) %>%
  group_by(day) %>%
  summarize(sum = sum(diff),
            num = n()) %>%
  mutate(variance = sum/(num-1)) -> variance

variance %>%
  select(day, variance) %>%
  left_join(varR.dataset) -> varR.dataset

cor.test(varR.dataset$antigenicTypes, varR.dataset$variance)


varR.dataset %>%
  mutate_at(-c(1), scale) %>%
  select(day, varR, antigenicTypes) %>%
  gather(key = variable, value = value, -day) %>%
  ggplot(aes(x = day, y = value, group = variable, color = variable)) + geom_line(size = 1.2) +
  scale_x_continuous(breaks = seq(1,25,2)) +
  scale_color_manual(values = c("orange", "purple"))
  
  #ggplot(aes(x = varR, y = variance)) + geom_point() + geom_smooth(method = "lm")
    
        

two.d.plot(data.set = dataPurged, variable.1 = "diversity_freq1", variable.2 = "varR_gp.2", num.bins = 15, digits = 4)





#################### ratio.meanR change
freq.3.subset %>%
  filter(name == "tropics_13") %>%
  select(antigentype,success) -> trial.antigens

viral.type.fitness[[3]] %>%
  select(day, meanR, antigenType) %>%
  filter(antigenType %in% trial.antigens$antigentype) %>%
  mutate_at("antigenType", as.character) %>%
  left_join(trial.antigens, by = c("antigenType" = "antigentype")) -> viral.types.R


freq.2.subset %>%
  filter(name == "tropics_13") %>%
  select(antigentype, day) -> day.start; colnames(day.start)[2] = "day.start"
  
freq.3.subset %>%
  filter(name == "tropics_13") %>%
  select(antigentype, day) -> day.finish ; colnames(day.finish)[2] = "day.finish"

viral.types.R %>%
  left_join(day.start, by = c("antigenType" = "antigentype")) %>%
  left_join(day.finish, by = c("antigenType" = "antigentype")) %>%
  mutate_at("antigenType", as.factor) %>%
  group_by(antigenType) %>%
  filter(day >= day.start[1] & day <= day.finish[1]) %>%
  mutate(time.point = row_number(),
         diff = day.finish-day.start) %>%
  slice(c(1, n())) %>%
  summarize(mean.R.diff = exp(meanR[2])-exp(meanR[1]),
            success = success[1],
            diff = diff[1],
            finishing.meanR = exp(meanR[2])) %>%
  ggplot(aes(x = finishing.meanR, y = mean.R.diff, group = antigenType, color = success)) + geom_point()






viral.fitness[[3]] %>%
  select(day, varR) %>%
  rename(pop.varR = varR) -> pop.varR

################ Look at ratios of varR
### Identify successful types
freq.df.subset %>%
  filter(name == "tropics_13") %>%
  group_by(antigentype) %>%
  select(antigentype, success) -> successful.antigens


viral.type.fitness[[3]] %>%
  mutate_at("antigenType", as.character) %>%
  left_join(successful.antigens, by = c("antigenType" = "antigentype")) %>%
  filter(!is.na(success)) %>%
  left_join(pop.varR) %>%
  mutate(exp(varR))

time.series[[3]] %>%
  ggplot(aes(x = (date*365), y = antigenicDiversity)) + geom_line() -> antigenic.diversity.plot

entropy[[3]] %>%
  ggplot(aes(x = (day/365)-5, y = entropy)) + geom_line() -> entropy.plot

plot_grid(antigenic.diversity.plot, entropy.plot, ncol = 1)
  
 


















#########################################################################
## Poster Figures 

transient.full <- read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.trans.03.csv")
trans03.poplevel <- read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.t03.poplevel.csv")
trans03.poplevel03 <- read.csv("~/Dropbox/current_fluantigen/model_csvs/roc.transpop.03.csv")

transient.full$model = "full"
transient.full=transient.full[,-1]
trans03.poplevel$model = "pop"
trans03.poplevel03$model = "pop03"

trans03.poplevel = trans03.poplevel[,-1]
trans03.poplevel03=trans03.poplevel03[,-1]

trans.per = rbind(transient.full, trans03.poplevel, trans03.poplevel03)
trans.per %>%
  ggplot(aes(x = sen, y = spec, color = model, group = model)) + geom_line(size = 1.2) +
  labs(x = "True Positive Rate", y = "True Negative Rate", color = "Model Variables") + 
  scale_color_manual(values = c("orange", "purple", "grey54"), labels = c("Full Set", "Measurable Subset", "Measurable Subset: 3%")) +
  geom_abline(slope = -1, intercept = 1, alpha =.5) + 
  scale_y_continuous(breaks = seq(0,1,.2)) + scale_x_continuous(breaks = seq(0,1,.1)) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 11)) + 
  theme(axis.title = element_text(size = 10), axis.text = element_text(size =8)) -> comparative.roc


########## Model Peformance
results.03 = read.csv("../results/transient.03.csv")
results.03.pop = read.csv("../results/trans03.poplevel.csv")
results.03.freq.3 = read.csv("../results/trans03.popleve3onlyl.csv")
results.03 %>%
  select(-variable2, -time.2) -> results.03


compare.results  = rbind(data.frame(model = "pop", results.03.pop),
                         data.frame(model = "full", results.03),
                         data.frame(model = "pop.3", results.03.freq.3))

compare.results %>% 
  group_by(model) %>%
  filter(term.type == "single") %>%
  mutate(term.number = row_number()) %>%
  ggplot(aes(x = term.number, y = per, group = model, color = model)) + geom_line(size = 1.2) +
  labs(x = "Term", y = "Average AUC on 5-fold cross validation", color = "Model Variables") +
  scale_x_continuous(breaks = seq(1,18,2)) + guides(color = FALSE) + 
  scale_color_manual(values = c("orange", "purple", "grey54")) + 
  scale_y_continuous(breaks = seq(.5,.95,.05)) + 
  theme(axis.title = element_text(size = 10), axis.text = element_text(size =8)) -> comparative.permcurve

compare.summary.lines = plot_grid(comparative.permcurve, comparative.roc, rel_widths = c(.9, 1.3))

########### Rank plots 
compare.results %>%
  group_by(model) %>%
  filter(term.type == "single") %>% 
  mutate(rank = row_number()) -> compare.results.rank

my.colors = colorRampPalette(brewer.pal(8, "YlOrRd"))(max(compare.results.rank$rank))
compare.results.rank$variable = factor(compare.results.rank$variable, 
                              levels = rev(c("ratio.meanR", "varR", "ratio.mutation", "entropy", "ratio.varR", 
                                             "ratio.varSigma", "infected", "day", "meanR", "individual.varSigma", "diversity", 
                                             "antigenicDiversity", "accel.day*freq.1.infected", "gp.2day*freq.1.infected")))
compare.results.rank$time.1 = factor(compare.results.rank$time.1, levels = c("freq.1", "freq.2", "freq.3", "gp.1","gp.2", "accel"))
compare.results.rank$model = factor(compare.results.rank$model, levels = c("full", "pop", "pop.3" ),labels = c("Full Set", "Measurable Subset", "Measurable Subset: 3%"))

compare.results.rank %>% 
  #filter(model == "full") %>%
  ggplot(aes(x = time.1, y = variable, fill = as.factor(rank))) + geom_tile() +
  scale_fill_manual(values = rev(my.colors)) + 
   facet_wrap(~model) + 
  labs(x = "Data Time Point", y = "Variable", fill = "Term") +
  scale_x_discrete(labels = c("Freq 1%", "Freq 2%", "Freq 3%", "EGP 1", "EGP 2", "Diff EGP")) +
  geom_text(aes(label = rank)) + guides(fill = FALSE) + 
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) -> rank.compare.plot


save_plot(compare.summary.lines, filename = "~/Dropbox/current_fluantigen/model_plots/compare.summary.lines.pdf", base_aspect_ratio = 2.3)
save_plot(rank.compare.plot, filename = "~/Dropbox/current_fluantigen/model_plots/rank.compare.plot.pdf", base_height = 8,base_aspect_ratio = 1.5)
compare.summary.lines


compare.3.plot = plot_grid(compare.summary.lines, rank.compare.plot, ncol = 1, labels = "AUTO")
save_plot(plot = compare.3.plot, filename = "~/Dropbox/current_fluantigen/model_plots/compare.3.pdf", base_height = 8, base_aspect_ratio = 1.5)
