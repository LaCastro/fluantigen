## Exclusion 


# What I believe was taken out at this point 
# vif.excluded = c("freq.2.diversity", "freq.2.antigenicDiversity", "freq.2.entropy",
#                 "freq.3.diversity", "freq.3.antigenicDiversity", "freq.2.dominant.freq",
#                 "freq.2.infected", "freq.1.entropy", "freq2.meanR", "freq.2.meanR",
#                 "freq.2.totalS", "freq.3.dominant.freq", "freq.2.tmrca")


###### 1/22 ################################################### 
# Snapshot 
vif.excluded = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR", "individual.meanSigma",
                 "ratio.meanSigma", "meanSigma", "covBetaSigma", "meanBeta", "individual.varR", "individual.varBeta", "ratio.varR")

# Growth 
vif.excluded = c("individual.meanMut", "individual.varBeta", "individual.varR", "ratio.meanR", "meanR", "ratio.meanBeta",
                 "individual.meanR", "ratio.meanSigma", "varSigma", "individual.varMut")

# Combined 
vif.excluded = c("freq.2.ratio.meanBeta", "freq.1.ratio.meanBeta", "freq.2.meanR", "diff.1.meanSigma",
                 "freq.1.meanR", "diff.1.individual.meanR", "diff.1.individual.meanSigma", "diff.2.individual.meanSigma",
                 "freq.2.varSigma", "freq.2.ratio.meanR", "freq.3.varSigma", "diff.1.ratio.varR", "accel.ratio.mutation",
                 "accel.ratio.varR", "diff.1.individual.meanBeta", "accel.ratio.varBeta", "accel.individual.meanBeta",
                 "freq.1.ratio.mutation", "freq.2.ratio.mutation", "accel.meanLoad", "freq.3.ratio.mutation")



##### 1/24 ######################################################

## Exclusion based on collinearity after everything has been combined and up until VIF < 5 
vif.excluded = c("freq.1.individual.meanMut", "freq.1.individual.meanBeta", "freq.2.individual.meanMut", "freq.3.individual.meanMut",
                 "diff.1.individual.meanMut", "diff.1.individual.meanBeta", "diff.2.individual.meanMut", "accel.individual.meanMut",
                 "diff.1.individual.varMut", "diff.1.individual.varBeta", "diff.1.individual.varR", "accel.individual.varMut",
                 "accel.individual.varBeta", "accel.individual.varR", "freq.2.individual.varMut", "freq.3.individual.varMut", 
                 "freq.1.individual.varBeta", "freq.1.individual.varR", "diff.2.individual.varMut", "diff.1.ratio.meanR",
                 "accel.ratio.meanR", "freq.2.individual.meanBeta", "freq.3.individual.meanBeta", "diff.1.individual.meanR",
                 "diff.1.ratio.meanSigma", "accel.ratio.meanSigma", "freq.2.ratio.meanSigma", "freq.3.ratio.meanSigma",
                 "diff.1.meanR", "diff.2.individual.meanR", "freq.2.individual.meanR", "freq.1.individual.meanR", "accel.ratio.meanBeta",
                 "accel.meanR", "freq.2.meanR", "freq.3.meanR", "freq.1.individual.meanSigma", "freq.2.individual.meanSigma", 
                 "freq.3.individual.meanSigma", "diff.2.ratio.meanR", "freq.2.ratio.meanR", "accel.individual.meanBeta", "accel.meanSigma",
                 "freq.3.individual.meanR", "diff.2.ratio.meanBeta", "freq.1.ratio.meanBeta", "freq.3.ratio.meanBeta", "accel.individual.meanSigma",
                 "freq.1.ratio.meanR", "diff.2.meanR", "diff.2.ratio.meanSigma", "diff.1.varSigma", "freq.1.meanSigma", "freq.2.meanSigma",
                 "freq.3.meanSigma", "freq.1.ratio.meanSigma", "accel.varSigma", "freq.2.varSigma", "freq.3.varSigma", "accel.individual.meanR",
                 "accel.varR", "freq.1.covBetaSigma", "freq.2.covBetaSigma", "freq.3.covBetaSigma", "diff.2.varSigma", "freq.1.meanBeta",
                 "freq.2.meanBeta", "freq.3.meanBeta", "accel.ratio.varBeta", "freq.2.individual.varBeta", "freq.2.individual.varR",
                 "freq.3.individual.varBeta", "freq.3.individual.varR", "freq.1.individual.varMut", "diff.1.ratio.varR", "diff.2.individual.varBeta",
                 "diff.2.individual.varR", "accel.ratio.varR", "accel.ratio.mutation", "freq.2.ratio.varR", "freq.3.ratio.varR",
                 "freq.1.ratio.varR", 'diff.1.meanLoad', "freq.1.ratio.mutation", "freq.2.ratio.mutation", "freq.2.ratio.meanBeta",
                 "diff.2.ratio.mutation", "diff.2.ratio.varR", "diff.1.ratio.varSigma", "diff.1.ratio.mutation", "freq.1.entropy",
                 "freq.2.entropy", "freq.3.entropy", "accel.diversity", "accel.antigenicDiversity", "freq.1.varSigma", "accel.ratio.varSigma",
                 "diff.1.entropy", "diff.1.dominant.freq", "diff.1.antigenicTypes")

## Exclusion based on collinearity after everything has been combined and up until VIF < 5 and
# incorporating the metric of I/N









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
