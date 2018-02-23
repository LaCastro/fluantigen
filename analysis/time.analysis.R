library(plyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(viridis)
library(cowplot)


data.dir = "~/Dropbox/current_fluantigen/model_csvs/"

term.files = dir(paste0(data.dir, "recent/three.point.complete/"))
term.list = map(term.files, function(x) read.csv(paste0(data.dir, "recent/three.point.complete/", x), header = TRUE))
term.df = do.call("rbind", term.list)
term.df$name =  rep(c(term.files), sapply(term.list, nrow))

term.df %>%
  separate(col = name, into = c("rf", "file"),  sep = ".csv") %>%
  separate(col = rf, into = c("date", "rf"), sep = "021318complete.") %>%
  select(-date, -file) -> term.df


lead.files = dir(paste0(data.dir, "lead.time/upto10/"))
lead.list = map(lead.files, function(x) read.csv(paste0(data.dir, "lead.time/upto10/", x), header = TRUE))
lead.df = do.call("rbind", lead.list)
lead.df$name =  rep(c(lead.files), sapply(lead.list, nrow))
lead.df %>%
  separate(col = name, into = c("rf", "file"),  sep = ".csv") %>%
  separate(col = rf, into = c("date", "rf"), sep = "021418.upto10.lead.") %>%
  select(-date, -file) -> lead.df


roc.files = dir(paste0(data.dir, "time.roc/three.point.complete/"))
roc.list = map(roc.files, function(x) read.csv(paste0(data.dir, "time.roc/three.point.complete/", x), header = TRUE))
roc.df = do.call("rbind", roc.list)
roc.df$name =  rep(c(roc.files), sapply(roc.list, nrow))
roc.df %>%
  separate(col = name, into = c("rf", "file"),  sep = ".csv") %>%
  separate(col = rf, into = c("date", "rf"), sep = "021318complete.roc.") %>%
  select(-date, -file) -> roc.df



################ First Figure
## Calculating average lead times

lead.df %>%
  group_by(rf) %>%
  summarize(low.qr = quantile(lead.time, probs = .25),
            median = median(lead.time),
            high.qr = quantile(lead.time, probs = .75)) -> lead.summary

## Getting Performance for final and two.two days
comparison = c(2,5,10)

term.df$rf = factor(term.df$rf, levels = c(2,3,4,5,6,7,8,9,10))

term.df %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  filter(rf %in% comparison) -> term.plot

term.df %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  filter(term == 1 | term == 2 | term ==3) -> first.three

term.plot %>%
  ggplot(aes(x = term, y = per, color = as.factor(rf))) + geom_line(size = 2) + 
  geom_ribbon(aes(ymin = min.per, ymax = max.per, 
                  color = NA, fill = as.factor(rf)), alpha = .2) + 
        scale_color_viridis(discrete = TRUE) +  
  scale_fill_viridis(discrete = TRUE) + 
  scale_x_continuous(breaks = seq(1,16,2)) +
  scale_y_continuous(breaks = seq(.75,.95, .02)) + 
  guides(fill = FALSE) + 
  labs(x = "Term", y = "AUC 5 Fold CV", color = "Rel Freq. \n Thres") -> performance.terms

save_plot(performance.terms, filename = "exploratory.figs/performance.terms.pdf", base_height = 8)

term.df %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  filter(term == max(term) | term == 2 ) -> term.subset

term.subset$rf = factor(term.subset$rf, levels = c(2,3,4,5,6,7,8,9,10))
lead.summary$rf = factor(lead.summary$rf, levels = c(2,3,4,5,6,7,8,9,10))

term.subset %>%
  filter(term != 2) %>%
  left_join(lead.summary) %>%
  ggplot(aes(x=median, y = per, color = rf)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin = min.per, ymax = max.per, xmin = low.qr, xmax = high.qr)) + 
  scale_color_viridis(discrete = TRUE) + labs(x = "Days Ahead", y = "AUC", color = "Rel. Freq. \n Thres") -> lead.vs.per


save_plot("exploratory.figs/lead.vs.per.pdf", lead.vs.per, base_height = 8)

####### SEcond figure, comparing secodn term to final 
term.subset %>%
  mutate(term.type =ifelse(term =="2", "second", "max")) %>%
  select(per, term.type, rf) %>%
  ggplot(aes(x=rf, y = per, color = rf)) +geom_point(size = 3) + geom_line(size = 2) + 
  scale_color_viridis(discrete = TRUE) + labs(x = "Ref. Freq. Thres", y = "AUC") +
  guides(color = FALSE) -> improvement.per


lead.df$rf = factor(lead.df$rf, levels = c(2,4,6,8,10))
lead.df %>%
#  filter(rf != 10) %>%
  ggplot(aes(x = rf, y = lead.time)) + geom_boxplot() +
  scale_y_continuous(breaks = seq(0, 1500, 100)) +
  geom_hline(yintercept = 365, color = "purple") + geom_hline(yintercept = 185, color = "orange") +
  labs(x = "Relative Frequency Threshold", y ="Lead Time") -> lead.time.plot
save_plot(filename = "exploratory.figs/lead.time.pdf", lead.time.plot, base_height = 8)


number.points = data.frame(points = c(5001, 3874, 3224,2782,2466,2249,2090,1932,1816))

n.suc = 1218
n.trans = c(7151,4124,2951,2268,1798,1465,1238,1076,909,783)
calculate_prop = function(n.suc, n.trans) {
  n.suc/(n.suc + n.trans)
}

prop.sucess = data.frame(prop = calculate_prop(n.suc, n.trans))

term.df %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  filter(term == max(term)) %>%
  ungroup() %>%
  mutate_at("rf", as.numeric) %>%
  arrange(rf) %>%
  bind_cols(number.points) -> term.subset

term.subset$prop = prop.sucess$prop[-1]

term.subset %>%
ggplot(aes(x = points, y = per, color = as.factor(rf))) + geom_point(size = 3) + scale_color_viridis(discrete = TRUE) +
  geom_errorbar(aes(ymin = min.per, ymax = max.per)) +
  labs(y = "5 Fold CV AUC", x = "Size of Data Set", color = "Rel. Freq Thres") +
  scale_x_reverse() -> number.points.per

save_plot(number.points.per, filename = "exploratory.figs/number.points.per.pdf", base_height = 8)




term.subset %>%
  ggplot(aes(x = rf, y = per, color = as.factor(rf))) + geom_point(size = 3) + scale_color_viridis(discrete = TRUE) +
  geom_errorbar(aes(ymin = min.per, ymax = max.per)) +
  geom_line(aes(x = rf, y = prop), color = "grey") + 
  labs(y = "5 Fold CV AUC", x = "Size of Data Set", color = "Rel. Freq Thres")

############################## On its own 

term.files = dir(paste0(data.dir, "recent/up.to.10"))
term.list = map(term.files, function(x) read.csv(paste0(data.dir, "recent/up.to.10/", x), header = TRUE))
term.df = do.call("rbind", term.list)
term.df$name =  rep(c(term.files), sapply(term.list, nrow))

term.df %>%
  separate(col = name, into = c("rf", "file"),  sep = ".csv") %>%
  separate(col = rf, into = c("date", "rf"), sep = "021418") %>%
  separate(col = rf, into = c("type", "rf"), sep = "upto10.") %>%
  select(-date, -file, -type) -> term.df

single = c("2", "4", "6", "8", "10")
term.df %>%
  filter(rf %in% single) %>%
  mutate(type = "single") -> term.single

term.df %>%
    filter(!(rf %in% single)) %>%
  mutate(type = "through") -> term.through

term.full = rbind(term.through, term.single)
term.full$rf = factor(term.full$rf, levels = c("2", "4", "6", "8", "10", "through4", "through6", "through8", "through10"))

term.full %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  ggplot(aes(x = term, y = per, color = rf)) + facet_wrap(~type) + geom_line(size = 2) +
  scale_color_viridis(discrete = TRUE) +
  scale_x_continuous(breaks = seq(2, 18, 2)) + labs(x = "Term", y = "Average AUC")



########################### UP to 10 

term.files = dir(paste0(data.dir, "recent/up.to.10"))
term.list = map(term.files, function(x) read.csv(paste0(data.dir, "recent/up.to.10/", x), header = TRUE))
term.df = do.call("rbind", term.list)
term.df$name =  rep(c(term.files), sapply(term.list, nrow))

term.df %>%
  separate(col = name, into = c("rf", "file"),  sep = ".csv") %>%
  separate(col = rf, into = c("date", "rf"), sep = "021418upto10.") %>%
  select(-date, -file) -> term.df

single = c("2", "4", "6", "8", "10")
term.df %>%
  filter(rf %in% single) %>%
  mutate(type = "single") -> term.single

term.df %>%
  filter(!(rf %in% single)) %>%
  mutate(type = "through") -> term.through

term.through %>%
  filter(rf == "through10")

lead.df %>%
  group_by(rf) %>%
  summarize(low.qr = quantile(lead.time, probs = .25),
            median = median(lead.time),
            high.qr = quantile(lead.time, probs = .75)) -> lead.summary


term.single$rf = factor(term.single$rf, levels = c(2,4,6,8,10))
term.through$rf = factor(term.through$rf, levels = c("through4", "through6", "through8",as.factor(term.through$rf)[1]))
lead.summary$rf = factor(lead.summary$rf, levels = c(2,4,6,8,10))

term.single %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  left_join(lead.summary) %>%
  filter(term == max(term) | term == 2) %>%
  mutate(term.type = ifelse(term == "2", "Terms 1 + 2", "Full Model")) %>%
  ggplot(aes(x=median, y = per, color = rf)) + geom_point(size = 3) + 
  facet_wrap(~term.type) + 
  geom_errorbar(aes(ymin = min.per, ymax = max.per)) +
  geom_errorbarh(aes(xmin = low.qr, xmax = high.qr)) + 
  guides(colour = guide_legend(reverse=T)) + 
  scale_color_viridis(discrete = TRUE) + labs(x = "Days Ahead", y = "5-Fold CV AUC", color = "Rel. Freq. \n Thres") -> lead.vs.per


term.through %>%
  group_by(rf) %>%
  mutate(term = row_number()) %>%
  left_join(lead.summary) %>%
  filter(term == max(term) | term == 2) %>%
  mutate(term.type = ifelse(term == "2", "Terms 1 + 2", "Full Model")) %>%
  ggplot(aes(x=median, y = per, color = rf)) + geom_point(size = 3) + 
  facet_wrap(~term.type) + 
  geom_errorbar(aes(ymin = min.per, ymax = max.per)) +
  geom_errorbarh(aes(xmin = low.qr, xmax = high.qr)) + 
  guides(colour = guide_legend(reverse=T)) + 
  scale_color_viridis(discrete = TRUE) + labs(x = "Days Ahead", y = "5-Fold CV AUC", color = "Rel. Freq. \n Thres") -> lead.vs.per



save_plot(lead.vs.per, filename = "exploratory.figs/single10.lead.vs.per.pdf", base_height = 8)







roc.files = dir(paste0(data.dir, "time.roc/three.point.complete/"))
roc.list = map(roc.files, function(x) read.csv(paste0(data.dir, "time.roc/three.point.complete/", x), header = TRUE))
roc.df = do.call("rbind", roc.list)
roc.df$name =  rep(c(roc.files), sapply(roc.list, nrow))
roc.df %>%
  separate(col = name, into = c("rf", "file"),  sep = ".csv") %>%
  separate(col = rf, into = c("date", "rf"), sep = "021318complete.roc.") %>%
  select(-date, -file) -> roc.df


rf.desired = c(2,5,9)
roc.df %>% 
  filter(rf %in% rf.desired) %>%
  group_by(rf, thres) -> roc.df.sub

roc.df.sub %>%
  mutate_at("fold", as.numeric) %>%
  group_by(rf, thres) %>%
  ggplot(aes(x =(1-spec), y = sen, color = rf, group=rf)) + geom_line(alpha = .6, size  = 2) +
  scale_color_viridis(discrete = TRUE) +
  labs( x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)", color = "Rel. Freq. \n Thres.") -> perf.details

save_plot(perf.details, filename = "exploratory.figs/perf.details.pdf", base_height = 8)


term.df %>%
  group_by(rf) %>%
  mutate(term.number = row_number()) %>%
  mutate(per.full = per/max(per),
         added.value = ifelse(term.number == 1, per.full, per.full-lag(per.full, default = first(per.full))),
         percentage = added.value * 100) %>%
  filter(term.number  < 6) -> term.number.sub
  
term.number.sub$rf = factor(term.number.sub$rf, levels = c(2,3,4,5,6,7,8,9,10))

term.number.sub %>%
  ggplot(aes(x = rf, y = percentage, fill = as.factor(term.number))) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_viridis(discrete = TRUE) + coord_cartesian(ylim=c(.2,0)) +
  labs(x = "Rel. Freq. Thres", y = "Contribution to final AUC", fill = "Term Number") -> value.added

save_plot("exploratory.figs/value.added.pdf", value.added, base_height = 8)


########################### Looking at slice
data.dir = "~/Dropbox/current_fluantigen/current_csvs/"
slice.files = dir(paste0(data.dir, "slice/"))
slice.list = map(slice.files, function(x) read.csv(paste0(data.dir, "slice/", x), header = TRUE))
freq = c("5", "10", "full")


slice.terms = do.call("rbind", slice.list[c(1,2,3)])
slice.terms$name =  rep(c(freq), sapply(slice.list[1:3], nrow))
slice.roc = do.call("rbind", slice.list[c(4,5,6)])
slice.roc$name =rep(c(freq), sapply(slice.list[4:6], nrow))


slice.terms$name = factor(slice.terms$name, levels = c("5", "10", "full"), labels = c("5%", "10%", ">10%"))
slice.terms %<>%
  group_by(name) %>%
  mutate(term  = row_number()) 
slice.terms$term = factor(slice.terms$term, levels = c(1,2,3,4,5,6,7,8,9), 
                          labels = c("Relative R", "Frequency", "Variance in R", "mean Pop R", "Relative Variance Beta", "Relative Deleterious Mutation",
                                     "Antigenic Diversity", "8", "9"))
slice.terms %>%
  ggplot(aes(x = term, y = per, color = name, group = name)) + geom_line(size = 2) +
  scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + 
  geom_ribbon(aes(ymin = min.per, ymax = max.per, fill = name), alpha = .3) + labs(x = "", y = "AUC", color = "") + 
  guides(fill = FALSE) + theme(legend.position = c(.8,.5)) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) -> performance.slice
performance.slice

######## ROC 
slice.roc$name = factor(slice.roc$name, levels = c("5", "10", "full"), labels = c("5%", "10%", ">10%"))
slice.roc %>%
  mutate_at("fold", as.factor) %>%
#  filter(fold == 1) %>%
  ggplot(aes(x =(1-spec), y = sen, color = name, group=as.factor(fold))) + geom_point(alpha = .4, size  = 2) +
  scale_color_viridis(discrete = TRUE) + guides(color = FALSE) + scale_x_continuous(breaks = seq(0,1,.2)) + 
  scale_y_continuous(breaks = seq(0,1,.2)) + 
  labs( x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)", color = "") -> sen.spec.slice



############# Percentage that first three make up
top.3 = c("Relative R", "Frequency", "Variance in R")

slice.terms %>%
  group_by(name) %>%
  mutate(value = per/max(per),
         percent.add = value-lag(x = value,1)) %>%
  mutate(percent = ifelse(is.na(percent.add), value, percent.add)) %>%
  filter(term %in% top.3) %>%
  ggplot(aes(x = name, y = percent, fill = term)) + 
    geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) + 
    scale_fill_manual(values = c("grey64", "grey37", "black")) + labs(x = "Frequency Cut-Off",
                                                                 y = "Proportion of Final AUC",
                                                                 fill = "") -> top.three.slice
column1 = plot_grid(sen.spec.slice, top.three.slice, ncol = 1)
slice.summary = plot_grid(column1, performance.slice, nrow = 1) 
save_plot(slice.summary, filename = "exploratory.figs/summary.slice.pdf", base_height = 8, base_aspect_ratio = 1.8)
