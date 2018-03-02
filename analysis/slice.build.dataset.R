##### Analysis Code for Random Point Surveillance 
rm(list=ls())
source('loadLibraries.R')
source('analysis_functions.R')
source('plotting_functions.R')


# Step 1. Put together list trials of transient and successful antigens ; this will be based on max.freq, and days above
tropics.folder = "../data/tropics/eligible/"
trial.dirs = dir(tropics.folder)

# This is a slow step
antigen.frequencies = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.antigenFrequencies.txt"), header = TRUE)) 
names(antigen.frequencies) = trial.dirs

######## Step 2: For each list entry need to calculate max frequency and days above
thres = .2

# Calculating the two criteria for determining success and if successful how long was it there 
antigen.frequencies %>%
  map(function(x) calculate_days_above(x, threshold = .2)) -> days.above

antigen.frequencies %>%
  map(find_max_frequency) %>%
  map(function(x) mutate_at(x, "antigentype", as.character)) -> max.frequencies

# Creating the combined data frame and using the criteria to assign a label 
rel.frequency.thres = .01 # Threshold at which we start monitoring/expect to pick up in surveillance 
full.data = map2(max.frequencies, days.above, left_join) %>%
  map(replace_na_zeros) %>% 
  map(determine_success_labels, max.rf = rel.frequency.thres)

####### Step 3: Subsetting for just transient and yes; and getting rid of zeros 
subset.data = map(full.data, filter_out_loss) 
subset.df = do.call("rbind", subset.data)
subset.df$name = rep(trial.dirs, sapply(subset.data, nrow))

subset.df %>%
  filter(antigentype != 0) -> subset.analyze

######### Step 4 - Calculate Entropy
entropy = map(antigen.frequencies, calculate_entropy)
names(entropy) = trial.dirs

######### Step 5 Generate Sample Times
set.seed(103114)
possible.sample = seq(from = 1833, to = 10947, by = 14)
sample.time.one = sample(possible.sample, 10, replace=FALSE) 
sample.time.two = adply(sample.time.one, .margins = 1, generate_second_time,.id = NULL)
sample.times = data.frame(t1 = sample.time.one, t2= sample.time.two$V1) %>% arrange(t1)

######################## Function to collect population level at those times
 # For trial time 1 
data.1.list = map(trial.dirs, create_trial_dataset, sample.times$t1)
data.1.df = map(data.1.list, create_time_variable) %>% 
  map(function(x) x %>% mutate(time.stamp = 't1')) %>%
  dplyr::bind_rows(.id = 'list')

# When just doing one 
#remove.variables = c("date", "day", "day.success", "stage", "simDay")
#data.1.df %<>%
#  filter(!(variable %in% remove.variables))

# What's present at 1 
d1.combos.present = data.1.df %>%
  select(antigentype, success, group, trial) %>%
  distinct(antigentype, success, group, trial) %>%
  mutate(time.1 = "t1")

###############################
data.2.list = map(trial.dirs, create_trial_dataset, sample.times$t2)
data.2.df = map(data.2.list, create_time_variable) %>% 
  map(function(x) x %>% mutate(time.stamp = 't2')) %>%
  dplyr::bind_rows(.id = 'list')

# What's present at time 2
d2.combos.present = data.2.df %>%
  select(antigentype, success,group, trial) %>%
  distinct(antigentype, success,group,trial) %>%
  mutate(time.2 = "t2")

# Select only what is present at both 
full_join(d1.combos.present, d2.combos.present) %>%
  filter(!(is.na(time.1) | is.na(time.2))) %>%
  mutate(id = paste0(antigentype, "_", trial)) -> full.data.subset

# Filter data set based on what is present in both 
data.time.combined = rbind(data.1.df, data.2.df)
data.time.combined %>%
  mutate(id = paste0(antigentype, "_", trial)) %>%
  filter(id %in% full.data.subset$id) -> data.present

########### Place back into list format based on trial  
# Back into list format based on trial 
full.data.list = dlply(data.present, "trial", identity)

remove.variables = c("day.success", "stage", "simDay", "date")

## Calculate diff for each trial   
list.calculate.diff = map(full.data.list, calculate_diff)

## Combine 
full.data.df = do.call("rbind", list.calculate.diff)
