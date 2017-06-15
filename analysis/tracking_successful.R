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
