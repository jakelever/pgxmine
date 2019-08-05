# Solution from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them

list.of.packages <- c('plotly', 'shiny', 'DT', 'plyr', 'dplyr', 'reshape2', 'RColorBrewer', 'data.table', 'stringr')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
