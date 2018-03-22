install.packages("dplyr")
install.packages("tidyr")
library(data.table)
library(ggplot2)
library(tidyr)
library("stringi")

mercount_files <- list.files('output', 'mercount.hist', recursive = TRUE, full.names = TRUE)

names(mercount_files) <- stri_replace_all_regex(mercount_files, '(.+)(/)(.+)(/)(.+)(/)(.+)(/)(.+)(/)(.+)(/)meraculous_mercount.+$', replacement = '$5 $7 $9 $11'  )

mercount_list <- lapply(mercount_files, fread, col.names = c('depth', 'kmers'))

mercount_plot <- rbindlist(mercount_list, idcol = 'param')

mercount_plot_param <- split(mercount_plot, mercount_plot$param)

for (assemblies in mercount_plot_param)

MA3_trim_decon_k_127_diplo1 <- data.frame(mercount_plot_param["MA3 trim_decon k_127 diplo_1"])

names(MA3_trim_decon_k_127_diplo1) <- c("param", "depth", "kmers")

MA3_trim_decon_k_127_diplo1

ggplot(MA3_trim_decon_k_127_diplo1, aes(x = depth, y = kmers, colour = params))+
  geom_path()+
  xlim(c(0, 200))+
  scale_y_log10()#+
  #facet_wrap(~params, scales = 'free_y')+
  #geom_vline(xintercept = c(19, 38))


mercount_plot_param[param == 'trim_decon', diff(diff(kmers) >0) !=0]

mercount_plot[readset == 'norm', diff(diff(kmers) >0) !=0]

x <- c(19, 38)
y <- c(151)
x * (y/(y-71+1))
x/(71/(71-31+1))
