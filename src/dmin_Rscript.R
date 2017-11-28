library(data.table)
library(ggplot2)

mercount_files <- list.files('output', 'mercount.hist', recursive = TRUE, full.names = TRUE)

names(mercount_files) <- gsub(pattern = '^.+MA3/(.+)/meraculous_mercount.+$', replacement = '\\1' , mercount_files)

mercount_list <- lapply(mercount_files, fread, col.names = c('depth', 'kmers'))

mercount_plot <- rbindlist(mercount_list, idcol = 'readset')

ggplot(mercount_plot, aes(x = depth, y = kmers, colour = readset))+
  geom_path()+
  xlim(c(0, 200))+
  scale_y_log10()+
  facet_wrap(~readset, scales = 'free_y')+
  geom_vline(xintercept = c(19, 38))

mercount_plot[readset == 'trim_decon', diff(diff(kmers) >0) !=0]

mercount_plot[readset == 'norm', diff(diff(kmers) >0) !=0]

x <- c(19, 38)
y <- c(151)
x * (y/(y-71+1))
x/(71/(71-31+1))
