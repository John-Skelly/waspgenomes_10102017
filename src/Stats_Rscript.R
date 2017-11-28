library(data.table)
library(ggplot2)
library(bit64)

stats <- fread('output/assembly_stats/stats.txt')

stats[ ,assembly_ID:= gsub(".*/meraculous/(.+)/meraculous_gap_closure.*", "\\1", filename)]

stats[ , c("strain", "read_set", "k", "diploid_mode") := tstrsplit(assembly_ID, "/")]

plot_data <- melt(stats,id.vars = c("strain", "read_set", "k", "diploid_mode") , measure.vars = c("n_scaffolds", "scaf_bp", "scaf_N50"))

ggplot(plot_data, aes(x = diploid_mode, y = value, fill = read_set))+
  facet_grid(strain+variable ~ k, scales = "free_y")+
  geom_col(position = "dodge")
