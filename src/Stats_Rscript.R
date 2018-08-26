library(data.table)
library(ggplot2)
library(bit64)

stats <- fread('output/assembly_stats/stats.txt')

stats[ ,assembly_ID:= gsub(".*/meraculous/(.+)/meraculous_gap_closure.*", "\\1", filename)]

stats[ , c("strain", "read_set", "k", "diploid_mode") := tstrsplit(assembly_ID, "/")]

setkey(stats, strain, read_set, k, diploid_mode)

stats_with_NA <- stats[CJ(unique(strain), unique(read_set), unique(k), unique(diploid_mode))]

stats_with_NA[,k:=factor(k,levels=c("k_31", "k_63", "k_67", "k_71", "k_75", "k_79", "k_127")) ]

plot_data_all <- melt(stats_with_NA,id.vars = c("strain", "read_set", "k", "diploid_mode") , measure.vars = c("n_scaffolds", "scaf_bp", "scaf_N50"))

plot_data_n <- melt(stats_with_NA,id.vars = c("strain", "read_set", "k", "diploid_mode") , measure.vars = c("n_scaffolds"), fill = TRUE)

plot_data_bp <- melt(stats_with_NA,id.vars = c("strain", "read_set", "k", "diploid_mode") , measure.vars = c("scaf_bp"))

plot_data_N50 <- melt(stats_with_NA,id.vars = c("strain", "read_set", "k", "diploid_mode") , measure.vars = c( "scaf_N50"))

plot_data_L50 <- melt(stats_with_NA,id.vars = c("strain", "read_set", "k", "diploid_mode") , measure.vars = c( "scaf_L50"))

ggplot(plot_data_all, aes(x = diploid_mode, y = value, fill = read_set, label=value))+
  facet_grid(strain+variable ~ k, scales = "free_y")+
  geom_col(position = "dodge")+
  geom_text(size = 3, check_overlap = TRUE, position = position_stack(0.5)
            )
  

ggplot(plot_data_n, aes(x = diploid_mode, y = value, fill = read_set, label=value))+
  facet_grid(strain+variable ~ k, scales = "free_y")+
  geom_col(position = "dodge")+
  geom_text(size = 3, check_overlap = TRUE, position = position_stack(0.1)
            )
  

ggplot(plot_data_bp, aes(x = diploid_mode, y = value, fill = read_set, label=scaf_bp, label=value))+
  facet_grid(strain+variable ~ k, scales = "free_y")+
  geom_col(position = "dodge")+
  geom_text(size = 3, check_overlap = TRUE, position = position_stack(0.5)
          )

ggplot(plot_data_N50, aes(x = diploid_mode, y = value, fill = read_set, label=scaf_N50, label=value))+
  facet_grid(strain+variable ~ k, scales = "free_y")+
  geom_col(position = "dodge")+
  geom_text(size = 3, check_overlap = TRUE, position = position_stack(0.5)
          )

ggplot(plot_data_L50, aes(x = diploid_mode, y = value, fill = read_set, label=scaf_L50, label=value))+
  facet_grid(strain+variable ~ k, scales = "free_y")+
  geom_col(position = "dodge")+
  geom_text(size = 3, check_overlap = TRUE, position = position_stack(0.1)
          )

