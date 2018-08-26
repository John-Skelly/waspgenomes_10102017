#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)


# busco reading function
ReadBuscoResults <- function(filepath){
                         fread(filepath,
                         header = TRUE,
                         skip = 4,
                         fill = TRUE,
                         na.strings = c(""))
}

output_rds <- snakemake@output[["rds"]]

output_plot <- snakemake@output[["plot"]]

busco_result_files <- snakemake@input[["busco_targets"]]

names(busco_result_files) <- gsub(".*/busco/(.+)/run_busco/.*", "\\1", busco_result_files)

busco_result_list <- lapply(busco_result_files, ReadBuscoResults)

busco_results_combined <- rbindlist(busco_result_list, idcol = "assembly")

# count the number in each status
busco_stats <- busco_results_combined[, .(
  n_buscos = length(unique(`# Busco id`))),
  by = .(Status, assembly)]

busco_stats[,total_busco:=sum(n_buscos), by= assembly]

busco_stats[,percent_complete:=n_buscos*100/total_busco, by= .(Status, assembly)]

busco_stats[,c("strain", "read_set", "k", "diploid_mode"):= tstrsplit(assembly, "/"), by= assembly]

#reorder table

busco_stats[,k:=factor(k,levels=c("k_63", "k_67", "k_71", "k_75", "k_79")) ]

#plot

gp <- ggplot(busco_stats,aes(x=k, y=percent_complete, fill=Status))+
  geom_col(position = "dodge")+
  facet_grid(read_set+diploid_mode~strain)

saveRDS(busco_stats, file = output_rds)

ggsave(filename = output_plot,plot = gp, width = 10, height = 7.5, units = 'in')
ggsave(filename = "output_plot.png",plot = gp, width = 10, height = 7.5, units = 'in', device = "png")
