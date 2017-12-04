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

#Globals

busco_result_files <- list("output/busco/IE/norm/k_127/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/IE/norm/k_31/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/IE/norm/k_71/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/IE/norm/k_71/diplo_1/run_busco/full_table_busco.tsv",
"output/busco/IE/trim_decon/k_127/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/IE/trim_decon/k_31/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/IE/trim_decon/k_31/diplo_1/run_busco/full_table_busco.tsv",
"output/busco/IE/trim_decon/k_71/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/MA3/norm/k_71/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/MA3/trim_decon/k_71/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/FR1/norm/k_31/diplo_1/run_busco/full_table_busco.tsv",
"output/busco/FR1/norm/k_71/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/FR1/norm/k_71/diplo_1/run_busco/full_table_busco.tsv",
"output/busco/FR1/trim_decon/k_31/diplo_0/run_busco/full_table_busco.tsv",
"output/busco/FR1/trim_decon/k_31/diplo_1/run_busco/full_table_busco.tsv")

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

busco_stats[,k:=factor(k,levels=c("k_31", "k_71", "k_127")) ]

#plot

ggplot(busco_stats,aes(x=k, y=percent_complete, fill=Status))+
  geom_col(position = "dodge")+
  facet_grid(read_set+diploid_mode~strain)
  