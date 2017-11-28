#!/usr/bin/env Rscript

library(data.table)
library(bit64)
library(ggplot2)

hist_files<- list.files("output/norm", pattern = "hist", full.names = TRUE)

names(hist_files) <- sub(".txt", "", basename(hist_files))
hist_data_list <- lapply(hist_files, fread, colClasses = c("integer","integer64", "integer64"))

hist_data <- rbindlist(hist_data_list, idcol = "filename")

hist_data[,strain := gsub("_.*", "" , filename)]

hist_data[grep("hist_out", filename), hist := "after"]

hist_data[grep("hist$", filename), hist := "before"]

ggplot(hist_data,mapping=aes(x = log(`#Depth`,4), y = Unique_Kmers, colour=hist))+
  facet_wrap(~strain)+
  scale_y_log10()+
  geom_path()
