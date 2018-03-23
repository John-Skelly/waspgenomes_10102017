#loads libraries
library(data.table)
library(ggplot2)
library("stringi")

#input and output
mercount_file <- snakemake@input[["mercount_file"]]
dmin_plot <- snakemake@output[["dmin_plot"]]
dmin_out <- snakemake@output[["dmin_out"]]
log_file <- snakemake@log[["log"]]

#debug
#mercount_file <- "output/meraculous/MA3/trim_decon/k_75/diplo_0/meraculous_mercount/mercount.hist"


#set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#adds identifier parameters to name 
names(mercount_file) <- stri_replace_all_regex(mercount_file, 
                                                '(.+)(/)(.+)(/)(.+)(/)(.+)(/)(.+)(/)(.+)(/)meraculous_mercount.+$', 
                                                replacement = '$5 $7 $9 $11'  )

#reads in file
mercount_data <- fread(mercount_file, col.names = c('depth', 'kmers'))

#adds name coloumn 
mercount_plot <- rbind2(mercount_data, idcol = 'param')

#returns dmin, constructs graph and indicates dmin
TRUE_kmers<- (which(((diff(diff(mercount_plot$kmers) >0) !=0))))
dmin <- TRUE_kmers[1]
dmin_plot <- ggplot(mercount_plot, aes(x = depth, y = kmers, colour = "red"))+
        geom_path()+
        xlim(c(0, 200))+
        scale_y_log10()+
        geom_vline(xintercept= dmin, linetype= "dashed", color= "blue")

#write output
ggsave(filename = "dmin_plot.pdf",
       plot = dmin_plot,
       device = cairo_pdf,
       width =10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
