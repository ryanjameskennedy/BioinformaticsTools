#!/usr/bin/env Rscript

rm(list=ls())
graphics.off()

install.packages("ggplot2", repos = "https://cloud.r-project.org")
install.packages("viridis", repos = "https://cloud.r-project.org")
install.packages("optparse", repos = "https://cloud.r-project.org")
library(ggplot2)
library(viridis)
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="dataset file path", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="out.png", help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).", call.=FALSE)
}

plot_ani_dddh <- function(ani_dddh, out_fpath) {
    ggplot(ani_dddh, aes(x=ANI, y=DDH.1, fill=Relationship)) +
        geom_point(aes(size=Relationship), stroke=0.3, color="black", shape = 21, position = position_jitter(width = 0.5, height = 0.5)) +
        geom_hline(yintercept=70, linetype="dashed", color = "grey") +
        geom_vline(xintercept=96, linetype="dashed", color = "grey") +
        geom_text(aes(96, 70, label="Species boundaries", hjust = 1.05, vjust = -0.5), color="grey", size=4.5) +
        geom_text(aes(min(ANI), 70, label="y = 70", vjust = -0.5), color="grey", size=4) +
        geom_text(aes(96, min(DDH.1), label="x = 96", angle = 90, vjust = -0.5), color="grey", size=4) +
        ylab("Digital DNA-DNA Hybridisation") +
        xlab("Average Nucleotide Identity") +
        scale_fill_manual(values=c('darkgrey','orange', 'lightblue', 'seagreen')) +
        #scale_fill_viridis(discrete = TRUE,  begin = 0.2, end = 0.85, direction = -1, option="viridis") +
        scale_size_manual(values=c(3,3,3,3)) +
        theme_bw() +
        theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            legend.title = element_text(size=13),
            legend.text=element_text(size=12))
    
    ggsave(out_fpath, width=9, height=6, dpi=600)
}

filepath <- opt$input
ani_dddh <- read.csv(filepath, header = TRUE, sep = ",", row.names = NULL)
out_fpath <- paste(tools::file_path_sans_ext(opt$output), ".png", sep="")
plot_ani_dddh(ani_dddh, out_fpath)

#dirs <- list.dirs(path=".", full.names=TRUE, recursive=FALSE)
#for (dir_fpath in dirs) {
#  isolate <- unlist(strsplit(dir_fpath, "/"))[-1]
#  filepath <- paste(dir_fpath,"/", isolate,".csv", sep="")
#  ani_dddh <- read.csv(filepath, header = TRUE, sep = ",", row.names = NULL)
#  out_fpath <- paste(dir_fpath,"/", isolate,".png", sep="")
#  plot_ani_dddh(ani_dddh, out_fpath) # apply function
#}

graphics.off()