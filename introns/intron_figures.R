# CLEAR ENVIRONMENT -----------------------------------------------------------------------------------------------------------------------

rm(list=ls())


# CLEAR GRAPHICS --------------------------------------------------------------------------------------------------------------------------

graphics.off()


# INSTALL PACKAGES ------------------------------------------------------------------------------------------------------------------------

install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)


# READ IN DATA ----------------------------------------------------------------------------------------------------------------------------

setwd("src/")
intron_counts <- read.csv("../data/stats/intron_counts.csv", header = TRUE, sep = ",", row.names = NULL)


# CREATE OUTPUT FOLDER --------------------------------------------------------------------------------------------------------------------

dir.create("../data/figures")


# TRANSFORM DATA -------------------------------------------------------------------------------------------------------------------------

intron_counts$Label <- paste0("n=",intron_counts$Count)


# PLOTS -----------------------------------------------------------------------------------------------------------------------------------

##Stacked barplot of seed counts
ggplot(intron_counts, aes(x=Species, y=Count, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=Label), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5)+
  ylab("Counts") +
  scale_fill_manual(values=c("#2F5597", "#81D4FA")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text=element_text(size=10))

ggsave("../data/figures/intron_counts.png", width=9, height=6, dpi=600)

dev.off()
