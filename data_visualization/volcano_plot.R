library(ggplot2)
library(dplyr)

# Data here: http://regulomics.mimuw.edu.pl/~ilona/metatranscriptomics_in_DS/DeSeq2/

#setwd('/mnt/chr3/People/Ilona/Analyses_microbiomes_gene_exprs/2021_10_02_volcano_plots/')

# WT vs DS volcano plot

data <- read.csv('wt_tri_pvals.csv', stringsAsFactors = F, sep = ',')
data$color <- 'not_significant'
data$color [data$padjFDR<0.05 & data$log2fold < 0] <- 'Ts65Dn_UpRegulated'
data$color [data$padjFDR<0.05 & data$log2fold > 0] <- 'WT_UpRegulated'
data$log10padjFDR <- log10(data$padjFDR)
ggplot(data, aes(x=log2fold, y=-log10padjFDR, color=color)) +
  geom_point(  ) +
  scale_color_manual("", values = c("not_significant" = 'grey', 
                                    "Ts65Dn_UpRegulated" = "#ff3030", 
                                    "WT_UpRegulated" = "#00bfff"),
                     labels = c("not significant","Ts65Dn up-regulated", "WT up-regulated")) +
  theme(panel.background = element_rect(fill = "white",colour = "white"), 
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) +
   ylab('-log10 (p adj FDR)')

ggsave('volcano_plot_WT_DS.png', width=4, height = 4)
ggsave('volcano_plot_WT_DS.svg', width=4, height = 4)

# T1 vs T2 volcano plot
colnames(data)

data <- read.csv('all4_all16_pvals.csv', stringsAsFactors = F, sep = ',')
data <- data %>% arrange(padjFDR)
data$color <- 'not_significant'
data$color [1:100] <- 'significant'
data$color [data$color=='significant' & data$log2fold > 0] <- 'T1_UpRegulated'
data$color [data$color=='significant' & data$log2fold < 0] <- 'T2_UpRegulated'

data$log_parpv <- log(data$DESeq2parpv)
ggplot(data, aes(x=log2fold, y=-log_parpv, color=color)) +
  geom_point(  ) +
  scale_color_manual("", values = c("not_significant" = 'grey', 
                                    "T1_UpRegulated" = "#ffb90f", 
                                    "T2_UpRegulated" = "#B8860B"),
                     labels = c("not significant","T1 up-regulated", "T2 up-regulated")) +
  theme(panel.background = element_rect(fill = "white",colour = "white"), 
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) +
  ylab('-log (p)') 

ggsave('volcano_plot_T1_T2.png', width=4, height = 4)
ggsave('volcano_plot_T1_T2.svg', width=4, height = 4)
