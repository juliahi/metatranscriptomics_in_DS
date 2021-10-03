library ('ggplot2')
require(gridExtra)
library(egg)
data <- read.csv('metagenemark_counts_filtered_genes_all_summed_normalized.csv', stringsAsFactors = F)
data$X<-data$X.1 <- data$FoldChange.WT.DN <- data$FoldChange.T1.T3 <- NULL

timepoints<-c('T2', 'T1', 'T2', 'T1', 'T2', 'T1', 'T2', 'T1', 'T2')
# Timepoints

rozm=2

rysuj_lewy <- function (gene) {
  new_data <- data.frame('gene' = t(data[data$gene==gene,c(2:ncol(data))]), 'timepoints' = timepoints) 
  colnames(new_data) <- c('expression', 'timepoints')
  ggplot(new_data, aes(x = timepoints, y = expression, color=timepoints)) + 
    geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
    scale_color_manual(values=c("#ffb90f", "#B8860B")) +
    geom_jitter(height=0.2, width = 0.2, size=rozm) +
    xlab(' ') + ylab('Read counts (normalized)') +ggtitle(paste0(gene,' gene')) +
    theme(legend.title  = element_blank(), axis.text=element_text(size=12)) +
    theme(legend.position = "none") +
    theme(panel.background = element_rect(fill = "white",colour = "white"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) 
}


rysuj_srodek <- function (gene) {
  new_data <- data.frame('gene' = t(data[data$gene==gene,c(2:ncol(data))]), 'timepoints' = timepoints) 
  colnames(new_data) <- c('expression', 'timepoints')
  ggplot(new_data, aes(x = timepoints, y = expression, color=timepoints)) + 
    geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
    scale_color_manual(values=c("#ffb90f", "#B8860B")) +
    geom_jitter(height=0.2, width = 0.2, size=rozm) +
    xlab(' ') + ylab(' ') +ggtitle(paste0(gene,' gene')) +
    theme(legend.title  = element_blank(), axis.text=element_text(size=12)) +
    theme(legend.position = "none") +
    theme(panel.background = element_rect(fill = "white",colour = "white"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) 
}

rysuj_prawy <- function (gene) {
  new_data <- data.frame('gene' = t(data[data$gene==gene,c(2:ncol(data))]), 'timepoints' = timepoints) 
  colnames(new_data) <- c('expression', 'timepoints')
  ggplot(new_data, aes(x = timepoints, y = expression, color=timepoints)) + 
    geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
    scale_color_manual(values=c("#ffb90f", "#B8860B")) +
    geom_jitter(height=0.2, width = 0.2, size=rozm) +
    xlab(' ') + ylab(' ') +ggtitle(paste0(gene,' gene')) +
    theme(legend.title  = element_blank(), axis.text=element_text(size=12)) +
    theme(panel.background = element_rect(fill = "white",colour = "white"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) 
}

p1 <- rysuj_lewy('LTAS')
p2 <- rysuj_srodek ('PROWX')
p3 <- rysuj_prawy('DHBF')



pdf('DE_genes_dotplots_timepoints.pdf', 6, 6, onefile=FALSE)
ggarrange(p1, p2, p3, ncol=3, bottom = 'Timepoints of HFD')
dev.off()

################################
genotypes<-c('WT', 'WT', 'WT', 'WT', 'WT', 'Ts65Dn', 'Ts65Dn', 'Ts65Dn', 'Ts65Dn')

# Genotypes
## Hypoxanthine
rozm <- 2

rysuj_lewy <- function (gene) {
  new_data <- data.frame('gene' = t(data[data$gene==gene,c(2:ncol(data))]), 'genotypes' = genotypes) 
  colnames(new_data) <- c('expression', 'genotypes')
  ggplot(new_data, aes(x = genotypes, y = expression, color=genotypes)) + 
    geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
    scale_color_manual(values=c("#ff3030", "#00bfff")) +
    geom_jitter(height=0.2, width = 0.2, size=rozm) +
    xlab(' ') + ylab('Read counts (normalized)') +ggtitle(paste0(gene,' gene')) +
    theme(legend.position = "none") +
    theme(legend.title  = element_blank(), axis.text=element_text(size=10)) +
    theme(panel.background = element_rect(fill = "white",colour = "white"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

rysuj_srodek <- function (gene) {
  new_data <- data.frame('gene' = t(data[data$gene==gene,c(2:ncol(data))]), 'genotypes' = genotypes) 
  colnames(new_data) <- c('expression', 'genotypes')
  ggplot(new_data, aes(x = genotypes, y = expression, color=genotypes)) + 
    geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
    scale_color_manual(values=c("#ff3030", "#00bfff")) +
    geom_jitter(height=0.2, width = 0.2, size=rozm) +
    xlab(' ') + ylab(' ') +ggtitle(paste0(gene,' gene')) +
    theme(legend.title  = element_blank(), axis.text=element_text(size=10)) +
    theme(legend.position = "none") +
    theme(panel.background = element_rect(fill = "white",colour = "white"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

rysuj_prawy <- function (gene) {
  new_data <- data.frame('gene' = t(data[data$gene==gene,c(2:ncol(data))]), 'genotypes' = genotypes) 
  colnames(new_data) <- c('expression', 'genotypes')
  ggplot(new_data, aes(x = genotypes, y = expression, color=genotypes)) + 
    geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
    scale_color_manual(values=c("#ff3030", "#00bfff")) +
    geom_jitter(height=0.2, width = 0.2, size=rozm) +
    xlab(' ') + ylab(' ') +ggtitle(paste0(gene,' gene')) +
    theme(legend.title  = element_blank(), axis.text=element_text(size=10)) +
    theme(panel.background = element_rect(fill = "white",colour = "white"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray90")) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
}

p1 <- rysuj_lewy ('HBAC')
p1
p2 <- rysuj_srodek('YGET')
p2
p3 <- rysuj_srodek('XDHC')
p3
p4 <- rysuj_srodek('ALLD')
p4
p5 <- rysuj_prawy('YLBA')
p5



pdf('DE_genes_dotplots_genotypes.pdf',6,6,onefile=FALSE)
ggarrange(p1, p2, p3, p4, p5, ncol=5, bottom = 'Mouse strain', top = 'DE genes involved in hypoxanthine metabolism')
dev.off()
