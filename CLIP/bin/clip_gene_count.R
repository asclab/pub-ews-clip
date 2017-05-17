# CLIP gene count plot
library(ggplot2)
bins <- read.table("results-refseq/EWS.genecounts.txt", header=F, comment.char = "#", sep = '\t')
colnames(bins) <- c('geneid', 'genename', 'ews1', 'ews2')

filtered <- bins[bins$ews1>0 & bins$ews2>0,]

ews1_total <- sum(filtered$ews1) / 1000000
ews2_total <- sum(filtered$ews2) / 1000000

filtered$ews1cpm <- log2((filtered$ews1) / ews1_total)
filtered$ews2cpm <- log2((filtered$ews2) / ews2_total)

cor(filtered$ews1cpm, filtered$ews2cpm)

data2.labels <- data.frame(
    x = c(12.5), 
    y = c(1), 
    label = c(paste0("R=",round(cor(filtered$ews1cpm, filtered$ews2cpm),5)))
)


scatterplot_p <- ggplot(filtered, aes(x=ews2cpm,y=ews1cpm)) + 
    geom_point() +
    theme_bw() +
    coord_cartesian(xlim = c(-0.5, 15), ylim = c(-0.5, 15)) +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="EWS-2 CPM (log2)", y="EWS-1 CPM (log2)") +
    geom_text(data=data2.labels, aes(x=x,y=y,label=label))
    
scatterplot_p
ggsave("figures/clip-gene-count.pdf", plot=scatterplot_p, height=6, width=7)







