FDR_thres <- 0.10
bins <- read.table("results-refseq/EWS.enriched.bins.txt", header=TRUE, comment.char = "#", sep = '\t')

nrow(bins)
sum(bins$EWS1.count)
sum(bins$EWS2.count)

bins <- bins[bins$region_used != 'INTERGENIC',]
bins <- bins[bins$region_used != 'MITOCHONDRIAL',]
bins <- bins[bins$region_used != '',]
bins<-bins[grep('_ANTI', bins$region_used, invert = T),]

# try a pre-filter?
nrow(bins)
sum(bins$EWS1.count)
sum(bins$EWS2.count)

prefiltered_bins <- bins

bins<-bins[bins$EWS1.count >=5 & bins$EWS2.count >=5,]

nrow(bins)
sum(bins$EWS1.count)
sum(bins$EWS2.count)


# expected values should be at least "1"
adj_min_one<-function(row, cond) {
    max(1,as.numeric(row[[cond]]))
}

bins$ews1_expected_region_alt <- apply(bins, 1, adj_min_one, 'ews1_expected_region')
bins$ews2_expected_region_alt <- apply(bins, 1, adj_min_one, 'ews2_expected_region')


# Calc p-values based on Poisson prop that *actual* is greater than expected
bins$ews1_region_pvalue<-1-ppois(bins$EWS1.count-1, bins$ews1_expected_region_alt)
bins$ews2_region_pvalue<-1-ppois(bins$EWS2.count-1, bins$ews2_expected_region_alt)

# Calc FDRs
bins$ews1_region_FDR<-p.adjust(bins$ews1_region_pvalue, 'BH')
bins$ews2_region_FDR<-p.adjust(bins$ews2_region_pvalue, 'BH')

# require a bin to be in both samples.
filtered_region<-bins[bins$ews1_region_FDR < FDR_thres & bins$ews2_region_FDR < FDR_thres,]
background<-bins[bins$ews1_region_FDR >= FDR_thres | bins$ews2_region_FDR >= FDR_thres,]

summary(bins$region_used)
nrow(bins)
length(unique(bins$gene_used))

summary(prefiltered_bins$region_used)
nrow(prefiltered_bins)
length(unique(prefiltered_bins$gene_used))

summary(background$region_used)
nrow(background)
length(unique(background$gene_used))


summary(filtered_region$region_used)
nrow(filtered_region)
length(unique(filtered_region$gene_used))

(summary(filtered_region$region_used)/nrow(filtered_region))/(summary(background$region_used)/nrow(background))

(summary(filtered_region$region_used)/nrow(filtered_region))/(summary(prefiltered_bins$region_used)/nrow(prefiltered_bins))


write.table(unique(filtered_region$gene_used), 'results-refseq/EWS.filtered.genes.region.txt', quote=F, sep='\t', row.names=F, na='', col.names=F)
write.table(filtered_region, 'results-refseq/EWS.filtered.bins.txt', quote=F, sep='\t', row.names=F, na='')
write.table(background, 'results-refseq/EWS.background.bins.txt', quote=F, sep='\t', row.names=F, na='')

