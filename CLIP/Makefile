#ENCODE = data/polyA.filtered.bin.counts.txt.gz
#TMP = tmp-polyA
#RESULTS = results-polyA

ENCODE = data/all.filtered.bin.counts.txt.gz

#TMP = tmp-all
#RESULTS = results-all
#GTF = data/ensembl.74.hg19.gtf.gz

TMP = tmp-refseq
RESULTS = results-refseq
GTF = data/refseq.nm.gtf.gz
HG19 = /Users/mbreese/Research/references/hg19/hg19.fa

READS = results-bam
JUNC_PERMUTATIONS = 1000

all: $(TMP)/genic_bins.bed.gz $(RESULTS)/EWS.genecounts.txt $(RESULTS)/EWS.genic_gene_counts.txt \
	 $(RESULTS)/EWS.filtered.bins.bed $(RESULTS)/background_regions.txt $(RESULTS)/enriched_regions.txt \
	 $(RESULTS)/EWS.filtered.gene_bin_counts.txt $(RESULTS)/background_regions_read_count.txt \
	 $(RESULTS)/enriched_regions_read_count.txt figures/clip-gene-count.pdf $(RESULTS)/EWS.filtered.blocks.bed\
	 $(RESULTS)/EWS.filtered.block_size.txt $(RESULTS)/EWS.enriched.reads_per_bin.txt \
	 $(RESULTS)/EWS.gene.bincount.txt \
	 $(RESULTS)/EWS.background.bins.bed \
	 $(RESULTS)/EWS.bins_present.txt \
	 $(RESULTS)/EWS.reads_enriched.txt \
	 $(RESULTS)/HeLa.expressed.regions.txt \
	 $(RESULTS)/EWS.CODING.blocks.fa \
	 $(RESULTS)/EWS.UTR5_INTRON.blocks.fa \
	 $(RESULTS)/EWS.UTR5.blocks.fa \
	 $(RESULTS)/EWS.UTR3_INTRON.blocks.fa \
	 $(RESULTS)/EWS.UTR3.blocks.fa \
	 $(RESULTS)/EWS.JUNCTION.blocks.fa \
	 $(RESULTS)/EWS.CODING_INTRON.blocks.fa \
	 $(RESULTS)/EWS.CODING.blocks.fa \
	 $(RESULTS)/EWS.NONINTRON.blocks.fa \
	 $(RESULTS)/EWS.filtered.blocks.fa \
	 $(READS)/readlists/EWS.failed.reads.1.txt \
	 $(READS)/EWS.filtered.reads.dist.txt \
	 $(READS)/readlists/EWS.failed.reads.1.txt.bed.dist \
	 figures/tlss-dist.pdf \
	 figures/donor-dist.pdf \
	 figures/splice-dist.pdf \
	 $(RESULTS)/EWS.genelist.UTR3.txt \
	 $(RESULTS)/EWS.genelist.txt \
	 $(RESULTS)/EWS.genelist.bg.txt \
	 $(RESULTS)/genic_bins_by_sample.txt \
	 $(RESULTS)/EWS.filtered.block_size_multiple_gene.txt \
	 $(RESULTS)/EWS.HeLa.expressed.bin.counts.txt \
	 $(READS)/gene-donor-acceptor-tlss-overlap.svg \
	 $(READS)/EWS.filtered.reads.donor.reads.fa \
	 $(READS)/EWS.filtered.reads.acceptor.reads.fa \
	 $(READS)/EWS.filtered.reads.tlss.reads.fa \
	 $(READS)/EWS.filtered.reads.bg.reads.fa \
	 $(READS)/EWS.filtered.reads.donor.gapcount.txt \
	 $(READS)/EWS.filtered.reads.acceptor.gapcount.txt \
	 $(READS)/EWS.filtered.reads.tlss.gapcount.txt \
	 $(READS)/EWS.filtered.reads.genes.txt \
	 $(RESULTS)/EWS.non-enriched-bins.txt


	 # $(RESULTS)/EWS.filtered.merged-da.txt $(RESULTS)/EWS.background.merged-da.txt \
	 # $(RESULTS)/EWS.filtered.confirmation.20.p3.txt \
	 # $(RESULTS)/EWS.genelist.donors.txt \
	 # $(RESULTS)/EWS.filtered.UTR3.fa \
	 # $(RESULTS)/EWS.filtered.donors-500.fa \
	 # $(RESULTS)/EWS.filtered.donors-500.regions.txt \
	 # $(RESULTS)/EWS.filtered.UTR3.regions.txt

clean:

# I don't want to keep this BAM file on my computer, so
# this is commented out. The file was generated on the server
# and copied over...

#PolyA.filtered.50.counts2.txt: PolyA.filtered.bam
#	ngsutilsj bam-bins --library-rf --silent --stranded --all $< | ./sample 100000 > $@


%.filtered.50.counts.txt.gz: %.filtered.bam
	# Count the number of reads in 50bp bins
	ngsutilsj bam-bins --bins 50 --library-fr --silent --stranded --all $< | gzip > $@

%.filtered.bed.gz: %.filtered.50.counts.txt.gz
	pv $< | gunzip | ../bin/bin_count_to_bed.py | gzip > $@

%.window.counts.txt.gz: ../data/%.filtered.bam $(TMP)/EWS.windows.bed
	# Count the number of reads in the extended windows
	ngsutilsj bam-count --bed $(TMP)/EWS.windows.bed --startonly --library-fr --silent $< | gzip > $@

$(TMP)/EWS.merged.50.counts.txt.gz: data/EWS1.filtered.50.counts.txt.gz data/EWS2.filtered.50.counts.txt.gz $(ENCODE) data/IgG.filtered.50.counts.txt.gz
	# Merge counts into one file
	tabutils merge 1,2n,3n,4 5 -names EWS1,EWS2,HeLa,IgG data/EWS1.filtered.50.counts.txt.gz data/EWS2.filtered.50.counts.txt.gz $(ENCODE) data/IgG.filtered.50.counts.txt.gz | bin/filter_bins.py | bin/sample | gzip > $(TMP)/EWS.merged.50.counts.txt.gz

$(TMP)/gene_sizes.txt.gz: $(GTF)
	ngsutilsj gtf-export --genes $(GTF) | gzip > $(TMP)/gene_sizes.txt.gz

$(TMP)/EWS.merged.50.gene.txt.gz: $(TMP)/EWS.merged.50.counts.txt.gz
	# Add gene annotation to merged counts file
	ngsutilsj annotate-gtf --col-end 3 --col-strand 4 --gtf $(GTF) --gene-id --gene-name --biotype --genic-region $(TMP)/EWS.merged.50.counts.txt.gz | gzip > $(TMP)/EWS.merged.50.gene.txt.gz

$(TMP)/hg19.bg.bins.bed.gz: data/hg19.fa.fai
	bin/fai_bins.py data/hg19.fa.fai | bin/sample | gzip > $(TMP)/hg19.bg.bins.bed.gz

$(TMP)/hg19.bg.bins.gene.txt.gz: $(TMP)/hg19.bg.bins.bed.gz
	ngsutilsj annotate-gtf --noheader --bed6 --gtf $(GTF) --gene-id --gene-name --biotype --genic-region $(TMP)/hg19.bg.bins.bed.gz | gzip > $(TMP)/hg19.bg.bins.gene.txt.gz

$(RESULTS)/hg19.bg.bins.gene.summary.txt: $(TMP)/hg19.bg.bins.gene.txt.gz
	gunzip -c tmp-refseq/hg19.bg.bins.gene.txt.gz | awk '{print $$8}' | bin/merge_genic_regions.py | sort | uniq -c > $(RESULTS)/hg19.bg.bins.gene.summary.txt

$(RESULTS)/HeLa.expressed.bins.txt: $(TMP)/EWS.merged.50.gene.txt.gz $(TMP)/gene_sizes.txt.gz bin/calc_enrichment_filter.py
	# Calculate bins enriched relative to ENCODE HeLa RNA control
	pv $(TMP)/EWS.merged.50.gene.txt.gz | gunzip | bin/calc_hela_expressed.py > $(RESULTS)/HeLa.expressed.bins.txt

$(RESULTS)/EWS.bins_present.txt: $(RESULTS)/EWS.enriched.bins.txt
	pv $(RESULTS)/EWS.enriched.bins.txt | bin/ews_bins_present.py 5 > $(RESULTS)/EWS.bins_present.txt

$(RESULTS)/EWS.reads_enriched.txt: $(RESULTS)/EWS.filtered.bins.txt
	pv $(RESULTS)/EWS.filtered.bins.txt | bin/ews_reads_enriched.py > $(RESULTS)/EWS.reads_enriched.txt

$(RESULTS)/HeLa.expressed.regions.txt: $(ENCODE)
	ngsutilsj annotate-gtf --col-end 3 --col-strand 4 --gtf $(GTF) --gene-id --gene-name --biotype --genic-region $(ENCODE) | bin/hela_expressed_regions.py 100 > $(RESULTS)/HeLa.expressed.regions.txt

$(RESULTS)/EWS.HeLa.expressed.bin.counts.txt: $(TMP)/EWS.merged.50.gene.txt.gz bin/calc_ews_expressed.py
	# Calculate bins enriched relative to ENCODE HeLa RNA control
	pv $(TMP)/EWS.merged.50.gene.txt.gz | gunzip | bin/calc_ews_expressed.py > $(RESULTS)/EWS.HeLa.expressed.bin.counts.txt

$(RESULTS)/EWS.enriched.bins.txt: $(TMP)/EWS.merged.50.gene.txt.gz $(TMP)/gene_sizes.txt.gz bin/calc_enrichment_filter.py
	# Calculate bins enriched relative to ENCODE HeLa RNA control
	bin/calc_enrichment_filter.py $(TMP)/EWS.merged.50.gene.txt.gz $(TMP)/gene_sizes.txt.gz > $(RESULTS)/EWS.enriched.bins.txt

$(RESULTS)/EWS.genecounts.txt: $(TMP)/EWS.merged.50.gene.txt.gz
	# Calculate the number of reads per gene
	bin/calc_genecount.py $(TMP)/EWS.merged.50.gene.txt.gz > $(RESULTS)/EWS.genecounts.txt

$(RESULTS)/EWS.genic_gene_counts.txt: $(TMP)/EWS.merged.50.gene.txt.gz
	# Calculate the number of reads per gene
	bin/calc_genic_gene_count.py $(TMP)/EWS.merged.50.gene.txt.gz > $(RESULTS)/EWS.genic_gene_counts.txt

$(TMP)/genic_bins.bed.gz: $(TMP)/EWS.merged.50.gene.txt.gz
	bin/genic_bin_to_bed.py $(TMP)/EWS.merged.50.gene.txt.gz | gzip > $(TMP)/genic_bins.bed.gz

$(RESULTS)/EWS.filtered.bins.txt: $(RESULTS)/EWS.enriched.bins.txt # bin/enriched_filter.R
	R --vanilla < bin/enriched_filter.R > $(RESULTS)/EWS.enriched.bins.log.txt

$(RESULTS)/EWS.non-enriched-bins.txt : $(RESULTS)/EWS.filtered.bins.txt $(RESULTS)/EWS.enriched.bins.txt
	bin/extract_non_enriched_bins.py $(RESULTS)/EWS.filtered.bins.txt $(RESULTS)/EWS.enriched.bins.txt > $(RESULTS)/EWS.non-enriched-bins.txt

$(RESULTS)/EWS.background.bins.bed: $(RESULTS)/EWS.enriched.bins.txt
	# Convert filtered bins to BED files
	bin/bg_bin_to_bed.py $(RESULTS)/EWS.enriched.bins.txt > $(RESULTS)/EWS.background.bins.bed
	bedutils sort $(RESULTS)/EWS.background.bins.bed

$(RESULTS)/EWS.filtered.bins.bed: $(RESULTS)/EWS.filtered.bins.txt
	# Convert filtered bins to BED files
	bin/bin_to_bed.py $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/EWS.filtered.bins.bed
	bedutils sort $(RESULTS)/EWS.filtered.bins.bed

$(RESULTS)/EWS.filtered.bins.genebins.txt: $(RESULTS)/EWS.filtered.bins.txt
	bin/bin_to_genebins.py $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/EWS.filtered.bins.genebins.txt

$(RESULTS)/background_regions.txt: $(RESULTS)/EWS.enriched.bins.txt
	bin/extract_genic.py < $(RESULTS)/EWS.enriched.bins.txt > $(RESULTS)/background_regions.txt

$(RESULTS)/enriched_regions.txt: $(RESULTS)/EWS.filtered.bins.txt
	bin/extract_genic.py < $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/enriched_regions.txt

$(RESULTS)/background_regions_read_count.txt: $(RESULTS)/EWS.enriched.bins.txt
	bin/extract_genic_read_count.py < $(RESULTS)/EWS.enriched.bins.txt > $(RESULTS)/background_regions_read_count.txt

$(RESULTS)/enriched_regions_read_count.txt: $(RESULTS)/EWS.filtered.bins.txt
	bin/extract_genic_read_count.py < $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/enriched_regions_read_count.txt

$(RESULTS)/EWS.filtered.gene_bin_counts.txt: $(RESULTS)/EWS.filtered.bins.txt
	cat  $(RESULTS)/EWS.filtered.bins.txt | awk '{print $$12}' | sort | uniq -c | grep -v 'gene_used' | sort -n >  $(RESULTS)/EWS.filtered.gene_bin_counts.txt

$(RESULTS)/EWS.filtered.blocks.bed: $(RESULTS)/EWS.filtered.bins.txt
	bin/merge_bins.py < $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/EWS.filtered.blocks.bed

$(RESULTS)/EWS.CODING.blocks.bed: $(RESULTS)/EWS.filtered.blocks.bed
	# Split BED file by region
	bin/bed_region_split.py $(RESULTS)/EWS.filtered.blocks.bed $(RESULTS)/EWS.%.blocks.bed

$(RESULTS)/EWS.filtered.block_size_multiple_gene.txt: $(RESULTS)/EWS.filtered.blocks.bed
	pv $(RESULTS)/EWS.filtered.blocks.bed | bin/calc_block_size_multiple_gene.py > $(RESULTS)/EWS.filtered.block_size_multiple_gene.txt

$(RESULTS)/EWS.%.blocks.fa: $(RESULTS)/EWS.%.blocks.bed
	bedutils extend -5 50 -3 50 $< > $<.extended.bed 
	bedutils tofasta $<.extended.bed $(HG19) > $@

$(RESULTS)/EWS.NONINTRON.blocks.fa: $(RESULTS)/EWS.CODING.blocks.fa $(RESULTS)/EWS.JUNCTION.blocks.fa $(RESULTS)/EWS.UTR3.blocks.fa $(RESULTS)/EWS.UTR5.blocks.fa

	cat $< > $@

$(RESULTS)/genic_bins_by_sample.txt: $(RESULTS)/EWS.enriched.bins.txt 
	pv $(RESULTS)/EWS.enriched.bins.txt| bin/extract_genic_by_sample.py > $(RESULTS)/genic_bins_by_sample.txt


figures/clip-gene-count.pdf: bin/clip_gene_count.R $(RESULTS)/EWS.genecounts.txt
	R --vanilla < bin/clip_gene_count.R

$(RESULTS)/EWS.filtered.block_size.txt: $(RESULTS)/EWS.filtered.blocks.bed
	bin/block_size.py $(RESULTS)/EWS.filtered.blocks.bed > $(RESULTS)/EWS.filtered.block_size.txt

$(RESULTS)/EWS.gene.bincount.txt: $(RESULTS)/EWS.filtered.bins.txt
	bin/calc_bins_per_gene.py < $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/EWS.gene.bincount.txt

$(RESULTS)/EWS.enriched.reads_per_bin.txt: $(RESULTS)/EWS.filtered.bins.txt
	bin/calc_reads_per_bin.py < $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/EWS.enriched.reads_per_bin.txt
#	bin/calc_reads_per_bin.py < $(RESULTS)/EWS.background.bins.txt > $(RESULTS)/EWS.background.reads_per_bin.txt
#	bin/calc_reads_per_bin.py < $(RESULTS)/EWS.enriched.bins.txt > $(RESULTS)/EWS.all.reads_per_bin.txt

data/refseq.nm.introns.bed: $(GTF)
	ngsutilsj gtf-export --introns $(GTF) > data/refseq.nm.introns.bed

data/refseq.nm.exons.bed: $(GTF)
	ngsutilsj gtf-export --exons $(GTF) > data/refseq.nm.exons.bed

data/refseq.nm.utr3.bed: $(GTF)
	ngsutilsj gtf-export --utr3 $(GTF) > data/refseq.nm.utr3.bed

# $(RESULTS)/EWS.background.merged-da.txt: $(RESULTS)/EWS.background.bins.bed data/refseq.nm.introns.bed data/refseq.nm.exons.bed bin/intron_donor_acceptor.py bin/exon_donor_acceptor.py
# 	bin/intron_donor_acceptor.py data/refseq.nm.introns.bed results-refseq/EWS.background.bins.bed > $(RESULTS)/EWS.background.intron-da.txt
# 	bin/exon_donor_acceptor.py data/refseq.nm.exons.bed results-refseq/EWS.background.bins.bed > $(RESULTS)/EWS.background.exon-da.txt
# 	cat $(RESULTS)/EWS.background.intron-da.txt $(RESULTS)/EWS.background.exon-da.txt > $(RESULTS)/EWS.background.merged-da.txt

# $(RESULTS)/EWS.filtered.merged-da.txt: $(RESULTS)/EWS.filtered.blocks.bed data/refseq.nm.introns.bed data/refseq.nm.exons.bed bin/intron_donor_acceptor.py bin/exon_donor_acceptor.py
# 	bin/intron_donor_acceptor.py data/refseq.nm.introns.bed results-refseq/EWS.filtered.blocks.bed > $(RESULTS)/EWS.filtered.intron-da.txt
# 	bin/exon_donor_acceptor.py data/refseq.nm.exons.bed results-refseq/EWS.filtered.blocks.bed > $(RESULTS)/EWS.filtered.exon-da.txt
# 	cat $(RESULTS)/EWS.filtered.intron-da.txt $(RESULTS)/EWS.filtered.exon-da.txt > $(RESULTS)/EWS.filtered.merged-da.txt

# $(RESULTS)/EWS.filtered.confirmation.20.p3.txt: confirmation-genelist.txt $(RESULTS)/EWS.filtered.bins.bed
# 	bin/bed_filter.py confirmation-genelist.txt $(RESULTS)/EWS.filtered.bins.bed > $(TMP)/EWS.filtered.confirmation.full.bed
# 	ngsutilsj bed-reduce --extend 250 $(TMP)/EWS.filtered.confirmation.full.bed > $(RESULTS)/EWS.filtered.confirmation.bed
# 	bin/bed_prioritize.py 20 $(RESULTS)/EWS.filtered.confirmation.bed | ngsutilsj bed-resize -5 250 -3 250 --max 500 - >  $(RESULTS)/EWS.filtered.confirmation.20.bed
# 	ngsutilsj bed-tofasta $(RESULTS)/EWS.filtered.confirmation.20.bed $(HG19) > $(RESULTS)/EWS.filtered.confirmation.20.fa
# 	bin/fasta_to_primer3.py <  $(RESULTS)/EWS.filtered.confirmation.20.fa >  $(RESULTS)/EWS.filtered.confirmation.20.p3
# 	bin/primer3_core -p3_settings_file=primer3.qpcr.config -format_output $(RESULTS)/EWS.filtered.confirmation.20.p3 > $(RESULTS)/EWS.filtered.confirmation.20.p3.out.txt
# 	bin/primer3_core -p3_settings_file=primer3.qpcr.config $(RESULTS)/EWS.filtered.confirmation.20.p3 > $(RESULTS)/EWS.filtered.confirmation.20.p3.out
# 	bin/primer3_to_fasta.py < $(RESULTS)/EWS.filtered.confirmation.20.p3.out > $(RESULTS)/EWS.filtered.confirmation.20.p3.fa
# 	bedutils fromprimers -db hg19 -fasta $(RESULTS)/EWS.filtered.confirmation.20.p3.fa > $(RESULTS)/EWS.filtered.confirmation.20.p3.primers.bed
# 	bin/primer3_to_tab.py < $(RESULTS)/EWS.filtered.confirmation.20.p3.out > $(RESULTS)/EWS.filtered.confirmation.20.p3.txt


$(READS)/EWS.filtered.bins.bam $(READS)/EWS.failed.bins.noM.bam: data/EWS1.filtered.bam data/EWS2.filtered.bam $(RESULTS)/EWS.filtered.bins.bed
	mkdir -p $(READS)
	ngsutilsj bam-filter --bed-include $(RESULTS)/EWS.filtered.bins.bed --failed $(READS)/EWS1.failed.bins.bam data/EWS1.filtered.bam $(READS)/EWS1.filtered.bins.bam
	ngsutilsj bam-filter --bed-include $(RESULTS)/EWS.filtered.bins.bed --failed $(READS)/EWS2.failed.bins.bam data/EWS2.filtered.bam $(READS)/EWS2.filtered.bins.bam

	samtools merge $(READS)/EWS.filtered.bins.bam $(READS)/EWS1.filtered.bins.bam $(READS)/EWS2.filtered.bins.bam
	samtools merge $(READS)/EWS.failed.bins.bam $(READS)/EWS1.failed.bins.bam $(READS)/EWS2.failed.bins.bam

	ngsutilsj bam-filter --ref-exclude chrM $(READS)/EWS.failed.bins.bam $(READS)/EWS.failed.bins.noM.bam 

$(READS)/readlists/EWS.failed.reads.1.txt $(READS)/readlists/EWS.failed.reads.1.txt.bed: $(READS)/EWS.filtered.bins.bam $(READS)/EWS.failed.bins.noM.bam
	mkdir -p $(READS)/readlists
	ngsutilsj bam-sample --reads $$(ngsutilsj bam-stats $(READS)/EWS.filtered.bins.bam | head -n 1 | awk '{print $$2}') --lists $(JUNC_PERMUTATIONS) --output $(READS)/readlists/EWS.failed.reads $(READS)/EWS.failed.bins.noM.bam
	for f in $$(find $(READS)/readlists -name '*txt'); do ngsutilsj bam-tobed --whitelist $$f $(READS)/EWS.failed.bins.bam > $$f.bed; done

$(READS)/EWS.filtered.reads.bed: $(READS)/EWS.filtered.bins.bam
	ngsutilsj bam-tobed $(READS)/EWS.filtered.bins.bam > $(READS)/EWS.filtered.reads.bed

$(READS)/junc.donors.bed:
	ngsutilsj gtf-export --donors $(GTF) > $(READS)/junc.donors.bed

$(READS)/junc.acceptors.bed:
	ngsutilsj gtf-export --acceptors $(GTF) > $(READS)/junc.acceptors.bed

#$(READS)/tlss.bed:
#	ngsutilsj gtf-export --tlss $(GTF) > $(READS)/tlss.bed

$(READS)/utr3.bed:
	ngsutilsj gtf-export --utr3 $(GTF) > $(READS)/utr3.bed

$(READS)/genes.bed:
	ngsutilsj gtf-export --genes $(GTF) > $(READS)/genes.bed

$(READS)/readlists/EWS.failed.reads.1.txt.bed.dist: bin/donor_acceptor_dist.py $(READS)/readlists/EWS.failed.reads.1.txt.bed $(READS)/EWS.filtered.bins.bam $(READS)/junc.donors.bed $(READS)/junc.acceptors.bed $(READS)/utr3.bed
	for f in $$(find $(READS)/readlists -name '*bed'); do echo $$f; bin/donor_acceptor_dist.py $(READS)/junc.donors.bed $(READS)/junc.acceptors.bed $(READS)/utr3.bed $$f > $$f.dist; done

$(READS)/EWS.filtered.reads.dist.txt: bin/donor_acceptor_dist.py $(READS)/EWS.filtered.reads.bed $(READS)/junc.donors.bed $(READS)/junc.acceptors.bed $(READS)/utr3.bed
	bin/donor_acceptor_dist.py $(READS)/junc.donors.bed $(READS)/junc.acceptors.bed $(READS)/utr3.bed $(READS)/EWS.filtered.reads.bed > $(READS)/EWS.filtered.reads.dist.txt

$(READS)/EWS.filtered.reads.tlss-dist.txt: bin/tlss_dist.py $(READS)/EWS.filtered.reads.bed $(READS)/utr3.bed $(READS)/genes.bed 
	bin/tlss_dist.py $(READS)/utr3.bed $(READS)/genes.bed $(READS)/EWS.filtered.reads.bed > $(READS)/EWS.filtered.reads.tlss-dist.txt

$(READS)/readlists/EWS.failed.reads.1.txt.bed.tlss-dist: bin/tlss_dist.py $(READS)/readlists/EWS.failed.reads.1.txt.bed $(READS)/EWS.filtered.bins.bam $(READS)/utr3.bed $(READS)/genes.bed
	for f in $$(find $(READS)/readlists -name '*bed'); do echo $$f; bin/tlss_dist.py $(READS)/utr3.bed $(READS)/genes.bed $$f > $$f.tlss-dist; done

$(READS)/EWS.filtered.reads.splice-dist.txt: bin/splicing_dist.py $(READS)/EWS.filtered.reads.bed EWSR1-star-combined.events.txt.junctions.bed $(READS)/genes.bed 
	bin/splicing_dist.py EWSR1-star-combined.events.txt.junctions.bed $(READS)/genes.bed $(READS)/EWS.filtered.reads.bed > $(READS)/EWS.filtered.reads.splice-dist.txt

$(READS)/readlists/EWS.failed.reads.1.txt.bed.splice-dist: bin/splicing_dist.py $(READS)/readlists/EWS.failed.reads.1.txt.bed $(READS)/EWS.filtered.bins.bam EWSR1-star-combined.events.txt.junctions.bed $(READS)/genes.bed
	for f in $$(find $(READS)/readlists -name '*bed'); do echo $$f; bin/splicing_dist.py EWSR1-star-combined.events.txt.junctions.bed $(READS)/genes.bed $$f > $$f.splice-dist; done

figures/splice-dist.pdf: bin/splice_plot.R $(READS)/EWS.filtered.reads.splice-dist.txt $(READS)/readlists/EWS.failed.reads.1.txt.bed.splice-dist 
	R --vanilla < bin/splice_plot.R

figures/tlss-dist.pdf: bin/tlss_plot.R $(READS)/EWS.filtered.reads.tlss-dist.txt $(READS)/readlists/EWS.failed.reads.1.txt.bed.tlss-dist 
	R --vanilla < bin/tlss_plot.R

figures/donor-dist.pdf figures/acceptor-dist.pdf: bin/donor_acceptor_plot.Rsh $(READS)/EWS.filtered.reads.dist.txt $(READS)/readlists/EWS.failed.reads.1.txt.bed.dist 
	bin/donor_acceptor_plot.Rsh donor
	bin/donor_acceptor_plot.Rsh acceptor

$(READS)/EWS.filtered.reads.donor.reads.txt $(READS)/EWS.filtered.reads.donor.genes.txt $(READS)/EWS.filtered.reads.donor.annotated.bed: $(READS)/EWS.filtered.reads.dist.txt $(READS)/EWS.filtered.reads.bed bin/donor_acceptor_dist_filter.py
	bin/donor_acceptor_dist_filter.py DONOR -500 2000 $(READS)/EWS.filtered.reads.dist.txt $(READS)/EWS.filtered.reads.bed > $(READS)/EWS.filtered.reads.donor.annotated.bed
	ngsutilsj annotate-gtf --gtf $(GTF) --noheader --bed6 $(READS)/EWS.filtered.reads.donor.annotated.bed --gene-name | grep -v ',' > $(READS)/EWS.filtered.reads.donor.annotated.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.donor.annotated.txt | awk '{print $$8}' | sort | uniq > $(READS)/EWS.filtered.reads.donor.genes.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.donor.annotated.txt | awk '{print $$4}' > $(READS)/EWS.filtered.reads.donor.reads.txt

$(READS)/EWS.filtered.reads.acceptor.reads.txt $(READS)/EWS.filtered.reads.acceptor.genes.txt $(READS)/EWS.filtered.reads.acceptor.annotated.bed: $(READS)/EWS.filtered.reads.dist.txt $(READS)/EWS.filtered.reads.bed bin/donor_acceptor_dist_filter.py
	bin/donor_acceptor_dist_filter.py ACCEPTOR -500 2000 $(READS)/EWS.filtered.reads.dist.txt $(READS)/EWS.filtered.reads.bed > $(READS)/EWS.filtered.reads.acceptor.annotated.bed
	ngsutilsj annotate-gtf --gtf $(GTF) --noheader --bed6  $(READS)/EWS.filtered.reads.acceptor.annotated.bed --gene-name | grep -v ',' > $(READS)/EWS.filtered.reads.acceptor.annotated.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.acceptor.annotated.txt | awk '{print $$8}' | sort | uniq > $(READS)/EWS.filtered.reads.acceptor.genes.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.acceptor.annotated.txt | awk '{print $$4}' > $(READS)/EWS.filtered.reads.acceptor.reads.txt

$(READS)/EWS.filtered.reads.tlss.reads.txt $(READS)/EWS.filtered.reads.tlss.genes.txt $(READS)/EWS.filtered.reads.tlss.annotated.bed: $(READS)/EWS.filtered.reads.tlss-dist.txt $(READS)/EWS.filtered.reads.bed  bin/donor_acceptor_dist_filter.py
	bin/donor_acceptor_dist_filter.py TLSS -500 2000 $(READS)/EWS.filtered.reads.tlss-dist.txt $(READS)/EWS.filtered.reads.bed > $(READS)/EWS.filtered.reads.tlss.annotated.bed
	ngsutilsj annotate-gtf --gtf $(GTF) --noheader --bed6 $(READS)/EWS.filtered.reads.tlss.annotated.bed --gene-name | grep -v ',' > $(READS)/EWS.filtered.reads.tlss.annotated.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.tlss.annotated.txt | awk '{print $$8}' | sort | uniq > $(READS)/EWS.filtered.reads.tlss.genes.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.tlss.annotated.txt | awk '{print $$4}' > $(READS)/EWS.filtered.reads.tlss.reads.txt

$(READS)/EWS.filtered.reads.donor.reads.fa: $(READS)/EWS.filtered.reads.donor.reads.txt bin/reads_to_meme_seq.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.donor.reads.txt $(HG19) 50 > $(READS)/EWS.filtered.reads.donor.reads.fa
	samtools view data/EWS2.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.donor.reads.txt $(HG19) 50 >> $(READS)/EWS.filtered.reads.donor.reads.fa

$(READS)/EWS.filtered.reads.acceptor.reads.fa: $(READS)/EWS.filtered.reads.acceptor.reads.txt bin/reads_to_meme_seq.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.acceptor.reads.txt $(HG19) 50 > $(READS)/EWS.filtered.reads.acceptor.reads.fa
	samtools view data/EWS2.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.acceptor.reads.txt $(HG19) 50 >> $(READS)/EWS.filtered.reads.acceptor.reads.fa

$(READS)/EWS.filtered.reads.tlss.reads.fa: $(READS)/EWS.filtered.reads.tlss.reads.txt bin/reads_to_meme_seq.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.tlss.reads.txt $(HG19) 50 > $(READS)/EWS.filtered.reads.tlss.reads.fa
	samtools view data/EWS2.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.tlss.reads.txt $(HG19) 50 >> $(READS)/EWS.filtered.reads.tlss.reads.fa

$(READS)/EWS.filtered.reads.bg.reads.txt: $(READS)/EWS.filtered.reads.bed
	grep -v '^#' $(READS)/EWS.filtered.reads.bed | awk '{print $$4}' > $(READS)/EWS.filtered.reads.bg.reads.txt

$(READS)/EWS.filtered.reads.genes.txt: $(READS)/EWS.filtered.reads.bed
	ngsutilsj annotate-gtf --gtf $(GTF) --noheader --bed6  $(READS)/EWS.filtered.reads.bed --gene-name | grep -v ',' > $(READS)/EWS.filtered.reads.annotated.txt
	grep -v '^#' $(READS)/EWS.filtered.reads.annotated.txt | awk '{print $$8}' | sort | uniq > $(READS)/EWS.filtered.reads.genes.txt

$(READS)/EWS.filtered.reads.bg.reads.fa: $(READS)/EWS.filtered.reads.bg.reads.txt bin/reads_to_meme_seq.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.bg.reads.txt $(HG19) 50 bg > $(READS)/EWS.filtered.reads.bg.reads.fa.tmp
	samtools view data/EWS2.filtered.bam | bin/reads_to_meme_seq.py $(READS)/EWS.filtered.reads.bg.reads.txt $(HG19) 50 bg >> $(READS)/EWS.filtered.reads.bg.reads.fa.tmp
	mv $(READS)/EWS.filtered.reads.bg.reads.fa.tmp $(READS)/EWS.filtered.reads.bg.reads.fa

$(READS)/EWS.filtered.reads.donor.gapcount.txt: $(READS)/EWS.filtered.reads.donor.reads.txt bin/reads_to_gapcount.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_gapcount.py $(READS)/EWS.filtered.reads.donor.reads.txt > $(READS)/EWS.filtered.reads.donor.gapcount.txt
	samtools view data/EWS2.filtered.bam | bin/reads_to_gapcount.py $(READS)/EWS.filtered.reads.donor.reads.txt >> $(READS)/EWS.filtered.reads.donor.gapcount.txt

$(READS)/EWS.filtered.reads.acceptor.gapcount.txt: $(READS)/EWS.filtered.reads.acceptor.reads.txt bin/reads_to_gapcount.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_gapcount.py $(READS)/EWS.filtered.reads.acceptor.reads.txt > $(READS)/EWS.filtered.reads.acceptor.gapcount.txt
	samtools view data/EWS2.filtered.bam | bin/reads_to_gapcount.py $(READS)/EWS.filtered.reads.acceptor.reads.txt >> $(READS)/EWS.filtered.reads.acceptor.gapcount.txt

$(READS)/EWS.filtered.reads.tlss.gapcount.txt: $(READS)/EWS.filtered.reads.tlss.reads.txt bin/reads_to_gapcount.py
	samtools view data/EWS1.filtered.bam | bin/reads_to_gapcount.py $(READS)/EWS.filtered.reads.tlss.reads.txt > $(READS)/EWS.filtered.reads.tlss.gapcount.txt
	samtools view data/EWS2.filtered.bam | bin/reads_to_gapcount.py $(READS)/EWS.filtered.reads.tlss.reads.txt >> $(READS)/EWS.filtered.reads.tlss.gapcount.txt


$(READS)/gene-donor-acceptor-tlss-overlap.svg: $(READS)/EWS.filtered.reads.donor.genes.txt $(READS)/EWS.filtered.reads.acceptor.genes.txt $(READS)/EWS.filtered.reads.tlss.genes.txt
	txtutils venn --names "Donor,Acceptor,3'UTR" -o $(READS)/gene-donor-acceptor-tlss-overlap.svg --svg --title "Donor/Acceptor/3'UTR overlap" $(READS)/EWS.filtered.reads.donor.genes.txt $(READS)/EWS.filtered.reads.acceptor.genes.txt $(READS)/EWS.filtered.reads.tlss.genes.txt

# $(READS)/EWS.filtered.reads.acceptor.fa $(READS)/EWS.filtered.reads.donor.fa $(READS)/EWS.filtered.reads.tlss.fa: $(READS)/EWS.filtered.reads.acceptor.annotated.bed $(READS)/EWS.filtered.reads.donor.annotated.bed $(READS)/EWS.filtered.reads.tlss.annotated.bed
# 	ngsutilsj bed-resize -3 100 -5 100 $(READS)/EWS.filtered.reads.acceptor.annotated.bed | ngsutilsj bed-tofasta - $(HG19) > $(READS)/EWS.filtered.reads.acceptor.fa
# 	ngsutilsj bed-resize -3 100 -5 100 $(READS)/EWS.filtered.reads.donor.annotated.bed | ngsutilsj bed-tofasta - $(HG19) > $(READS)/EWS.filtered.reads.donor.fa
# 	ngsutilsj bed-resize -3 100 -5 100 $(READS)/EWS.filtered.reads.tlss.annotated.bed | ngsutilsj bed-tofasta - $(HG19) > $(READS)/EWS.filtered.reads.tlss.fa

# $(READS)/EWS.filtered.reads.UTR3.bed $(READS)/EWS.filtered.reads.DONOR.bed $(READS)/EWS.filtered.reads.ACCEPTOR.bed: $(READS)/EWS.filtered.reads.bed $(READS)/EWS.filtered.reads.dist.txt
# 	bin/merge_reads_to_bed.py utr3 $(READS)/EWS.filtered.reads.bed $(READS)/EWS.filtered.reads.dist.txt > tmp.bed
# 	ngsutilsj bed-resize -3 100 -5 100 tmp.bed > tmp1.bed
# 	ngsutilsj bed-reduce --rename --extend 100 tmp1.bed > $(READS)/EWS.filtered.reads.UTR3.bed

# 	bin/merge_reads_to_bed.py donor $(READS)/EWS.filtered.reads.bed $(READS)/EWS.filtered.reads.dist.txt > tmp.bed
# 	ngsutilsj bed-resize -3 100 -5 100 tmp.bed > tmp1.bed
# 	ngsutilsj bed-reduce --rename --extend 100 tmp1.bed > $(READS)/EWS.filtered.reads.DONOR.bed

# 	bin/merge_reads_to_bed.py acceptor $(READS)/EWS.filtered.reads.bed $(READS)/EWS.filtered.reads.dist.txt > tmp.bed
# 	ngsutilsj bed-resize -3 100 -5 100 tmp.bed > tmp1.bed
# 	ngsutilsj bed-reduce --rename --extend 100 tmp1.bed > $(READS)/EWS.filtered.reads.ACCEPTOR.bed

# 	rm tmp1.bed tmp.bed

# $(READS)/EWS.filtered.reads.UTR3.fa $(READS)/EWS.filtered.reads.DONOR.fa $(READS)/EWS.filtered.reads.ACCEPTOR.fa: $(READS)/EWS.filtered.reads.UTR3.bed $(READS)/EWS.filtered.reads.DONOR.bed $(READS)/EWS.filtered.reads.ACCEPTOR.bed
# 	ngsutilsj bed-tofasta $(READS)/EWS.filtered.reads.UTR3.bed $(HG19) > $(READS)/EWS.filtered.reads.UTR3.fa
# 	ngsutilsj bed-tofasta $(READS)/EWS.filtered.reads.DONOR.bed $(HG19) > $(READS)/EWS.filtered.reads.DONOR.fa
# 	ngsutilsj bed-tofasta $(READS)/EWS.filtered.reads.ACCEPTOR.bed $(HG19) > $(READS)/EWS.filtered.reads.ACCEPTOR.fa


$(RESULTS)/EWS.genelist.UTR3.txt: $(RESULTS)/EWS.filtered.bins.txt data/refseq.nm.utr3.bed
	bin/bin_bed_filter.py $(RESULTS)/EWS.filtered.bins.txt data/refseq.nm.utr3.bed | awk '{print $$12}' | sort | uniq | grep -v 'gene_used' > $(RESULTS)/EWS.genelist.UTR3.txt

$(RESULTS)/EWS.genelist.txt: $(RESULTS)/EWS.filtered.bins.txt
	bin/extract_gene_list.py $(RESULTS)/EWS.filtered.bins.txt > $(RESULTS)/EWS.genelist.txt

# $(RESULTS)/EWS.filtered.donors-500.txt: $(RESULTS)/EWS.filtered.merged-da.txt bin/da_donor_filter.py
# 	bin/da_donor_filter.py $(RESULTS)/EWS.filtered.merged-da.txt > $(RESULTS)/EWS.filtered.donors-500.txt

# $(RESULTS)/EWS.genelist.donors.txt: $(RESULTS)/EWS.filtered.donors-500.txt
# 	cat results-refseq/EWS.filtered.donors-500.txt | sed -e 's/;/ /' | awk '{print $$4}' | sort | uniq > $(RESULTS)/EWS.genelist.donors.txt 

# $(RESULTS)/EWS.filtered.UTR3.bed: $(RESULTS)/EWS.filtered.bins.txt data/refseq.nm.utr3.bed
# 	bin/bin_bed_filter.py $(RESULTS)/EWS.filtered.bins.txt data/refseq.nm.utr3.bed | bin/bin_to_bed.py - | ngsutilsj bed-resize -5 100 -3 100 -  | ngsutilsj bed-reduce - > $(RESULTS)/EWS.filtered.UTR3.bed

# $(RESULTS)/EWS.filtered.UTR3.fa: $(RESULTS)/EWS.filtered.UTR3.bed
# 	ngsutilsj bed-tofasta $(RESULTS)/EWS.filtered.UTR3.bed $(HG19) > $(RESULTS)/EWS.filtered.UTR3.fa

# $(RESULTS)/EWS.filtered.donors-500.bed: $(RESULTS)/EWS.filtered.bins.txt data/refseq.nm.utr3.bed
# 	bin/da_donor_filter.py -bed $(RESULTS)/EWS.filtered.merged-da.txt | ngsutilsj bed-resize -5 100 -3 100 -  | ngsutilsj bed-reduce - > $(RESULTS)/EWS.filtered.donors-500.bed

# $(RESULTS)/EWS.filtered.donors-500.fa: $(RESULTS)/EWS.filtered.donors-500.bed
# 	ngsutilsj bed-tofasta $(RESULTS)/EWS.filtered.donors-500.bed $(HG19) > $(RESULTS)/EWS.filtered.donors-500.fa

# $(RESULTS)/EWS.filtered.donors-500.regions.txt: $(RESULTS)/EWS.filtered.donors-500.bed
# 	bin/bed_to_region.py < $(RESULTS)/EWS.filtered.donors-500.bed > $(RESULTS)/EWS.filtered.donors-500.regions.txt

# $(RESULTS)/EWS.filtered.UTR3.regions.txt: $(RESULTS)/EWS.filtered.UTR3.bed
# 	bin/bed_to_region.py < $(RESULTS)/EWS.filtered.UTR3.bed > $(RESULTS)/EWS.filtered.UTR3.regions.txt

$(RESULTS)/EWS.genelist.bg.txt: bin/bg_genes.py $(RESULTS)/EWS.enriched.bins.txt
	bin/bg_genes.py $(RESULTS)/EWS.enriched.bins.txt > $(RESULTS)/clip_present_genes.txt
