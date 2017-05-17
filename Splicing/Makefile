all: EWSR1-star-combined.events.txt EWSR1-star-combined.gene.txt events-genelist.txt events.bed valid-ensembl.txt splicing-background-genes.txt splicing-clip-overlap.txt clip-splicing-overlap-genelist.txt EF-HBM-combined.events.txt EF-events-genelist.txt donor-acceptor-utrclip-splicing.svg

clean:
	rm EWSR1-star-combined.txt EWSR1-star-combined.events.txt genelist.txt

EF-HBM-combined.txt:
	cgsplice splice-diff -v --groups 1,2,2,2,1,1,1,2,2,2 --max-edit-distance 2 --min-total-count 50 -o EF-HBM-combined.txt raw/G{3,4,5,6}*counts2a.txt.gz raw/H{1,2,3,4,5,6}*counts2a.txt.gz

EWSR1-star-combined.txt:
	cgsplice splice-diff -v --groups 1,1,1,2,2,2,2 --max-edit-distance 2 --min-total-count 50 -o EWSR1-star-combined.txt raw/Control?.counts2a.txt.gz raw/EWS?.counts2a.txt.gz

%.events.txt: %.txt
	cgsplice combine-events --failed $@.failed.bed --bed $@.junctions.bed $< > $@.tmp
	ngsutilsj annotate-gtf --col-junction 2 --col-strand 3 --within-gene --gtf ~/Research/references/hg19/ensembl.74.hg19.gtf.gz --gene-name --biotype  $@.tmp > $@
	rm $@.tmp

%.gene.txt: %.txt
	ngsutilsj annotate-gtf --col-junction 1 --col-strand 2 --within-gene --gtf ~/Research/references/hg19/ensembl.74.hg19.gtf.gz --gene-name  $< > $@

EF-events-genelist.txt: EF-HBM-combined.events.txt
	bin/event_gene_list.py < EF-HBM-combined.events.txt  | sort | uniq > EF-events-genelist.txt

events-genelist.txt: EWSR1-star-combined.events.txt
	bin/event_gene_list.py < EWSR1-star-combined.events.txt  | sort | uniq > events-genelist.txt

events.bed: EWSR1-star-combined.events.txt
	bin/event_to_bed.py  < EWSR1-star-combined.events.txt  > events.bed

valid-junctions.txt: EWSR1-star-combined.events.txt
	bin/export_valid_junctions.py < EWSR1-star-combined.events.txt  | sort | uniq > valid-junctions.txt

splicing-background-genes.txt: EWSR1-star-combined.gene.txt bin/export_bg_genes.py
	bin/export_bg_genes.py < EWSR1-star-combined.gene.txt > splicing-background-genes.txt

refseq.junctions.txt: refseq.nm.gtf.gz
	ngsutilsj gtf-export --junctions refseq.nm.gtf.gz > refseq.junctions.txt

ensembl.74.junctions.txt: ~/Research/references/hg19/ensembl.74.hg19.gtf.gz
	ngsutilsj gtf-export --junctions ~/Research/references/hg19/ensembl.74.hg19.gtf.gz > ensembl.74.junctions.txt

valid-ensembl.txt: valid-junctions.txt ensembl.74.junctions.txt
	cat valid-junctions.txt ensembl.74.junctions.txt | sort | uniq > valid-ensembl.txt

splicing-clip-overlap.txt: events-genelist.txt
	bin/splicing_clip_overlap.py events-genelist.txt ~/Research/Projects/EWS/CLIP/results-refseq/EWS.filtered.bins.bed > splicing-clip-overlap.txt

clip-splicing-overlap-genelist.txt: events-genelist.txt
	txtutils overlap ~/Research/Projects/EWS/CLIP/results-refseq/EWS.genelist.txt  events-genelist.txt > clip-splicing-overlap-genelist.txt

Venn-overlap.svg: events-genelist.txt
	txtutils venn --svg -o Venn-overlap.svg --names CLIP,Splicing,Exp-Up,Exp-Down ~/Research/Projects/EWS/CLIP/results-refseq/EWS.genelist.txt  events-genelist.txt ../RNASeq\ expression/up-genes.txt ../RNASeq\ expression/down-genes.txt

donor-acceptor-utrclip-splicing.svg: events-genelist.txt
	txtutils venn --names Spliced,CLIP-Donor,CLIP-Acceptor,CLIP-3\'UTR --svg -o donor-acceptor-utrclip-splicing.svg events-genelist.txt ~/Research/Projects/EWS/CLIP/results-bam/EWS.filtered.reads.donor.genes.txt ~/Research/Projects/EWS/CLIP/results-bam/EWS.filtered.reads.acceptor.genes.txt ~/Research/Projects/EWS/CLIP/results-bam/EWS.filtered.reads.tlss.genes.txt
