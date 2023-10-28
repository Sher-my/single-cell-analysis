#!/bin/bash
###This sctipt simply deals with single cell data
wkd=/home/lf-3/yt/rawdatas/sc_10/
sample=sample_name
###Get ref index from https://www.10xgenomics.com/support/software/space-ranger/downloads or made by yourself:
cellranger mkref --genome=${wkd}/ref/mm10_cellranger --nthreads=10 --fasta= ${wkd}/ref/mm10.fa --genes=${wkd}/ref/mm10.refGene.gtf
###Matric
cd ${wkd}
for s in ${sample}
do
	cellranger count --id=${s}_cellranger_count --fastqs=${wkd} --sample=${s} --transcriptome=${wkd}/ref/refdata-gex-mm10-2020-A --r1-length 26
done
