#!/bin/bash
genomedir=/disk/yt/single_cell/ref
datadir=/disk/yt/single_cell/01matrix
sample=sample_name
###step 1 download 10× ref from https://www.10xgenomics.com/support/software/space-ranger/downloads or made it by yourself
cellranger mkref --genome=${genomedir}/mm10_cellranger --nthreads=10 --fasta=${genomedir}/mm10.fa --genes=${genomedir}/mm10.refGene.gtf
###step2 create matrix
cd ${dataidr}
for s in ${sample}
do
	cellranger count --id=${s}_cellranger_count --fastqs=${datadir} --sample=${s} --transcriptome=${genomedir}/refdata-gex-mm10-2020-A --r1-length 26
done
