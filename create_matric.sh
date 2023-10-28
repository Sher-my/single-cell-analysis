#!/bin/bash
###This sctipt simply deals with single cell data
wkd=/home/lf-3/yt/rawdatas/sc_10/
cd ${wkd}
###Build ref index
STAR --runMode genomeGenerate --genomeDir index --runThreadN 16 --genomeFastaFiles /home/lf-3/yt/rawdatas/sc_10/raw_genome_gtf/GRCh38.P13.genome.fa --sjdbGTFfile /home/lf-3/yt/rawdatas/sc_10/raw_genome_gtf/gencode.v38.annotation.gtf --sjdbOverhang 149
###Mapping
STAR --runThreadN 16 --genomeDir /lf-3/yt/rawdatas/sc_10/ --readFilesIn /home/lf-3/yt/rawdatas/sc_10/samplename_S1_L004_R1_001.fastq.gz /home/lf-3/yt/rawdatas/sc_10/samplename_S1_L004_R2_001.fastq.gz --readFilesCommand gunzip -c --outFileNamePrefix /home/lf-3/yt/rawdatas/sc_10/raw_out/ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile /home/lf-3/yt/rawdatas/sc_10/raw-genes.gtf
###Matric
cellranger count --id = sample1 --transcriptome = /lf-3/yt/rawdatas/sc_10/ --fastqs = /home/lf-3/yt/rawdatas/sc_10/ --sample = samplename --expect-cells = 1000 --nosecondary
