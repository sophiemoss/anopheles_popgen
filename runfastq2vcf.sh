#!/bin/bash
#the reference file needs to be in the same directory!
# 2>&1 means "copy standard errors to standard output", and the > redicrects the standard output (now containing erros) to the file log.txt

ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > fastq2vcfsamples.txt

cat fastq2vcfsamples.txt | parallel -j 60 \
"fastq2vcf.py all --read1 {}_1.fastq.gz --read2 {}_2.fastq.gz \
--ref /mnt/storage11/sophie/bijagos_mosq_wgs/fastq2vcf_out_altogambiae_combinedreads/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
--prefix {}" > fastq2vcf_log.txt 2>&1
