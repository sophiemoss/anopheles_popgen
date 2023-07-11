## to make a genomics database of sample VCFs, use the following

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import --sample-file fastq2vcfsamples.txt --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa --prefix gambiae_nov2022 --vcf-dir .

## now merge VCF files

merge_vcfs.py genotype --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa --prefix gambiae_nov2022 > mergevcf_log.txt 2>&1

# resulting vcf should be called gambiae_nov2022.2023_07_05.genotyped.vcf.gz