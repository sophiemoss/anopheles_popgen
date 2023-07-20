## download additional samples from ncbi using fasterq-dump
conda activate fasterq-dump
fasterq-dump SRR000001

## to make a genomics database of sample VCFs, use the following

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import --sample-file fastq2vcfsamples.txt --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa --prefix gambiae_nov2022 --vcf-dir .

## now merge VCF files

merge_vcfs.py genotype --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa --prefix gambiae_nov2022 > mergevcf_log.txt 2>&1

# resulting vcf is called gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# concerned that this genotyping does not use a set of validated variants like I had for plasmodium
# this has a genotyping pipeline:
# https://github.com/malariagen/pipelines/blob/v0.0.4/docs/specs/snp-genotyping-vector.md
# show Jody pipeline from malariagen