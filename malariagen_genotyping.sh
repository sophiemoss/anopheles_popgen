# Genotyping BAM files using malariagen 'Mosquito SNP genotyping pipeline specification'
# https://github.com/malariagen/pipelines/blob/991b0328ea8027bc6b1137f893a5340e27c8e87c/docs/specs/snp-genotyping-vector.md
# to use this pipeline you need bam files which have been aligned to the reference genome
# malariagen have a pipeline for this (short-read-alignment-vector.md), however mine are already aligned and I have the bam files.

# I downloaded the variant file from ensembl http://ftp.ensemblgenomes.org/pub/metazoa/release-57/variation/vcf/anopheles_gambiae/anopheles_gambiae.vcf.gz
# check but I think this is probably the correct file - downloaded this (anopheles_gambiae.vcf.gz) and renamed this ensembl_agamP4_germlinevariations_agamp4.vcf.gz

java -jar GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I {sample BAM} \
    --alleles {alleles VCF} \
    -R {reference sequence} \
    --out {output VCF} \
    --genotype_likelihoods_model BOTH \
    --genotyping_mode GENOTYPE_GIVEN_ALLELES \ 
    --heterozygosity 0.015 \
    --heterozygosity_stdev 0.05 \
    --indel_heterozygosity 0.001 \
    --downsampling_type BY_SAMPLE \
    -dcov 250 \
    --output_mode EMIT_ALL_SITES \
    --min_base_quality_score 17 \
    -stand_call_conf 0.0 \
    -contamination 0.0 \
    -A DepthPerAlleleBySample \
    -XA RMSMappingQuality \
    -XA Coverage \
    -XA ExcessHet \
    -XA InbreedingCoeff \
    -XA MappingQualityZero \
    -XA HaplotypeScore \
    -XA SpanningDeletions \
    -XA FisherStrand \
    -XA StrandOddsRatio \
    -XA ChromosomeCounts \
    -XA BaseQualityRankSumTest \
    -XA MappingQualityRankSumTest \
    -XA QualByDepth \
    -XA ReadPosRankSumTest

    ## 

java -jar GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I bu1002_Combined.mkdup.bam \
    --alleles /mnt/storage11/sophie/bijagos_mosq_wgs/ensembl_agamP4_germlinevariations_agamp4.vcf.gz \
    -R Anopheles_gambiae.AgamP4.dna.toplevel.fa \
    --out bu1002_Combined \
    --genotype_likelihoods_model BOTH \
    --genotyping_mode GENOTYPE_GIVEN_ALLELES \ 
    --heterozygosity 0.015 \
    --heterozygosity_stdev 0.05 \
    --indel_heterozygosity 0.001 \
    --downsampling_type BY_SAMPLE \
    -dcov 250 \
    --output_mode EMIT_ALL_SITES \
    --min_base_quality_score 17 \
    -stand_call_conf 0.0 \
    -contamination 0.0 \
    -A DepthPerAlleleBySample \
    -XA RMSMappingQuality \
    -XA Coverage \
    -XA ExcessHet \
    -XA InbreedingCoeff \
    -XA MappingQualityZero \
    -XA HaplotypeScore \
    -XA SpanningDeletions \
    -XA FisherStrand \
    -XA StrandOddsRatio \
    -XA ChromosomeCounts \
    -XA BaseQualityRankSumTest \
    -XA MappingQualityRankSumTest \
    -XA QualByDepth \
    -XA ReadPosRankSumTest