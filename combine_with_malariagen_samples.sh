########## MALARIAGEN DATA - phase 3 ##################

## download vcfs from malariagen

# view species calls for malariagen samples

head ~/vo_agam_release/v3/metadata/species_calls_aim_20220528/AG1000G-GM-C/samples.species_aim.csv

awk -F ',' '$7 == "gambiae" && $7 != "intermediate_gambiae_coluzzii" {print $1, $7}' ~/vo_agam_release/v3/metadata/species_calls_aim_20220528/AG1000G-GM-C/samples.species_aim.csv

# view sample ID and snp_genotypes_vcf for Guinea-Bissau dataset

head ~/vo_agam_release/v3/metadata/general/AG1000G-GW/wgs_snp_data.csv | cut -f1,4 -d,

# made file with all sample IDs and now downloading using python script
# python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/download_malariagen_samples.py

# merge the malariagen VCF files
bcftools merge --output-type z --output malariagen_gambiae_GB_GM-ABC_merged.vcf.gz *.vcf.gz

# I got these errors, discuss with someone:

bcftools merge --output-type z --output malariagen
_gambiae_GB_GM-ABC_merged.vcf.gz *.vcf.gz
[W::bcf_sr_add_reader] No BGZF EOF marker; file 'AJ0023-C.vcf.gz' may be truncated
[W::hts_idx_load3] The index file is older than the data file: AJ0023-C.vcf.gz.tbi
[W::hts_idx_load3] The index file is older than the data file: AJ0023-C.vcf.gz.tbi
[W::vcf_parse] INFO 'Excess' is not defined in the header, assuming Type=String
Error: The INFO field is not defined in the header: Excess

## UP TO HERE. NEED TO COMBINE WITH MY SAMPLES.

## bcftools concatenate is used for when you want to merge 2 or more VCFs with the exact same samples, 
## but for example, each VCF is for a single chromosome. 
## bcftools merge is when you want to merge 2 or more VCFs with non-overlapping samples 

## merge my samples with malariagen samples vcf when they have been aligned to Agam_P3 and have been combined into one large VCF

bcftools merge 1.vcf.gz 2.vcf.gz -O z -o output.vcf.gz